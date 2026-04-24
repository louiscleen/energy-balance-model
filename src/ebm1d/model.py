"""
Core implementation of the 1D Energy Balance Model (EBM) based on Budyko's formulation.
"""

from typing import Any

import numpy as np
import numpy.typing as npt
import ebm1d.constants as cst
import ebm1d.config as cfg
import ebm1d.inputs as inp

import matplotlib.pyplot as plt

FloatArray = npt.NDArray[np.float64]


class EBM1DBudyko:
    YEAR_IN_DAYS = cst.YEAR_IN_DAYS

    def __init__(self, config: cfg.EBM1DConfig, input_data: inp.InputData, seasonal: bool = None):
        self.cfg = config
        self.input_data = input_data

        self.n_lat = config.model.n_lat
        self._seasonal = config.insolation.seasonal if seasonal is None else seasonal
        self.heat_capacity = config.heat_capacity.C_value
        self._S0 = config.insolation.S0
        self.K_meridional = config.transport.K_meridional

        time = config.time
        alb = config.albedo
        rad = config.radiation

        self.dt_days = time.dt
        self.n_years = time.n_years

        self.alpha_atmosphere = alb.alpha_atmosphere
        self.alpha_cloud = alb.alpha_cloud
        self.alpha_land = alb.alpha_land
        self.alpha_ocean = alb.alpha_ocean
        self.alpha_snow = alb.alpha_snow
        self.alpha_ice = alb.alpha_ice
        self.T_ice_min = alb.T_ice_min
        self.T_ice_max = alb.T_ice_max
        self.T_snow_min = alb.T_snow_min
        self.T_snow_max = alb.T_snow_max

        self.p_CO2 = rad.p_CO2
        self.A1 = rad.A1
        self.A2 = rad.A2
        self.B1 = rad.B1
        self.B2 = rad.B2
        self.C = rad.C
        self.D = rad.D

        print (f"Initializing the 1D EBM with {self.n_lat} latitude band and seasonal insolation = {self._seasonal}")

        # Precompute latitude-dependent parameters such as area weights and latitudinal coordinates
        lat_edges_rad = np.linspace(-np.pi / 2, np.pi / 2, self.n_lat + 1)
        self.area_weights = np.sin(lat_edges_rad[1:]) - np.sin(lat_edges_rad[:-1])
        self.area = 2 * np.pi * cst.EARTH_RADIUS**2 * self.area_weights
        self.f_area = self.area / np.sum(self.area)
        self.lat_centers_rad = (lat_edges_rad[:-1] + lat_edges_rad[1:]) / 2
        self.lat_centers_deg = np.degrees(self.lat_centers_rad)

        # Set the method to compute the solar flux based on whether we are using seasonal insolation or not
        if self._seasonal:
            self.get_solar_flux = self._compute_solar_flux_seasonal
        else:
            # For the non-seasonal case, we precompute the solar flux once since it does not depend on time, to save computation during integration.
            self._precompute_solar_flux = self._compute_solar_flux_annual_mean()
            self.get_solar_flux = lambda t: self._precompute_solar_flux

    @property
    def S0(self) -> float:
        return self._S0
        
    @S0.setter
    def S0(self, value: float):
        self._S0 = value
        if not self._seasonal:
            self._precompute_solar_flux = self._compute_solar_flux_annual_mean()

    def compute_fraction_ice(self, T: FloatArray) -> FloatArray:
        """ Compute the fraction of ice as a function of temperature, with a linear transition between T_ice_min and T_ice_max."""
        return np.clip((self.T_ice_max - T) / (self.T_ice_max - self.T_ice_min), 0.0, 1.0)

    def compute_fraction_snow(self, T: FloatArray) -> FloatArray:
        """ Compute the fraction of snow as a function of temperature, with a linear transition between T_snow_min and T_snow_max."""
        return np.clip((self.T_snow_max - T) / (self.T_snow_max - self.T_snow_min), 0.0, 1.0)
    
    def compute_albedo_surf(self, T: FloatArray) -> FloatArray:
        """ Compute the surface albedo as a function of temperature, considering the fractions of snow and ice, and the land/ocean distribution."""
        f_land = self.input_data.f_land
        self.f_ice = self.compute_fraction_ice(T)
        self.f_snow = self.compute_fraction_snow(T)

        albedo_continent = self.f_snow * self.alpha_snow + (1 - self.f_snow) * self.alpha_land
        albedo_ocean = self.f_ice * self.alpha_ice + (1 - self.f_ice) * self.alpha_ocean
        albedo_surf = f_land * albedo_continent + (1 - f_land) * albedo_ocean
        return albedo_surf

    def compute_albedo_clear_sky(self, T: FloatArray) -> FloatArray:
        """ Compute the clear sky albedo as a function of temperature, considering the surface albedo and atmospheric effects."""
        albedo_surf = self.compute_albedo_surf(T)

        # The calculation take the form of a multiple reflection series which we truncate at the second order term
        albedo_clear_sky = self.alpha_atmosphere + (1 - self.alpha_atmosphere)**2 * albedo_surf
        return albedo_clear_sky
    
    def compute_albedo_cloudy_sky(self, T: FloatArray) -> FloatArray:
        """ Compute the cloudy sky albedo as a function of temperature, considering the cloud albedo and atmospheric effects."""
        #The tops of the clouds are located higher up in the atmosphere, so we arbitrarily divide by 2.
        albedo_cloudy_sky = self.alpha_atmosphere / 2.0 + (1 - self.alpha_atmosphere / 2.0) * self.alpha_cloud
        return albedo_cloudy_sky

    def compute_albedo(self, T: FloatArray) -> FloatArray:
        """
        Calculation of albedo as a function of temperature.
        T : temperature table by latitude (in K)
        Returns a table of the same format with the albedo (dimensionless).
        """
        albedo_clear_sky = self.compute_albedo_clear_sky(T)
        albedo_cloudy_sky = self.compute_albedo_cloudy_sky(T)

        # Final albedo: weighted average between clear sky and cloudy sky
        albedo = self.input_data.f_cloud * albedo_cloudy_sky + (1 - self.input_data.f_cloud) * albedo_clear_sky
        return albedo

    def _compute_solar_flux_annual_mean(self) -> FloatArray:
        """
        Calculation of the annual mean solar flux using Legendre polynomials
        Returns an array of size n_lat with the flux in W m^-2.
        """
        # Simplified annual mean insolation: Q = S0/4 * (1 - 0.241 * (3 * sin^2(lat) - 1))
        x = np.sin(self.lat_centers_rad)
        insolation_weight = (1.0 - 0.241 * (3 * x**2 - 1))
        S = (self.S0 / 4.0) * insolation_weight
        return S

    def _compute_solar_flux_seasonal(self, t_days) -> FloatArray:
        """ Calculation of the seasonal solar flux as a function of time
            t_days : time in days
            Returns an array of size n_lat with the flux in W m^-2.
        """

        day = t_days % self.YEAR_IN_DAYS

        omega = 2 * np.pi * day / self.YEAR_IN_DAYS
        nu = omega - 2 * np.pi * cst.PERIHELION_DAY / self.YEAR_IN_DAYS
        lambda_value = omega - 2 * np.pi * cst.VERNAL_EQUINOX_DAY / self.YEAR_IN_DAYS
        
        # Calculation of the solar declination angle
        delta = np.arcsin(cst.SIN_EPSILON * np.sin(lambda_value))    
        
        # Calculation of the hour angle at sunset
        cos_H0 = - np.tan(delta) * np.tan(self.lat_centers_rad)
        H0 = np.arccos(np.clip(cos_H0, -1, 1))

        # Calculation of the distance correction factor due to Earth's eccentricity
        dist_factor = 1 + 2*cst.EARTH_ECCENTRICITY * np.cos(nu)

        daily_mean_sin_component = H0 * np.sin(self.lat_centers_rad) * np.sin(delta)
        daily_mean_cos_component = np.cos(self.lat_centers_rad) * np.cos(delta) * np.sin(H0)

        S = self.S0 * dist_factor * (daily_mean_sin_component + daily_mean_cos_component) / np.pi
        return S
    
    def compute_outgoing_IR_flux(self, T: FloatArray) -> FloatArray:
        """
        Calculation of outgoing IR radiation flux using a linear parameterization.
        Returns an array of size n_lat with the flux in W m^-2.
        """
        A = self.A1 + self.A2 * np.log(self.p_CO2 / cst.P_CO2_0)
        B = self.B1 + self.B2 * np.log(self.p_CO2 / cst.P_CO2_0)
        outgoing_IR_flux = (A + B * (T - 273.15)) - (self.C + self.D * (T - 273.15)) * self.input_data.f_cloud
        return outgoing_IR_flux

    def compute_meridional_heat_transport(self, T: FloatArray) -> FloatArray:
        """
        Calculation of meridional heat transport using Budyko's parameterization.
        Returns an array of size n_lat with the flux in W m^-2.
        """
        # Area-weighted global temperature
        T_global = self.global_temperature(T)  

        # Heat transport proportional to the temperature difference
        return self.K_meridional * (T - T_global)
    
    def compute_temperature_tendency(self, t_days, T: FloatArray) -> FloatArray: #compute_delta_T
        """
        Calculation of the temperature tendency dT/dt according to the EBM equation.
        T : array of temperatures by latitude (in K)
        t_days : time in days
        Returns an array of the same shape with the temperature tendency (in K/s).
        """
        
        S = self.get_solar_flux(t_days)
        alpha = self.compute_albedo(T)
        absorbed_solar_flux = S * (1 - alpha)
        outgoing_IR_flux = self.compute_outgoing_IR_flux(T)
        meridional_heat_transport = self.compute_meridional_heat_transport(T)

        # Calculation of the temperature tendency according to the energy balance equation
        dT_dt = (absorbed_solar_flux - outgoing_IR_flux - meridional_heat_transport) / self.heat_capacity
        return dT_dt
    
    def global_temperature(self, T: FloatArray) -> float:
        """
        Calculation of the area-weighted global temperature.
        T : array of temperatures by latitude (in K)
        Returns the global temperature (in K).
        """
        return np.dot(T, self.f_area)

    def integrate(
        self, 
        T0: FloatArray = None, 
        dt_days: float = None, 
        n_years: int = None, 
        stop_at_convergence=False, 
        tol = 1e-4, 
        min_years = 5.0
    ) -> tuple[FloatArray, FloatArray, dict[str, Any]]:
        """
        Integration of the EBM over time using an explicit Euler method.
        
        Args:
            T0 : initial temperature profile by latitude (in K). If None, use the default from the input dataset.
            dt_days : time step in days. If None, use the default from config.
            n_years : maximum number of years to integrate. If None, use the default from config.
            stop_at_convergence : if True, stop the integration when the maximum temperature change is below 'tol' after 'min_years'. Default is False.
            tol : threshold for convergence in K. Default is 1e-4.
            min_years : minimum number of years to integrate before checking for convergence. Default is 5.0.
        
        Returns: 
            tuple (t_days, T_history, info) where:
            t_days : array of time points in days            
            T_history : array of temperature profiles by latitude at each time point
            info : dictionary with additional information about the integration (converged or not, number of steps per year and the total time in years)
        """

        T = T0.copy() if T0 is not None else self.input_data.T0.copy()
        dt_days = dt_days if dt_days is not None else self.dt_days
        n_years = n_years if n_years is not None else self.n_years
        
        if (self.YEAR_IN_DAYS % dt_days) != 0:
            raise ValueError(f"dt_days must divide evenly into the number of days in a year ({self.YEAR_IN_DAYS}) to ensure consistent time steps. Got dt_days={dt_days}.")
        
        dt = dt_days * 24.0 * 3600.0
        steps_per_year = int(round(self.YEAR_IN_DAYS / dt_days))
        max_steps = n_years * steps_per_year
        min_steps = min_years * steps_per_year
        

        T_history = np.zeros((max_steps + 1, self.n_lat), dtype=T.dtype)
        t_days = np.zeros(max_steps + 1, dtype=float)

        T_history[0] = T
        t_days[0] = 0.0

        converged = False
        last_index = max_steps


        for n in range(max_steps):
            dTdt = self.compute_temperature_tendency(n * dt_days, T)
            T_new = T + dt * dTdt

            T_history[n + 1] = T_new
            t_days[n + 1] = (n + 1) * dt_days

            if  (n + 1) % steps_per_year == 0 and (n + 1) >= min_steps:
                past_year_idx = (n + 1) - steps_per_year
                max_T_change = np.max(np.abs(T_new - T_history[past_year_idx]))

                if max_T_change < tol:
                    converged = True
                    if stop_at_convergence:
                        last_index = n + 1
                        T = T_new
                        break

            T = T_new

        info = {
            "converged": converged,
            "steps_per_year": steps_per_year,
            "t_final_years": t_days[last_index] / self.YEAR_IN_DAYS
        }
        return t_days[:last_index + 1], T_history[:last_index + 1], info
    

if __name__ == "__main__":
    print("This is the EBM1DBudyko model implementation. Please run the 'run_model.py' script to execute the model.")