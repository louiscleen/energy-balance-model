import numpy as np
import numpy.typing as npt
import ebm1d.constants as cst

import matplotlib.pyplot as plt

# pas oublier de pondérer avec surface
# vérif floatarry

FloatArray = npt.NDArray[np.float64]

class EBM1DBudyko:
    def __init__(self, n_latitudes=18, seasonal=False):
        self.n_lat = n_latitudes
        self.lat_width_deg = 180.0 / self.n_lat
        lat_edges_rad = np.linspace(-np.pi / 2, np.pi / 2, self.n_lat + 1)
        self.lat_centers_rad = (lat_edges_rad[:-1] + lat_edges_rad[1:]) / 2
        self.lat_deg = np.degrees(self.lat_centers_rad)
        
        self.seasonal = seasonal
        if self.seasonal:
            self._get_solar_flux = self._compute_solar_flux_seasonal
        else:
            self._precompute_solar_flux = self._compute_solar_flux_annual_mean()
            self._get_solar_flux = lambda day_of_year: self._precompute_solar_flux

        self.C = 2.0e8                # capacité thermique surfacique [J m^-2 K^-1]

    def _compute_fraction_ice(self, T: FloatArray) -> FloatArray:
        return np.clip((cst.TEMP_SUP_ICE - T) / (cst.TEMP_SUP_ICE - cst.TEMP_INF_ICE), 0.0, 1.0)

    def _compute_fraction_snow(self, T: FloatArray) -> FloatArray:
        return np.clip((cst.TEMP_SUP_SNOW - T) / (cst.TEMP_SUP_SNOW - cst.TEMP_INF_SNOW), 0.0, 1.0)

    def _compute_albedo(self, T: FloatArray) -> FloatArray:
        """
        Calculation of albedo as a function of temperature.
        T : temperature table by latitude (in K)
        Returns a table of the same format with the albedo (dimensionless).
        """
        self.f_snow = self._compute_fraction_snow(T)
        self.f_ice = self._compute_fraction_ice(T)

        #
        albedo_continent = self.f_snow * cst.ALBEDO_SNOW + (1 - self.f_snow) * cst.ALBEDO_LAND
        albedo_ocean = self.f_ice * cst.ALBEDO_ICE + (1 - self.f_ice) * cst.ALBEDO_OCEAN
        albedo_surf = np.array(cst.FRACTION_LAND) * albedo_continent + (1 - np.array(cst.FRACTION_LAND)) * albedo_ocean

        # Albedo under clear sky: considering the two first terms of the multiple reflection series between 
        # the surface and the atmosphere (through Rayleigh scattering)
        albedo_clear_sky = cst.ALBEDO_ATMOSPHERE + (1 - cst.ALBEDO_ATMOSPHERE)**2 * albedo_surf
        # Albedo under cloudy sky: considering that the top of clouds are higher in the atmosphere, 
        # we assume that the albedo is a sum of a fraction of the atmospheric albedo and the cloud albedo.
        albedo_cloudy_sky = cst.ALBEDO_ATMOSPHERE / 2.0 + (1 - cst.ALBEDO_ATMOSPHERE / 2.0) * cst.ALBEDO_CLOUD
        
        # Final albedo: weighted average between clear sky and cloudy sky
        albedo = np.array(cst.FRACTION_CLOUD) * albedo_cloudy_sky + (1 - np.array(cst.FRACTION_CLOUD)) * albedo_clear_sky
        return albedo
    
    def _compute_solar_flux_annual_mean(self) -> FloatArray:
        """
        Calcul du flux solaire incident moyen annuel par latitude.
        Retourne un tableau de taille n_lat avec le flux en W m^-2.
        """
        # Insolation moyenne annuelle simplifiée : Q = S0/4 * (1 + 0.241 * (1.5 * sin^2(lat) - 0.5))
        #x = np.sin(self.lat_centers_rad)
        #Q = (cst.SOLAR_CONSTANT / 4.0) * (1.0 + 0.241 * (1.5 * x**2 - 0.5))

        Q = np.array(cst.S_s4)   # Utilisation des valeurs pré-calculées de S * s/4 pour chaque bande de latitude

        return Q

    # def _compute_solar_flux_seasonal(self, day_of_year) -> FloatArray:
    #     return
    
    def _compute_outgoing_IR_flux(self, T: FloatArray) -> FloatArray:
        """
        Calcul du flux de rayonnement infrarouge sortant selon une paramétrisation linéaire.
        T : tableau des températures par latitude (en K)
        Retourne un tableau du même format avec le flux IR sortant (en W m^-2).
        """

        A = 243.39 - 4.48 * np.log(cst.P_CO2 / cst.P_CO2_0)
        B = 2.07 - 0.0514 * np.log(cst.P_CO2 / cst.P_CO2_0)
        outgoing_IR_flux = (A + B * (T - 273.15)) - (cst.IR_C + cst.IR_D * (T - 273.15)) * np.array(cst.FRACTION_CLOUD)
    
        return outgoing_IR_flux

    def _compute_meridional_heat_transport(self, T: FloatArray) -> FloatArray:
        """
        Calcul du transport de chaleur meridional selon la paramétrisation de Budyko.
        T : tableau des températures par latitude (en K)
        Retourne un tableau du même format avec le transport de chaleur (en W m^-2).
        """
        # Différence de température entre les latitudes adjacentes
        T_mean = np.dot(T, cst.FRACTION_AREA)  # Température moyenne globale pondérée par l'aire
        # Transport de chaleur proportionnel à la différence de température
        transport = cst.K_T * (T - T_mean)
        # On ajoute une valeur nulle aux extrémités pour correspondre à la taille de T

        # Vérification de la conservation de l'énergie : le transport total doit être nul
        # à faire
        return transport
    
    def compute_temperature_tendency(self, t, T: FloatArray) -> FloatArray: #compute_delta_T
        """
        Calcul de la tendance de température dT/dt selon l'équation de l'EBM.
        T : tableau des températures par latitude (en K)
        Retourne un tableau du même format avec la tendance de température (en K/s).
        """
        
        Q = self._get_solar_flux(t)
        alpha = self._compute_albedo(T)
        absorbed_solar_flux = Q * (1 - alpha)
        outgoing_IR_flux = self._compute_outgoing_IR_flux(T)
        meridional_heat_transport = self._compute_meridional_heat_transport(T)

        # Équation de l'EBM : C * dT/dt = absorbed_solar_flux - outgoing_IR_flux + divergence of heat transport
        # On peut utiliser np.diff pour calculer la différence entre les latitudes adjacentes, en ajoutant une valeur nulle aux extrémités.

        dT_dt = (absorbed_solar_flux - outgoing_IR_flux - meridional_heat_transport) / self.C
        return dT_dt
    
    def global_mean_temperature(self, T: FloatArray) -> float:
        """
        Calcul de la température moyenne globale pondérée par l'aire.
        T : tableau des températures par latitude (en K)
        Retourne la température moyenne globale (en K).
        """
        return np.dot(T, cst.FRACTION_AREA)

    def integrate(self, T0, years=50, dt_days=5.0):
        dt = dt_days * 24.0 * 3600.0
        nsteps = int(years * 365.0 / dt_days)

        T = T0.copy()
        history = np.zeros((nsteps + 1, self.n_lat))
        t_years = np.zeros(nsteps + 1)

        history[0] = T

        for n in range(1, nsteps + 1):
            dTdt = self.compute_temperature_tendency(n * dt, T)
            T = T + dt * dTdt

            history[n] = T
            t_years[n] = n * dt_days / 365.0

        return t_years, history
    

if __name__ == "__main__":
    model = EBM1DBudyko(n_latitudes=18, seasonal=False)
    # print("Latitudes des bandes (en degrés) :")
    # for i, lat in enumerate(model.lat_centers_rad): 
    #     print(f"Band {i}: Latitude = {np.degrees(lat):.2f}°")

    # Intégration
    t, Thist = model.integrate(
        T0=np.array(cst.INITIAL_TEMPERATURE),
        years=500,
        dt_days=5.0
    )

    Tfinal = Thist[-1]
    Tbar_final = model.global_mean_temperature(Tfinal)

    print("Température finale par bande (°C) :")
    print(np.round(Tfinal, 2))
    print(f"Température moyenne globale finale = {Tbar_final:.2f} °C")

    # -----------------------------
    # Graphiques
    # -----------------------------
    plt.figure(figsize=(8, 5))
    plt.plot(model.lat_deg, np.array(cst.INITIAL_TEMPERATURE) - 273.15, label="Initiale")
    plt.plot(model.lat_deg, Tfinal - 273.15, label="Finale")
    plt.xlabel("Latitude (°)")
    plt.ylabel("Température (°C)")
    plt.title("EBM 1D de Budyko - profil latitudinal")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    area = np.array(cst.AREA)
    f_area = np.array(cst.FRACTION_AREA)

    # plt.figure(figsize=(8, 5))
    # gmst = np.sum(Thist * f_area[None, :], axis=1) - 273.15
    # plt.plot(t, gmst)
    # plt.xlabel("Temps (années)")
    # plt.ylabel("Température moyenne globale (°C)")
    # plt.title("Évolution de la température moyenne globale")
    # plt.grid(True)
    # plt.tight_layout()
    # plt.show()


    plt.figure(figsize=(8, 5))
    plt.plot(model.lat_deg, np.array(cst.FRACTION_CLOUD), label="Fraction de couverture nuageuse")
    plt.xlabel("Latitude (°)")
    plt.ylabel("Fraction de couverture nuageuse")
    plt.ylim(0, 1)   # force l'axe y entre 0 et 1
    plt.title("Fraction de couverture nuageuse par latitude")
    plt.legend()
    plt.tight_layout()
    plt.show()

    meridional_heat_transport = model._compute_meridional_heat_transport(Tfinal) * area  # Transport de chaleur par latitude, pondéré par l'aire pour obtenir le transport total en W
    plt.figure(figsize=(8, 5))
    plt.plot(model.lat_deg, meridional_heat_transport, label="Transport de chaleur meridional")
    plt.xlabel("Latitude (°)")
    plt.ylabel("Transport de chaleur (W)")
    plt.title("Transport de chaleur meridional par latitude")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 5))
    plt.plot

    print(f"Fraction de neige par latitude : {np.round(model.f_snow, 2)}")
    print(model.f_snow)

    print("Fraction de glace par latitude :")
    print(model.f_ice)
