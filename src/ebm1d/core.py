import numpy as np
import ebm1d.constants as cst

class EBM1DBudyko:
    def __init__(self, n_latitudes=18):
        self.n_lat = n_latitudes
        self.lat_width_deg = 180.0 / self.n_lat
        lat_edges_rad = np.linspace(-np.pi / 2, np.pi / 2, self.n_lat + 1)
        self.lat_centers_rad = (lat_edges_rad[:-1] + lat_edges_rad[1:]) / 2

        self.C = 2.0e8                # capacité thermique surfacique [J m^-2 K^-1]
        self.k_t = 1.0e6              # coefficient de transport de chaleur [W m^-2 K^-1]

    def compute_fraction_ice(self, T):
        return np.clip((cst.TEMP_SUP_ICE - T) / (cst.TEMP_SUP_ICE - cst.TEMP_INF_ICE), 0.0, 1.0)

    def compute_fraction_snow(self, T):
        return np.clip((cst.TEMP_SUP_SNOW - T) / (cst.TEMP_SUP_SNOW - cst.TEMP_INF_SNOW), 0.0, 1.0)

    def compute_albedo(self, T):
        """
        Calculation of albedo as a function of temperature.
        T : temperature table by latitude (in K)
        Returns a table of the same format with the albedo (dimensionless).
        """
        f_snow = self.compute_fraction_snow(T)
        f_ice = self.compute_fraction_ice(T)

        #
        albedo_continent = f_snow * cst.ALBEDO_SNOW + (1 - f_snow) * cst.ALBEDO_LAND
        albedo_ocean = f_ice * cst.ALBEDO_ICE + (1 - f_ice) * cst.ALBEDO_OCEAN
        albedo_surf = cst.FRACTION_LAND * albedo_continent + (1 - cst.FRACTION_LAND) * albedo_ocean

        # Albedo under clear sky: considering the two first terms of the multiple reflection series between 
        # the surface and the atmosphere (through Rayleigh scattering)
        albedo_clear_sky = cst.ALBEDO_ATMOSPHERE + (1 - cst.ALBEDO_ATMOSPHERE)**2 * albedo_surf
        # Albedo under cloudy sky: considering that the top of clouds are higher in the atmosphere, 
        # we assume that the albedo is a sum of a fraction of the atmospheric albedo and the cloud albedo.
        albedo_cloudy_sky = cst.ALBEDO_ATMOSPHERE / 2.0 + (1 - cst.ALBEDO_ATMOSPHERE / 2.0) * cst.ALBEDO_CLOUD
        
        # Final albedo: weighted average between clear sky and cloudy sky
        albedo = cst.FRACTION_CLOUD * albedo_cloudy_sky + (1 - cst.FRACTION_CLOUD) * albedo_clear_sky
        return albedo
    
    def compute_solar_flux_annual_mean(self):
        """
        Calcul du flux solaire incident moyen annuel par latitude.
        Retourne un tableau de taille n_lat avec le flux en W m^-2.
        """
        # Insolation moyenne annuelle simplifiée : Q = S0/4 * (1 + 0.241 * (1.5 * sin^2(lat) - 0.5))
        x = np.sin(self.lat_centers_rad)
        Q = (cst.SOLAR_CONSTANT / 4.0) * (1.0 + 0.241 * (1.5 * x**2 - 0.5))
        return Q

    def compute_solar_flux_seasonal(self, day_of_year):
        return
    
    def compute_outgoing_IR_flux(self, T):
        """
        Calcul du flux de rayonnement infrarouge sortant selon une paramétrisation linéaire.
        T : tableau des températures par latitude (en K)
        Retourne un tableau du même format avec le flux IR sortant (en W m^-2).
        """

        A = 243.39 - 4.48 * np.log(cst.P_CO2 / cst.P_CO2_0)
        B = 2.07 - 0.0514 * np.log(cst.P_CO2 / cst.P_CO2_0)
        outgoing_IR_flux = A + B * (T - 273.15) + cst.IR_C + cst.IR_D * (T - 273.15)
    
        return outgoing_IR_flux

    def compute_meridional_heat_transport(self, T):
        """
        Calcul du transport de chaleur meridional selon la paramétrisation de Budyko.
        T : tableau des températures par latitude (en K)
        Retourne un tableau du même format avec le transport de chaleur (en W m^-2).
        """
        # Différence de température entre les latitudes adjacentes
        T_mean = np.mean(T)
        # Transport de chaleur proportionnel à la différence de température
        transport = -self.k_t * (T - T_mean)
        # On ajoute une valeur nulle aux extrémités pour correspondre à la taille de T

        # Vérification de la conservation de l'énergie : le transport total doit être nul
        # à faire
        return transport
    