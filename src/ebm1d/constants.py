"""
Constants for the 1D Energy Balance Model.
All units are in SI (International System of Units) unless specified otherwise.
"""

# --- Albedo reference values ---
ALBEDO_PLANETARY = 0.3      # Average albedo of the Earth [-]    
ALBEDO_ATMOSPHERE = 0.07    # Atmospheric albedo (Rayleigh scattering) [-]
ALBEDO_CLOUD = 0.5          # Albedo of clouds [-]
ALBEDO_LAND = 0.25          # Albedo of land surfaces [-]
ALBEDO_OCEAN = 0.07         # Albedo of ocean surfaces [-]
ALBEDO_SNOW = 0.8           # Albedo of snow surfaces [-]
ALBEDO_ICE = 0.6            # Albedo of ice surfaces [-]


# --- Fractional coverage of surface types (for albedo calculations) ---
FRACTION_LAND = 0.3         # Fraction of Earth's surface covered by land [-]
FRACTION_CLOUD = 0.5        # Average cloud cover fraction [-]

# --- Temperature thresholds for phase changes ---
TEMP_INF_ICE = 272.15       # Temperature threshold for complete ice formation [K]
TEMP_SUP_ICE = 274.15       # Temperature threshold for full ice melting [K]
TEMP_INF_SNOW = 272.15      # Temperature threshold for complete snow formation [K]
TEMP_SUP_SNOW = 274.15      # Temperature threshold for full snow melting [K]

# --- Physical constants ---
SIGMA = 5.670373e-8          # Stefan-Boltzmann constant [W m^-2 K^-4]
SOLAR_CONSTANT = 1361.0      # Total solar irradiance [W m^-2]

# --- Outgoing IR flux parameters ---
P_CO2 = 430.0                # Current CO2 partial pressure [ppm] # This value is temporary and should be read from a configuration file
P_CO2_0 = 280.0              # Reference CO2 partial pressure (pre-industrial) [ppm]
IR_C = 54.13                 # Linearization coefficient for cloud outgoing IR flux [W m^-2]
IR_D = 0.58                  # Temperature dependence for cloud outgoing IR flux [W m^-2 K^-1]

# --- Time constants (Conversion) ---
SECONDS_IN_DAY = 86400       # [s]
DAYS_IN_YEAR = 365.25        # [days]
SECONDS_IN_YEAR = SECONDS_IN_DAY * DAYS_IN_YEAR  # [s]



