"""
Constants for the 1D Energy Balance Model.
All units are in SI (International System of Units) unless specified otherwise.
"""

# DEBUG
AREA = [
    3.8733E+12,
    1.15022E+13,
    1.87816E+13,
    2.54904E+13,
    3.14246E+13,
    3.6404E+13,
    4.02773E+13,
    4.29268E+13,
    4.4272E+13,
    4.4272E+13,
    4.29268E+13,
    4.02773E+13,
    3.6404E+13,
    3.14246E+13,
    2.54904E+13,
    1.87816E+13,
    1.15022E+13,
    3.8733E+12
]

FRACTION_AREA = [area / sum(AREA) for area in AREA]

# --- Albedo reference values ---
ALBEDO_PLANETARY = 0.3      # Average albedo of the Earth [-]    
#ALBEDO_ATMOSPHERE = 0.07    # Atmospheric albedo (Rayleigh scattering) [-]
ALBEDO_ATMOSPHERE = 0.06    # Atmospheric albedo (Rayleigh scattering) [-]
#ALBEDO_CLOUD = 0.5          # Albedo of clouds [-]
ALBEDO_CLOUD = 0.4081       # Albedo of clouds [-] ajusted

#ALBEDO_LAND = 0.25          # Albedo of land surfaces [-]
ALBEDO_LAND = 0.15          # Albedo of land surfaces [-]
#ALBEDO_OCEAN = 0.07         # Albedo of ocean surfaces [-]
ALBEDO_OCEAN = 0.08         # Albedo of ocean surfaces [-]
#ALBEDO_SNOW = 0.8           # Albedo of snow surfaces [-]
ALBEDO_SNOW = 0.7           # Albedo of snow surfaces [-]
#ALBEDO_ICE = 0.6            # Albedo of ice surfaces [-]
ALBEDO_ICE = 0.7            # Albedo of ice surfaces [-]


# --- Fractional coverage of surface types (for albedo calculations) ---
#FRACTION_LAND = 0.3         # Fraction of Earth's surface covered by land [-]
#FRACTION_CLOUD = 0.5        # Average cloud cover fraction [-]

# Fraction of land by latitude band (from 90°S to 90°N, in 10° increments) [-]
FRACTION_LAND = [1.00000, 0.73621, 0.09520, 0.00836, 0.03090, 0.11372, 0.23115, 0.21944, 0.23478, 0.22744, 0.26193, 0.37483, 0.42770, 0.52424, 0.57586, 0.70952, 0.29890, 0.09914]

# Monthly cloud cover fraction from CERES data (https://ceres.larc.nasa.gov/data/)
FRACTION_CLOUD = [0.47000, 0.58600, 0.74600, 0.76700, 0.65200, 0.53400, 0.47500, 0.47000, 0.48500, 0.48000, 0.44000, 0.46400, 0.55900, 0.63300, 0.67900, 0.67600, 0.65100, 0.64000]

# --- Temperature thresholds for phase changes ---
TEMP_INF_ICE = 257.0       # Temperature threshold for complete ice formation [K]
TEMP_SUP_ICE = 270.0       # Temperature threshold for full ice melting [K]
TEMP_INF_SNOW = 263.0      # Temperature threshold for complete snow formation [K]
TEMP_SUP_SNOW = 283.0      # Temperature threshold for full snow melting [K]

# --- Physical constants ---
SIGMA = 5.670373e-8          # Stefan-Boltzmann constant [W m^-2 K^-4]
#SOLAR_CONSTANT = 1361.0      # Total solar irradiance [W m^-2]
SOLAR_CONSTANT = 1368.0      # Total solar irradiance [W m^-2]

# --- Solar insolation parameters ---
S_s4 = [171.0, 181.602, 213.408, 263.34, 305.064, 349.182, 383.04, 406.638, 416.898, 416.898, 406.638, 383.04, 349.182, 305.064, 263.34, 213.408, 181.602, 171] # S * s/4

# --- Outgoing IR flux parameters ---
#P_CO2 = 430.0                # Current CO2 partial pressure [ppm] # This value is temporary and should be read from a configuration file
P_CO2 = 280.0                # Current CO2 partial pressure [ppm] (initial value, can be modified during the simulation)
P_CO2_0 = 280.0              # Reference CO2 partial pressure (pre-industrial) [ppm]
IR_C = 54.13                 # Linearization coefficient for cloud outgoing IR flux [W m^-2]
IR_D = 0.58                  # Temperature dependence for cloud outgoing IR flux [W m^-2 K^-1]

# --- Meridional heat transport parameters ---
K_T = 3.057                  # Meridional heat transport coefficient [W m^-2 K^-1]

# --- Initial conditions ---
# Initial temperature (global average) [K]
INITIAL_TEMPERATURE = [230.8, 243.6, 266.3, 276.0, 283.0, 289.6, 294.5, 297.7, 299.2, 299.6, 299.2, 296.0, 289.3, 281.9, 275.3, 268.0, 260.8, 256.2]


# --- Time constants (Conversion) ---
SECONDS_IN_DAY = 86400       # [s]
DAYS_IN_YEAR = 365.25        # [days]
SECONDS_IN_YEAR = SECONDS_IN_DAY * DAYS_IN_YEAR  # [s]




