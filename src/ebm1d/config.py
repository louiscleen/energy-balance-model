"""
Configuration management for the 1D Energy Balance Model (EBM).
This module defines data classes for each section of the configuration and a function to load the configuration from a TOML file.
"""

from dataclasses import dataclass
from pathlib import Path
import tomllib

# --------------------------------------------------
# Model configuration
# --------------------------------------------------

@dataclass(frozen=True)
class ModelConfig:
    n_lat: int

    def __post_init__(self):
        if self.n_lat < 1:
            raise ValueError("Number of latitude bands must be at least 1.") 


# --------------------------------------------------
# Time integration
# --------------------------------------------------

@dataclass(frozen=True)
class TimeConfig:
    dt: float
    n_years: int

    def __post_init__(self):
        if self.dt <= 0:
            raise ValueError("Time step must be positive.")
        if self.n_years < 1:
            raise ValueError("Number of years must be at least 1.")


# --------------------------------------------------
# Surface heat capacity
# --------------------------------------------------

@dataclass(frozen=True)
class HeatCapacityConfig:
    C_value: float

    def __post_init__(self):
        if self.C_value <= 0:
            raise ValueError("Heat capacity must be positive.")


# --------------------------------------------------
# Albedo parameters
# --------------------------------------------------

@dataclass(frozen=True)
class AlbedoConfig:
    alpha_atmosphere: float
    alpha_cloud: float
    alpha_land: float
    alpha_ocean: float
    alpha_snow: float
    alpha_ice: float
    T_ice_min: float
    T_ice_max: float
    T_snow_min: float
    T_snow_max: float
    
    def __post_init__(self):
        for attr in ["alpha_atmosphere", "alpha_cloud", "alpha_land", "alpha_ocean", "alpha_snow", "alpha_ice"]:
            value = getattr(self, attr)
            if not (0.0 <= value <= 1.0):
                raise ValueError(f"{attr} must be between 0 and 1.")
        
        if self.T_ice_min >= self.T_ice_max:
            raise ValueError("T_ice_min must be less than T_ice_max.")
        if self.T_snow_min >= self.T_snow_max:
            raise ValueError("T_snow_min must be less than T_snow_max.")
        

# --------------------------------------------------
# Insolation
# --------------------------------------------------

@dataclass(frozen=True)
class InsolationConfig:
    S0: float
    seasonal : bool

    def __post_init__(self):
        if self.S0 <= 0:
            raise ValueError("Solar constant must be positive.")


# --------------------------------------------------
# Radiation parameters
# --------------------------------------------------

@dataclass(frozen=True)
class RadiationConfig:
    p_CO2: float
    A1: float
    A2: float
    B1: float
    B2: float
    C: float
    D: float

    def __post_init__(self):
        if self.p_CO2 <= 0:
            raise ValueError("CO2 concentration must be positive.")


# --------------------------------------------------
# Heat transport
# --------------------------------------------------

@dataclass(frozen=True)
class TransportConfig:
    K_meridional: float

    def __post_init__(self):
        if self.K_meridional <= 0:
            raise ValueError("Meridional heat transport coefficient must be positive.")


# --------------------------------------------------
# Input dataset
# --------------------------------------------------

@dataclass(frozen=True)
class DataConfig:
    dataset: Path


# --------------------------------------------------
# Output configuration
# --------------------------------------------------

@dataclass(frozen=True)
class OutputConfig:
    save_every: int
    output_dir: Path

    def __post_init__(self):
        if self.save_every < 1:
            raise ValueError("save_every must be at least 1.")


# --------------------------------------------------
# Global configuration object
# --------------------------------------------------

@dataclass(frozen=True)
class EBM1DConfig:

    model: ModelConfig
    time: TimeConfig
    heat_capacity: HeatCapacityConfig
    albedo: AlbedoConfig
    insolation: InsolationConfig
    radiation: RadiationConfig
    transport: TransportConfig
    data: DataConfig
    output: OutputConfig


def _require_section(data: dict, section: str) -> dict:
    try:
        value = data[section]
    except KeyError as e:
        raise KeyError(f"Missing [{section}] section in config file.") from e

    if not isinstance(value, dict):
        raise TypeError(f"Section [{section}] must be a TOML table.")

    return value

def load_config(path: str | Path) -> EBM1DConfig:

    path = Path(path).resolve()

    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    
    with path.open("rb") as f:
        data = tomllib.load(f)

    base_dir = path.parents[1]

    model_data = _require_section(data, "model")
    time_data = _require_section(data, "time")
    heat_capacity_data = _require_section(data, "heat_capacity")
    albedo_data = _require_section(data, "albedo")
    insolation_data = _require_section(data, "insolation")
    radiation_data = _require_section(data, "radiation")
    transport_data = _require_section(data, "transport")
    data_data = _require_section(data, "data")
    output_data = _require_section(data, "output")

    cfg = EBM1DConfig(

        model=ModelConfig(**model_data),

        time=TimeConfig(**time_data),

        heat_capacity=HeatCapacityConfig(**heat_capacity_data),

        albedo=AlbedoConfig(**albedo_data),

        insolation=InsolationConfig(**insolation_data),

        radiation=RadiationConfig(**radiation_data),

        transport=TransportConfig(**transport_data),

        data=DataConfig(dataset=(base_dir / "data" / data_data["dataset"]).resolve()),

        output=OutputConfig(
            save_every=output_data["save_every"],
            output_dir=(base_dir / output_data["output_dir"]).resolve()
        ),
    )

    return cfg