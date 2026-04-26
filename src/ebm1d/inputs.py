"""
Input data management for the 1D Energy Balance Model (EBM).
This module defines the structure for input data and provides functions to load it from CSV files.
"""

from dataclasses import dataclass
import pandas as pd
import numpy as np
import numpy.typing as npt
from pathlib import Path

FloatArray = npt.NDArray[np.float64]


@dataclass(frozen=True)
class InputData:
    f_land: FloatArray
    f_cloud: FloatArray
    T0: FloatArray


def load_input_data(dataset_path: Path, n_lat: int) -> InputData:
    df = pd.read_csv(dataset_path, sep="\t")

    required = {"f_land", "f_cloud", "T0"}
    missing = required - set(df.columns)

    if missing: 
        raise ValueError(
            f"Missing columns in dataset: {sorted(missing)}"
        )

    if len(df) != n_lat:
        raise ValueError(
            f"Dataset length ({len(df)}) must match n_lat ({n_lat})"
        )

    return InputData(
        f_land=df["f_land"].to_numpy(dtype=np.float64),
        f_cloud=df["f_cloud"].to_numpy(dtype=np.float64),
        T0=df["T0"].to_numpy(dtype=np.float64),
    )