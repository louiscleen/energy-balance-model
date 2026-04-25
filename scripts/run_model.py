import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ebm1d import EBM1DBudyko, model
from ebm1d import load_config
from ebm1d import load_input_data

from pathlib import Path

def main():
    BASE_DIR = Path(__file__).resolve().parent.parent

    config = load_config(BASE_DIR/"configs/default.toml")
    inputs = load_input_data(config.data.dataset, config.model.n_lat)

    print("Configuration loaded successfully.")

    model = EBM1DBudyko(config, inputs, seasonal=True)



    t, Thist, info = model.integrate(
        # T0=np.array(cst.INITIAL_TEMPERATURE),
        stop_at_convergence=True,
        n_years=100,
        dt_days=1.0
    )


 



if __name__ == "__main__":
    main()