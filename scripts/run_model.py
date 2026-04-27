"""
Run the EBM1DBudyko model in its seasonal configuration.

This script integrates the model for 50 years using the default configuration
and saves the temperature history and figures in results/run_model.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    'figure.figsize': (8, 5),
    'lines.linewidth': 2,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'figure.constrained_layout.use': True,
    'savefig.dpi': 300,     
    'savefig.bbox': 'tight'
})

from pathlib import Path

from ebm1d import EBM1DBudyko
from ebm1d import load_config
from ebm1d import load_input_data

PROJECT_ROOT = Path(__file__).resolve().parents[1]

def main():  
    # Create results directory to save the figures
    results_dir = PROJECT_ROOT / "results" / "run_model"
    results_dir.mkdir(parents=True, exist_ok=True)

    # Load configuration and input data
    config = load_config(PROJECT_ROOT / "configs" / "default.toml")
    inputs = load_input_data(config.input.dataset, config.model.n_lat)

    # Creating an instance of the seasonal EBM
    ebm_seasonal = EBM1DBudyko(config, inputs, seasonal=True)

    # Integrating the model over 50 years
    t_days, T_history, info = ebm_seasonal.integrate(
        n_years=50,
        dt_days=1
    )

    # Let's display the convergence information
    if info['converged']:
        print(f"Seasonal EBM has converged. Final state reached at t = {info['t_final_years']} " 
              f"years after {info['steps_per_year'] * info['t_final_years']} steps.")
    else:
        print(f"Seasonal EBM has not converged !")


    lat_centers_deg = ebm_seasonal.lat_centers_deg

    # Create a DataFrame to store the temperature history for easier analysis
    T_history_df = pd.DataFrame(
        T_history,
        columns=[f"T_lat_{lat:.1f}" for lat in lat_centers_deg]
    )
    T_history_df.insert(0, "t_days", t_days)
    T_history_df.to_csv(results_dir / "temperature_history.csv", index=False)

    # Global temperature history
    T_global = ebm_seasonal.global_temperature(T_history)


    # Extracting the temperature history of last year
    steps_per_year = info['steps_per_year']
    T_last_year = T_history[-steps_per_year:]

    # Mean temperature profile of last year 
    T_last_year_mean = T_last_year.mean(axis=0)


    # Plotting the initial and final temperature profiles in °C
    fig, ax = plt.subplots()
    ax.plot(lat_centers_deg, ebm_seasonal.input_data.T0 - 273.15, label="Initial condition", ls='--', marker='o', color='gray')
    ax.plot(lat_centers_deg, T_last_year_mean - 273.15, label="Seasonal model", marker='o')
    ax.set_xlabel("Latitude (°)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("Latitudinal temperature profile at equilibrium")
    ax.legend()
    fig.savefig(results_dir / "temperature_profile.png")
    plt.close(fig)


    # Converting time from days to years for better readability on the plot
    YEAR_IN_DAYS = ebm_seasonal.YEAR_IN_DAYS 
    # Note that .YEAR_IN_DAYS is a class attribute of EBM1DBudyko (i.e., not instance-specific), so it can also be accessed directly from the class.
    t_years = t_days / YEAR_IN_DAYS

    # Plotting the global temperature history in °C
    fig, ax = plt.subplots()
    ax.plot(t_years, T_global - 273.15)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Global temperature (°C)")
    ax.set_title("Global temperature history – seasonal EBM")
    fig.savefig(results_dir / "global_temperature_history.png")
    plt.close(fig)


    # Calculating and displaying the global mean temperature of the last year
    T_global_mean_last_year = T_global[-steps_per_year:].mean()
    print(f"Global mean temperature of the last year: {T_global_mean_last_year - 273.15:.2f} °C")
    print(f"All results have been saved in: {results_dir}")

if __name__ == "__main__":
    main()