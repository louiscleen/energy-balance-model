"""
Script to run the calculation of entropy production as a function of the meridional heat transport coefficient (kt).
The script integrates the seasonal EBM for a range of kt values, computes the entropy production for each case, and plots the results.
The maximum entropy production give a theoretical estimate of the optimal kt value for the model (see Paltridge 1975).
"""

import numpy as np
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
    print("--- Calculation of entropy production as a function of kt ---")

    results_dir = PROJECT_ROOT / "results" / "entropy_production"
    results_dir.mkdir(parents=True, exist_ok=True)

    config = load_config(PROJECT_ROOT / "configs" / "default.toml")
    inputs = load_input_data(config.input.dataset, config.model.n_lat)

    ebm_seasonal = EBM1DBudyko(config, inputs, seasonal=True)


    kt = np.arange(0, 6.1, 0.1, dtype=float)
    area = ebm_seasonal.area
    sigma = np.zeros(len(kt))

    
    fig, ax = plt.subplots()

    # Loop over kt values to compute entropy production and plot temperature profiles
    for i, kt_value in enumerate(kt):
        ebm_seasonal.K_meridional = kt_value
        t_days, history, info = ebm_seasonal.integrate(
            stop_at_convergence=True,
            n_years=100,
            dt_days=1.0,
            min_years=10.0
        )
        T_last_year = history[-365:]
        T_global_last_year = ebm_seasonal.global_temperature(T_last_year)

        # Compute the entropy production for the last year using the formula: sigma = kt * sum(area * (T_global - T) / T)
        # Since the model is seasonal, we compute the entropy production for each day and then average over the last year
        T_ratio_last_year = ((T_global_last_year[:, np.newaxis] - T_last_year) / T_last_year)
        sigma_last_year = kt_value * (T_ratio_last_year @ area)
        sigma[i] = sigma_last_year.mean()

        # Print progress every ~10%
        progress = (i + 1) / len(kt) * 100
        if (progress % 10) < (100 / len(kt)):
            print(f"Progress: {i+1}/{len(kt)} ({progress:.1f}%)")

        # Plot the latitudinal temperature profile for the last year for integer kt values (so as not to overcrowd the plot)
        if (kt_value % 1.0) == 0.0:
            ax.plot(ebm_seasonal.lat_centers_deg, T_last_year.mean(axis=0) - 273.15, label=fr"$k_t = {kt_value:.1f}$", marker='o')

    ax.set_xlabel("Latitude (°)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("Latitudinal temperature profiles for different kt values - Seasonal EBM")
    ax.legend(loc='lower center')
    plt.savefig(results_dir / "temperature_profiles_kt.pdf")
    plt.savefig(results_dir / "temperature_profiles_kt.png")
    plt.close(fig)

    # Find the maximum entropy production and corresponding kt
    idx_sigma_max = np.argmax(sigma)
    kt_optimal = kt[idx_sigma_max]
    print("-"*25)
    print(f"The maximum entropy production is: {sigma[idx_sigma_max]:.2E} W/K for kt = {kt_optimal:.3f} W/(m²K)")

    # Plot entropy production as a function of kt
    fig, ax = plt.subplots()
    ax.plot(kt, sigma, marker='o')
    ax.axvline(kt_optimal, color='tab:gray', linestyle='--')
    ax.text(kt_optimal + 0.1, 0.0, f"$k_t = {kt_optimal:.2f}$",
        ha='left', va='bottom')
    ax.set_xlabel(r"$k_t$ (W/(m²K))")
    ax.set_ylabel(r"$\sigma$ (W/K)")
    ax.set_title(r"Entropy production as a function of $k_t$ - Seasonal EBM")
    plt.savefig(results_dir / "entropy_production_kt.png")
    plt.savefig(results_dir / "entropy_production_kt.pdf")
    plt.close(fig)

    print(f"Figures saved in {results_dir} directory.")


if __name__ == "__main__":
    main()