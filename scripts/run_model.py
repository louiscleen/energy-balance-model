import numpy as np
import matplotlib.pyplot as plt
import tomllib

from ebm1d import EBM1DBudyko
from ebm1d import load_config
from ebm1d import load_input_data

def main():
    print("Hello from energy-balance-model!")

    config = load_config("configs/default.toml")
    inputs = load_input_data(config.data.dataset, config.model.n_lat)

    print("Configuration loaded:")

    model = EBM1DBudyko(config, inputs)


    # model = EBM1DBudyko(n_latitudes=18, seasonal=False)
    # print("Latitudes des bandes (en degrés) :")
    # for i, lat in enumerate(model.lat_centers_rad): 
    #     print(f"Band {i}: Latitude = {np.degrees(lat):.2f}°")





    t, Thist = model.integrate(
        # T0=np.array(cst.INITIAL_TEMPERATURE),
        # years=500,
        # dt_days=5.0
    )

    Tfinal = Thist[-1]
    Tbar_final = model.global_mean_temperature(Tfinal)

    print(f"Température moyenne globale finale = {Tbar_final:.2f} °C")

    # -----------------------------
    # Graphiques
    # -----------------------------
    plt.figure(figsize=(8, 5))
    plt.plot(model.lat_centers_deg, model.input_data.T0 - 273.15, label="Initiale")
    plt.plot(model.lat_centers_deg, Tfinal - 273.15, label="Finale")
    plt.xlabel("Latitude (°)")
    plt.ylabel("Température (°C)")
    plt.title("EBM 1D de Budyko - profil latitudinal")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()






if __name__ == "__main__":
    main()