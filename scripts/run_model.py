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

    temperature_last_year_per_day = Thist[-365:]
    temperature_last_year = (18*[0])
 

    i = 0

    for day_temp in temperature_last_year_per_day:
        temperature_last_year = temperature_last_year + day_temp
        i = i + 1

    print (i)
    temperature_last_year = temperature_last_year / 365

    planetary_temperature_last_year = model.global_mean_temperature(temperature_last_year)

    print(f"Température moyenne globale sur la dernière année = {planetary_temperature_last_year - 273.15:.2f} °C")

    # -----------------------------
    # Graphiques
    # -----------------------------
    plt.figure(figsize=(8, 5))
    plt.plot(model.lat_centers_deg, model.input_data.T0 - 273.15, label="Initiale")
    plt.plot(model.lat_centers_deg, temperature_last_year - 273.15, label="Finale")
    plt.xlabel("Latitude (°)")
    plt.ylabel("Température (°C)")
    plt.title("EBM 1D de Budyko - profil latitudinal")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    # plt.show()




    temp_lat_45 = [Thist[time][13] - 273.15 for time in range(len(Thist))]
    temp_lat_5 = [Thist[time][9] - 273.15 for time in range(len(Thist))]
    temp_lat_85 = [Thist[time][17] - 273.15 for time in range(len(Thist))]

    plt.figure(figsize=(8, 5))
    plt.plot((t[-730:] - t[-730])*365 , temp_lat_45[-730:], label="Température à 45° de latitude") 
    plt.ylim(-10, 20)
    plt.xlabel("Temps (années)")
    plt.ylabel("Température à 45° de latitude (°C)")
    plt.title("EBM 1D de Budyko - évolution temporelle 2 dernières années")
    plt.grid(True)
    plt.tight_layout()

    plt.figure(figsize=(8, 5))
    plt.plot((t[-730:] - t[-730])*365 , temp_lat_5[-730:], label="Température à 5° de latitude") 
    plt.ylim(-20, 50)
    plt.xlabel("Temps (années)")
    plt.ylabel("Température à 5° de latitude (°C)")
    plt.title("EBM 1D de Budyko - évolution temporelle 2 dernières années")
    plt.grid(True)
    plt.tight_layout()

    plt.figure(figsize=(8, 5))
    plt.plot((t[-730:] - t[-730])*365 , temp_lat_85[-730:], label="Température à 85° de latitude") 
    plt.ylim(-35, 5)
    plt.xlabel("Temps (années)")
    plt.ylabel("Température à 85° de latitude (°C)")
    plt.title("EBM 1D de Budyko - évolution temporelle 2 dernières années")
    plt.grid(True)
    plt.tight_layout()



    temp_lat_45_last_year = temp_lat_45[-365:]
    max_temp_lat_45_last_year = max(temp_lat_45_last_year)
    day_at_max_temp_lat_45_last_year = temp_lat_45_last_year.index(max_temp_lat_45_last_year)
    print(f"Température maximale à 45° de latitude lors de la dernière année : {max_temp_lat_45_last_year:.2f} °C")
    print(f"Jour de l'année où la température à 45° de latitude est maximale : {day_at_max_temp_lat_45_last_year} (0 = 1er janvier)")

    plt.show()

if __name__ == "__main__":
    main()