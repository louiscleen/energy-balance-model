import numpy as np

from ebm1d.core import EBM1DBudyko

def main():
    print("Hello from energy-balance-model!")

    model = EBM1DBudyko(n_latitudes=18, seasonal=False)
    print("Latitudes des bandes (en degrés) :")
    for i, lat in enumerate(model.lat_centers_rad): 
        print(f"Band {i}: Latitude = {np.degrees(lat):.2f}°")


if __name__ == "__main__":
    main()