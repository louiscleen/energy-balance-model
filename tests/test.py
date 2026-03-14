import numpy as np
from ebm1d import EBM1DBudyko

ebm = EBM1DBudyko(n_latitudes=18)

def print_latitudes():
    for i, lat in enumerate(ebm.lat_centers_rad):
        print(f"Band {i}: Latitude = {np.degrees(lat):.2f}°")

if __name__ == "__main__":
    print("Latitudes des bandes (en degrés) :")
    print_latitudes()