import numpy as np
import pandas as pd
from pyproj import Transformer

zx1_15_path = "data/ZX1_2005-11-15_1105-1135_GGA-only.TXT"
zx1_16_path = "data/ZX1_2005-11-16_1100-1130_GGA-only.TXT"
zx2_15_path = "data/ZX1_2005-11-15_1105-1135_GGA-only.TXT"
zx2_16_path = "data/ZX1_2005-11-16_1100-1130_GGA-only.TXT"

def read_nmea(nmea_file):
    output = []

    with open(nmea_file, "r") as file:
        for line in file:
            line = line.strip()

            values = line.split(",")

            output.append([
                values[0],   # id
                values[1],   # time
                values[2],   # latitude
                values[4],   # longitude
                values[6],   # quality
                values[7],   # nsats
                values[8],   # hdop
                values[9],   # geoidheight
                values[11],  # geoidundulation
                values[13],  # dgps
                values[14]   # check
            ])

    df = pd.DataFrame(output, columns=[
        "id",
        "time",
        "latitude",
        "longitude",
        "quality",
        "nsats",
        "hdop",
        "geoidheight",
        "geoidundulation",
        "dgps",
        "check"
    ])

    df["latitude"] = df["time"].astype(float)
    df["latitude"] = df["latitude"].astype(float)
    df["latitude"] = df["longitude"].astype(float)
    df["latitude"] = df["quality"].astype(float)
    df["latitude"] = df["nsats"].astype(float)
    df["latitude"] = df["hdop"].astype(float)
    df["latitude"] = df["geoidheight"].astype(float)
    df["latitude"] = df["geoidundulation"].astype(float)

    return df

zx1_15_geo = read_nmea(zx1_15_path)
zx1_16_geo = read_nmea(zx1_16_path)
zx2_15_geo = read_nmea(zx2_15_path)
zx2_16_geo = read_nmea(zx2_16_path)

def geo_to_ecef(df):
    transformer = Transformer.from_crs(
        "EPSG:4979",
        "EPSG:4978",
        always_xy=True
    )

    h = df.geoidheight + df.geoidundulation

    x, y, z = transformer.transform(
        df.longitude.values,
        df.latitude.values,
        h.values
    )

    df = df.copy()
    df["X"] = x
    df["Y"] = y
    df["Z"] = z

    return df

def calc_deg(deg_min):
    min = (deg_min % 100) / 60
    deg = (deg_min // 100)
    deg_dec = min + deg
    return deg_dec


