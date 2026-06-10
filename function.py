import numpy as np
import pandas as pd

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

    df["time"] = df["time"].astype(float)
    df["latitude"] = df["latitude"].astype(float)
    df["longitude"] = df["longitude"].astype(float)
    df["quality"] = df["quality"].astype(float)
    df["nsats"] = df["nsats"].astype(float)
    df["hdop"] = df["hdop"].astype(float)
    df["geoidheight"] = df["geoidheight"].astype(float)
    df["geoidundulation"] = df["geoidundulation"].astype(float)

    return df

def calc_height(df):
    df['height'] = df['geoidheight'] + df['geoidundulation']
    return df

def calc_deg(deg_min):
    minutes = (deg_min % 100) / 60
    deg = (deg_min // 100)
    deg_dec = minutes + deg

    return np.radians(deg_dec)

def df_deg_dec(df):
    df['latitude'] = calc_deg(df['latitude'])
    df['longitude'] = calc_deg(df['longitude'])
    return df

def philamh_to_xyz(phi, lam, h, a, b):
  """
          Function to convert geographic coordinates (phi, lambda, h) to ECEF (X,Y,Z) coordinates.
          Args:
            phi: phi coordinate of point
            lam: lambda coordinate of point
            h: ellipsoidal height of point
            a: major semi axis of ellipsoid
            b: minor semi axis of ellipsoid

          Returns: numpy array with shape (X, Y, Z)
      """
  c = (a**2)/b
  e_strich_sq = (a**2 - b**2)/b**2
  V = np.sqrt(1 + e_strich_sq*np.cos(phi)**2)

  X = (c/V + h) * np.cos(phi) * np.cos(lam)
  Y = (c/V + h) * np.cos(phi) * np.sin(lam)
  Z = (b/V + h) * np.sin(phi)
  xyz_array = np.array([X, Y, Z])

  return xyz_array

def df_add_xyz(df, xyz_array):
    df['X'] = xyz_array[0, :]
    df['Y'] = xyz_array[1, :]
    df['Z'] = xyz_array[2, :]
    return df


def calc_mean_philam(df):
    phi_mean = df['latitude'].mean()
    lam_mean = df['longitude'].mean()
    return phi_mean, lam_mean

def calc_mean_xyz(df):
    x_mean = df['X'].mean()
    y_mean = df['Y'].mean()
    z_mean = df['Z'].mean()
    xyz_mean = np.array([x_mean, y_mean, z_mean])
    return xyz_mean


def rot_ned_to_ecef(phi, lam):
    """Rotation matrix from gnss ws2025/26"""
    return np.array([
        [-np.sin(phi) * np.cos(lam), -np.sin(lam), -np.cos(phi) * np.cos(lam)],
        [-np.sin(phi) * np.sin(lam), np.cos(lam), -np.cos(phi) * np.sin(lam)],
        [np.cos(phi), 0.0, -np.sin(phi)]
    ])




###########################################################

def ecef_to_neu(df):
    phi_mean, lam_mean = calc_mean_philam(df)

    phi_mean = np.radians(phi_mean)
    lam_mean = np.radians(lam_mean)

    xyz_mean = calc_mean_xyz(df)
    xyz = df[["X", "Y", "Z"]].to_numpy()
    dxyz = xyz - xyz_mean

    R_ecef_to_ned = rot_ned_to_ecef(phi_mean, lam_mean).T
    ned = dxyz @ R_ecef_to_ned.T

    north = ned[:, 0]
    east = ned[:, 1]
    up = -ned[:, 2]   # D -> U

    return north, east, up
