import numpy as np
import pandas as pd

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

    df["time"] = df["time"].astype(float)
    df["latitude"] = df["latitude"].astype(float)
    df["longitude"] = df["longitude"].astype(float)
    df["quality"] = df["quality"].astype(float)
    df["nsats"] = df["nsats"].astype(float)
    df["hdop"] = df["hdop"].astype(float)
    df["geoidheight"] = df["geoidheight"].astype(float)
    df["geoidundulation"] = df["geoidundulation"].astype(float)

    return df

zx1_15_geo = read_nmea(zx1_15_path)
zx1_16_geo = read_nmea(zx1_16_path)
zx2_15_geo = read_nmea(zx2_15_path)
zx2_16_geo = read_nmea(zx2_16_path)



def calc_deg(deg_min):
    min = (deg_min % 100) / 60
    deg = (deg_min // 100)
    deg_dec = min + deg
    return deg_dec

def geo_to_ecef(df):
    h = df.geoidheight + df.geoidundulation

    df.latitude = calc_deg(df.latitude)

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

def df_deg_dec(df, xyz_array):
    df['latitude'] = calc_deg(df['latitude'])
    df['longitude'] = calc_deg(df['longitude'])
    df['X'] = xyz_array[:, 0]
    df['Y'] = xyz_array[:, 1]
    df['Z'] = xyz_array[:, 2]
    return df

print(zx1_15_geo)

def create_rotation_matrix(lam, phi):
  '''
  Create rotation matrix to go from ECEF to N-E-D at a specific position
      Args:
          pos_xyz: ECEF coordinates of the position for which the rotation matrix should be calculated

      Returns: numpy array of the rotation matrix
  '''
  lam = np.squeeze(lam)
  phi = np.squeeze(phi)
  lam = np.radians(lam)
  phi = np.radians(phi)
  R_NED = np.array([[-np.sin(phi) * np.cos(lam), -np.sin(lam), np.cos(phi) * -np.cos(lam)],
                        [-np.sin(phi) * np.sin(lam), np.cos(lam), np.cos(phi) * -np.sin(lam)],
                        [np.cos(phi), 0, -np.sin(phi)]])
  return R_NED.T

def calc_height(df):
    df['height'] = df['geoidheight'] + df['geoidundulation']
    return df

xyz_array = philamh_to_xyz(zx2_15_geo['latitude'], zx2_15_geo['longitude'], zx2_15_geo['height'], )