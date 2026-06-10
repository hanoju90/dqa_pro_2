import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

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

def convert_time(df):
    df["time"] = pd.to_datetime(df["time"], format='%H%M%S')
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

    xyz_mean = calc_mean_xyz(df)
    xyz = df[["X", "Y", "Z"]].to_numpy()
    dxyz = xyz - xyz_mean

    R_ecef_to_ned = rot_ned_to_ecef(phi_mean, lam_mean).T
    ned = dxyz @ R_ecef_to_ned.T

    north = ned[:, 0]
    east = ned[:, 1]
    up = -ned[:, 2]   # D -> U

    return north, east, up


def plot_neu_scatter(df, title):
    plt.figure(figsize=(6, 6))

    plt.xlabel("East [m]")
    plt.ylabel("North [m]")

    plt.suptitle("Horizontal Scatter North-East-Up", fontsize=14)
    plt.title(title, fontsize=12)

    plt.scatter(df["east"], df["north"], s=50, edgecolors="black", linewidths=0.3, color="blue", label="Measurements")
    plt.scatter(df["east"].mean(), df["north"].mean(), s=15, color="red", label="Mean Position")

    plt.xticks(np.linspace(-2, 2, 9))
    plt.yticks(np.linspace(-2, 2, 9))

    plt.xlim(-2.4, 2.4)
    plt.ylim(-2.4, 2.4)

    plt.grid(True)
    plt.legend()
    plt.show()


def geometric_distance(df1, df2):
    df_pair = pd.merge(
        df1[["time", "X", "Y", "Z"]],
        df2[["time", "X", "Y", "Z"]],
        on="time",
        suffixes=("_1", "_2")
    )

    df_pair["distance"] = np.sqrt(
        (df_pair["X_2"] - df_pair["X_1"])**2 +
        (df_pair["Y_2"] - df_pair["Y_1"])**2 +
        (df_pair["Z_2"] - df_pair["Z_1"])**2
    )

    df_distance = df_pair[["time", "distance"]].copy()

    return df_distance

def plot_distance_timeseries(df, title):
    plt.figure(figsize=(10, 3))
    plt.plot(df["time"], df["distance"], color="blue")

    plt.xlabel("Time")
    plt.ylabel("Distance [m]")

    plt.title(f"Distance Time Series\n{title}", fontsize=12)

    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_neu_diff_timeline(df, receiver, day):
  fig, ax = plt.subplots()
  ax.plot(df['time'], df['north'],
                 label='North', color='red')
  ax.plot(df['time'], df['east'],
                 label='East', color='green')
  ax.plot(df['time'], df['up'],
                 label='Up', color='blue')
  ax.set_ylabel('Difference [m]')
  ax.set_xlabel('Time')
  ax.legend(loc='upper right')
  ax.grid(True)

  locator = mdates.AutoDateLocator()
  ax.xaxis.set_major_locator(locator)
  ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax.set_yticks(np.arange(0, 5))
  fig.suptitle(f'Difference North-East-Up', fontsize=14)
  plt.title(f'{day} - {receiver}', fontsize=12)
  plt.show()

def plot_hdop_timeline(df, receiver, day):
  fig, ax = plt.subplots()
  ax.plot(df['time'], df['hdop'], color='red', label=f'HDOP')
  ax.set_xlabel("Time")
  ax.set_ylabel("HDOP")

  locator = mdates.AutoDateLocator()
  ax.xaxis.set_major_locator(locator)
  ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  fig.suptitle(f"HDOP :", fontsize=14)
  plt.title(f'{day} - {receiver}', fontsize=12)
  ax.legend(loc='upper right')
  plt.legend()
  plt.grid()
  plt.show()

def export_df_as_csv(df, filename):
    df["time"] = df["time"].dt.strftime("%H:%M:%S")
    df.to_csv(f'{filename}.csv', sep=';', encoding='utf-8', index=False, header=True)