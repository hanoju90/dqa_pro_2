import function as f

a_wgs84 = 6378137.00000  # m
b_wgs84 = 6356752.31425  # m

day1 = '15.11.2005'
day2 = '16.11.2005'

zx1_15_path = "data/ZX1_2005-11-15_1105-1135_GGA-only.TXT"
zx1_16_path = "data/ZX1_2005-11-16_1100-1130_GGA-only.TXT"
zx2_15_path = "data/ZX2_2005-11-15_1105-1135_GGA-only.TXT"
zx2_16_path = "data/ZX2_2005-11-16_1100-1130_GGA-only.TXT"


# data preprocessing
zx1_15_geo = f.read_nmea(zx1_15_path)
zx1_16_geo = f.read_nmea(zx1_16_path)
zx2_15_geo = f.read_nmea(zx2_15_path)
zx2_16_geo = f.read_nmea(zx2_16_path)

zx1_15_geo = f.convert_time(zx1_15_geo)
zx1_16_geo = f.convert_time(zx1_16_geo)
zx2_15_geo = f.convert_time(zx2_15_geo)
zx2_16_geo = f.convert_time(zx2_16_geo)

zx1_15_geo = f.calc_height(zx1_15_geo)
zx1_16_geo = f.calc_height(zx1_16_geo)
zx2_15_geo = f.calc_height(zx2_15_geo)
zx2_16_geo = f.calc_height(zx2_16_geo)

zx1_15_geo = f.df_deg_dec(zx1_15_geo)
zx1_16_geo = f.df_deg_dec(zx1_16_geo)
zx2_15_geo = f.df_deg_dec(zx2_15_geo)
zx2_16_geo = f.df_deg_dec(zx2_16_geo)

zx1_15_xyz = f.philamh_to_xyz(zx1_15_geo['latitude'], zx1_15_geo['longitude'], zx1_15_geo['height'], a_wgs84, b_wgs84)
zx1_16_xyz = f.philamh_to_xyz(zx1_16_geo['latitude'], zx1_16_geo['longitude'], zx1_16_geo['height'], a_wgs84, b_wgs84)
zx2_15_xyz = f.philamh_to_xyz(zx2_15_geo['latitude'], zx2_15_geo['longitude'], zx2_15_geo['height'], a_wgs84, b_wgs84)
zx2_16_xyz = f.philamh_to_xyz(zx2_16_geo['latitude'], zx2_16_geo['longitude'], zx2_16_geo['height'], a_wgs84, b_wgs84)

df_zx1_15 = f.df_add_xyz(zx1_15_geo, zx1_15_xyz)
df_zx1_16 = f.df_add_xyz(zx1_16_geo, zx1_16_xyz)
df_zx2_15 = f.df_add_xyz(zx2_15_geo, zx2_15_xyz)
df_zx2_16 = f.df_add_xyz(zx2_16_geo, zx2_16_xyz)


# 1 a)
df_zx1_15["north"], df_zx1_15["east"], df_zx1_15["up"] = f.ecef_to_neu(df_zx1_15)
df_zx1_16["north"], df_zx1_16["east"], df_zx1_16["up"] = f.ecef_to_neu(df_zx1_16)
df_zx2_15["north"], df_zx2_15["east"], df_zx2_15["up"] = f.ecef_to_neu(df_zx2_15)
df_zx2_16["north"], df_zx2_16["east"], df_zx2_16["up"] = f.ecef_to_neu(df_zx2_16)

f.plot_neu_scatter(df_zx1_15, "15.11.2005 - Receiver 1")
f.plot_neu_scatter(df_zx1_16, "16.11.2005 - Receiver 1")
f.plot_neu_scatter(df_zx2_15, "15.11.2005 - Receiver 2")
f.plot_neu_scatter(df_zx2_16, "16.11.2005 - Receiver 2")


 # 1 b) i)
f.plot_neu_diff_timeline(df_zx1_15, 'Receiver 1' , day1)
f.plot_neu_diff_timeline(df_zx1_16, 'Receiver 1' ,day2)
f.plot_neu_diff_timeline(df_zx2_15, 'Receiver 2' ,day1)
f.plot_neu_diff_timeline(df_zx2_16, 'Receiver 2' ,day2)

# 1 b) ii)
f.plot_hdop_timeline(df_zx1_15, 'Receiver 1' , day1)
f.plot_hdop_timeline(df_zx1_16, 'Receiver 1' ,day2)
f.plot_hdop_timeline(df_zx2_15, 'Receiver 2' ,day1)
f.plot_hdop_timeline(df_zx2_16, 'Receiver 2' ,day2)


# 2 b)
zx1_15_xyz_mean = f.calc_mean_xyz(df_zx1_15)
zx1_16_xyz_mean = f.calc_mean_xyz(df_zx1_16)
zx2_15_xyz_mean = f.calc_mean_xyz(df_zx2_15)
zx2_16_xyz_mean = f.calc_mean_xyz(df_zx2_16)
print(f"mean ecef-coordinates (15.11.2005 - receiver 1)", zx1_15_xyz_mean)
print(f"mean ecef-coordinates (16.11.2005 - receiver 1)", zx1_16_xyz_mean)
print(f"mean ecef-coordinates (15.11.2005 - receiver 2)", zx2_15_xyz_mean)
print(f"mean ecef-coordinates (16.11.2005 - receiver 2)", zx2_16_xyz_mean)


# 2 c)
df_dist_zx1_zx2_15 = f.geometric_distance(df_zx1_15, df_zx2_15)
df_dist_zx1_zx2_16 = f.geometric_distance(df_zx1_16, df_zx2_16)
df_dist_zx1_15_zx1_16 = f.geometric_distance(df_zx1_15, df_zx1_16)
df_dist_zx2_15_zx2_16 = f.geometric_distance(df_zx2_15, df_zx2_16)

df_dist_zx1_zx2_15 = f.convert_time(df_dist_zx1_zx2_15)
df_dist_zx1_zx2_16 = f.convert_time(df_dist_zx1_zx2_16)
df_dist_zx1_15_zx1_16 = f.convert_time(df_dist_zx1_15_zx1_16)
df_dist_zx2_15_zx2_16 = f.convert_time(df_dist_zx2_15_zx2_16)

f.plot_distance_timeseries(df_dist_zx1_zx2_15, "15.11.2005 - Receiver 1 & 2")
f.plot_distance_timeseries(df_dist_zx1_zx2_16, "16.11.2005 - Receiver 1 & 2")
f.plot_distance_timeseries(df_dist_zx1_15_zx1_16, "15.-16.11.2005 - Receiver 1")
f.plot_distance_timeseries(df_dist_zx2_15_zx2_16, "15.-16.11.2005 - Receiver 2")

lags_zx1_zx2_15, acf_zx1_zx2_15 = f.calc_sample_autocorrelation(df_dist_zx1_zx2_15["distance"], max_lag=400)
lags_zx1_zx2_16, acf_zx1_zx2_16 = f.calc_sample_autocorrelation(df_dist_zx1_zx2_16["distance"], max_lag=400)
lags_zx1_15_zx1_16, acf_zx1_15_zx1_16 = f.calc_sample_autocorrelation(df_dist_zx1_15_zx1_16["distance"], max_lag=400)
lags_zx2_15_zx2_16, acf_zx2_15_zx2_16 = f.calc_sample_autocorrelation(df_dist_zx2_15_zx2_16["distance"], max_lag=400)

f.plot_sample_autocorrelation(lags_zx1_zx2_15, acf_zx1_zx2_15, "15.11.2005 - Receiver 1 & 2")
f.plot_sample_autocorrelation(lags_zx1_zx2_16, acf_zx1_zx2_16, "16.11.2005 - Receiver 1 & 2")
f.plot_sample_autocorrelation(lags_zx1_15_zx1_16, acf_zx1_15_zx1_16, "15.-16.11.2005 - Receiver 1")
f.plot_sample_autocorrelation(lags_zx2_15_zx2_16, acf_zx2_15_zx2_16, "15.-16.11.2005 - Receiver 2")