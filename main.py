import function as f

a_wgs84 = 6378137.00000  # m
b_wgs84 = 6356752.31425  # m

day1 = '15.11.2005'
day2 = '16.11.2005'

zx1_15_path = "data/ZX1_2005-11-15_1105-1135_GGA-only.TXT"
zx1_16_path = "data/ZX1_2005-11-16_1100-1130_GGA-only.TXT"
zx2_15_path = "data/ZX2_2005-11-15_1105-1135_GGA-only.TXT"
zx2_16_path = "data/ZX2_2005-11-16_1100-1130_GGA-only.TXT"


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


#######################################################
df_zx1_15["north"], df_zx1_15["east"], df_zx1_15["up"] = f.ecef_to_neu(df_zx1_15)
df_zx1_16["north"], df_zx1_16["east"], df_zx1_16["up"] = f.ecef_to_neu(df_zx1_16)
df_zx2_15["north"], df_zx2_15["east"], df_zx2_15["up"] = f.ecef_to_neu(df_zx2_15)
df_zx2_16["north"], df_zx2_16["east"], df_zx2_16["up"] = f.ecef_to_neu(df_zx2_16)

#f.plot_neu_scatter(df_zx1_15, "Antenna 1, 15.11.2005")
#f.plot_neu_scatter(df_zx1_16, "Antenna 1, 16.11.2005")
#f.plot_neu_scatter(df_zx2_15, "Antenna 2, 15.11.2005")
#f.plot_neu_scatter(df_zx2_16, "Antenna 2, 16.11.2005")


print(df_zx1_15[["north", "east", "up"]].mean())
print(df_zx1_15[["north", "east", "up"]].std())

zx1_15_xyz_mean = f.calc_mean_xyz(df_zx1_15)
zx1_16_xyz_mean = f.calc_mean_xyz(df_zx1_16)
zx2_15_xyz_mean = f.calc_mean_xyz(df_zx2_15)
zx2_16_xyz_mean = f.calc_mean_xyz(df_zx2_16)
print(f"mean ecef-coordinates (antenna 1, 15.11.2005)", zx1_15_xyz_mean)
print(f"mean ecef-coordinates (antenna 1, 16.11.2005)", zx1_16_xyz_mean)
print(f"mean ecef-coordinates (antenna 2, 15.11.2005)", zx2_15_xyz_mean)
print(f"mean ecef-coordinates (antenna 2, 16.11.2005)", zx2_16_xyz_mean)


df_dist_zx1_zx2_15 = f.geometric_distance(df_zx1_15, df_zx2_15)
df_dist_zx1_zx2_16 = f.geometric_distance(df_zx1_16, df_zx2_16)
df_dist_zx1_15_zx1_16 = f.geometric_distance(df_zx1_15, df_zx1_16)
df_dist_zx2_15_zx2_16 = f.geometric_distance(df_zx2_15, df_zx2_16)
print(df_dist_zx1_15_zx1_16.to_string())

f.plot_distance_timeseries(df_dist_zx1_zx2_15, "Antenna 1 & 2, 15.11.2005")
f.plot_distance_timeseries(df_dist_zx1_zx2_16, "Antenna 1 & 2, 16.11.2005")
f.plot_distance_timeseries(df_dist_zx1_15_zx1_16, "Antenna 1 & 1, 15.-16.11.2005")
f.plot_distance_timeseries(df_dist_zx2_15_zx2_16, "Antenna 1 & 1 15.-16.11.2005")
'''
plt.figure(figsize=(6, 6))
plt.scatter(df_zx1_15["east"], df_zx1_15["north"], s=5)
plt.xlabel("East [m]")
plt.ylabel("North [m]")
plt.title("Horizontal scatter in local NEU system")
plt.axis("equal")
plt.grid(True)
plt.show()'''
print(df_zx1_15)
f.plot_neu_diff_timeline(df_zx1_15, 'Receiver 1' , day1)
f.plot_neu_diff_timeline(df_zx1_16, 'Receiver 1' ,day2)
f.plot_neu_diff_timeline(df_zx2_15, 'Receiver 2' ,day1)
f.plot_neu_diff_timeline(df_zx2_16, 'Receiver 2' ,day2)

f.plot_hdop_timeline(df_zx1_15, 'Receiver 1' , day1)
f.plot_hdop_timeline(df_zx1_16, 'Receiver 1' ,day2)
f.plot_hdop_timeline(df_zx2_15, 'Receiver 2' ,day1)
f.plot_hdop_timeline(df_zx2_16, 'Receiver 2' ,day2)

#12:06