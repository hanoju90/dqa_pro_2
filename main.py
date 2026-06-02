import numpy as np
import function as f

a_wgs84 = 6378137.00000  # m
b_wgs84 = 6356752.31425  # m

zx1_15_path = "data/ZX1_2005-11-15_1105-1135_GGA-only.TXT"
zx1_16_path = "data/ZX1_2005-11-16_1100-1130_GGA-only.TXT"
zx2_15_path = "data/ZX1_2005-11-15_1105-1135_GGA-only.TXT"
zx2_16_path = "data/ZX1_2005-11-16_1100-1130_GGA-only.TXT"

zx1_15_geo = f.read_nmea(zx1_15_path)
zx1_16_geo = f.read_nmea(zx1_16_path)
zx2_15_geo = f.read_nmea(zx2_15_path)
zx2_16_geo = f.read_nmea(zx2_16_path)

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

#mean phi lam coords
phi_mean_zx1_15, lam_mean_zx1_15 = f.calc_mean_philam(df_zx1_15)
phi_mean_zx1_16, lam_mean_zx1_16 = f.calc_mean_philam(df_zx1_16)
phi_mean_zx2_15, lam_mean_zx2_15 = f.calc_mean_philam(df_zx2_15)
phi_mean_zx2_16, lam_mean_zx2_16 = f.calc_mean_philam(df_zx2_16)

xyz_mean_zx1_15 = f.calc_mean_xyz(df_zx1_15)
xyz_mean_zx1_16 = f.calc_mean_xyz(df_zx1_16)
xyz_mean_zx2_15 = f.calc_mean_xyz(df_zx2_15)
xyz_mean_zx2_16 = f.calc_mean_xyz(df_zx2_16)

'''zx1_15_R = f.create_rotation_matrix(phi_mean_zx1_15, lam_mean_zx1_15)
zx1_16_R = f.create_rotation_matrix(phi_mean_zx1_16, lam_mean_zx1_16)
zx2_15_R = f.create_rotation_matrix(phi_mean_zx2_15, lam_mean_zx2_15)
zx2_16_R = f.create_rotation_matrix(phi_mean_zx2_16, lam_mean_zx2_16)

zx1_15_ll = zx1_15_R @  zx1_15_xyz
zx_16_ll =  zx1_16_R @ zx1_16_xyz
zx2_15_ll = zx2_15_R @ zx2_15_xyz
zx2_16_ll = zx2_16_R @ zx2_16_xyz'''

zx1_15_ll = f.dX_local_level(zx1_15_xyz, phi_mean_zx1_15, lam_mean_zx1_15, xyz_mean_zx1_15)
zx1_16_ll = f.dX_local_level(zx1_16_xyz, phi_mean_zx1_16, lam_mean_zx1_16, xyz_mean_zx1_16)
zx_2_15_ll = f.dX_local_level(zx2_15_xyz, phi_mean_zx2_15, lam_mean_zx2_15, xyz_mean_zx2_15)
zx2_16_ll = f.dX_local_level(zx2_16_xyz, phi_mean_zx2_16, lam_mean_zx2_16, xyz_mean_zx2_16)

print(zx1_15_ll.shape)




