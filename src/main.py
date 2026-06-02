import numpy as np


def calc_deg(deg_min):
    min = (deg_min % 100) / 60
    deg = (deg_min // 100)
    deg_dec = min + deg
    return deg_dec

lat = 4703.85386
lon = 01527.195522

dec_lat = calc_deg(lat)
dec_lon = calc_deg(lon)

print(dec_lat, dec_lon)

