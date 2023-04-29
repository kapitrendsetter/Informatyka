<<<<<<< HEAD
from math import *

o = object()

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2

    def hirvonen(self,X,Y,Z):
        '''
        

        Parameters
        ----------
        X : TYPE
            DESCRIPTION.
        Y : TYPE
            DESCRIPTION.
        Z : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        p = sqrt(X**2 + Y**2)
        f = atan(Z / (p * (1 - self.e2)))
        while True:
            N = self.a / sqrt(1 - self.e2 * sin(f)**2)
            h = (p/cos(f)) - N
            fp = f
            f = atan(Z/(p*(1 - self.e2 * N/ (N +h))))
            if abs(fp - f) < (0.000001/206265):
                break
        l = atan2(Y , X)
        return(f,l,h)
    
    
    
    

    
    
=======
import argparse
import numpy as np


# inicjalizacja parsera
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="nazwa pliku txt z danymi wejściowymi")

# parsowanie argumentów
args = parser.parse_args()

# otwieranie pliku i odczytanie zawartości
with open(args.filename, 'r') as file:
    for line in file:  

        # przetwarzanie danych z pliku linia po linii
        # i wywołanie odpowiedniej funkcji transformacji


class CoordinatesTransformations:
    def __init__(self, ellipsoid='GRS80'):
        self.a, self.f = self._get_ellipsoid_params(ellipsoid)

    def _get_ellipsoid_params(self, ellipsoid):
        ellipsoids = {
            'GRS80': (6378137.0, 1/298.257222101),
            'WGS84': (6378137.0, 1/298.257223563),
            'Krasowski': (6378245.0, 1/298.3)
        }
        if ellipsoid not in ellipsoids:
            raise ValueError(f'Unsupported ellipsoid: {ellipsoid}')
        return ellipsoids[ellipsoid]

    def _deg2rad(self, degrees):
        return degrees * np.pi / 180

    def _rad2deg(self, radians):
        return radians * 180 / np.pi

    def xyz2blh(self, x, y, z):
        p = np.sqrt(x**2 + y**2)
        theta = np.arctan2(z * self.a, p * self.a * (1 - self.f))
        phi = np.arctan2(y, x)
        N = self.a / np.sqrt(1 - self.f * (2 - self.f) * np.sin(theta)**2)
        h = p / np.cos(theta) - N
        return self._rad2deg(phi), self._rad2deg(theta), h

    def blh2xyz(self, phi, theta, h):
        N = self.a / np.sqrt(1 - self.f * (2 - self.f) * np.sin(self._deg2rad(theta))**2)
        x = (N + h) * np.cos(self._deg2rad(theta)) * np.cos(self._deg2rad(phi))
        y = (N + h) * np.cos(self._deg2rad(theta)) * np.sin(self._deg2rad(phi))
        z = (N * (1 - self.f)**2 + h) * np.sin(self._deg2rad(theta))
        return x, y, z

    def xyz2neup(self, x, y, z, lat0, lon0):
        phi0 = self._deg2rad(lat0)
        lambda0 = self._deg2rad(lon0)
        T = np.array([[-np.sin(lambda0), np.cos(lambda0), 0],
                      [-np.sin(phi0)*np.cos(lambda0), -np.sin(phi0)*np.sin(lambda0), np.cos(phi0)],
                      [np.cos(phi0)*np.cos(lambda0), np.cos]])
>>>>>>> d1afa99aef527b8991df6f1209c5173e01927ade
