
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
    
    def x_y_2000(self, X, Y, Z):
        '''
        

        Parameters
        ----------
        x : TYPE
            DESCRIPTION.
        y : TYPE
            DESCRIPTION.
        z : TYPE
            DESCRIPTION.

        Returns
        -------
        x_00 : TYPE
            DESCRIPTION.
        y_00 : TYPE
            DESCRIPTION.

        '''
        a = self.a
        e2 = self.e2
        fi = self.hirvonen(X,Y,Z)[0]
        lam = self.hirvonen(X,Y,Z)[1]
        if abs(degrees(lam) - 15) <= 1.5:
            l0_deg = 15
        elif abs(degrees(lam) - 18) < 1.5:
            l0_deg = 18
        elif abs(degrees(lam) - 21) <= 1.5:
            l0_deg = 21
        else:
            l0_deg = 24
        l0 = radians(l0_deg)
        a2 = a**2
        b2 = a2 * (1 - e2)
        e_2 = (a2 - b2)/b2
        dl = lam - l0
        dl2 = dl**2
        dl4 = dl**4
        t = tan(fi)
        t2 = t**2
        t4 = t**4
        n2 = e_2 * (cos(fi)**2)
        n4 = n2 ** 2
        N = self.a / sqrt(1 - self.e2 * sin(fi)**2)
        e4 = e2**2
        e6 = e2**3
        A0 = 1 - (e2/4) - ((3*e4)/64) - ((5*e6)/256)
        A2 = (3/8) * (e2 + e4/4 + (15*e6)/128)
        A4 = (15/256) * (e4 + (3*e6)/4)
        A6 = (35*e6)/3072
        sigma = a * ((A0 * fi) - A2 * sin(2*fi) + A4 * sin(4*fi) - A6 * sin(6*fi))
        xgk = sigma + ((dl**2)/2) * N * sin(fi) * cos(fi) * (1 + ((dl**2)/12)*(cos(fi)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (cos(fi)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
        ygk = dl * N * cos(fi) * (1 + (dl2/6) * (cos(fi)**2) * (1 - t2 + n2) + (dl4/120) * (cos(fi)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
        strefa = int(l0 * 180/pi)/3
        x_00 = xgk * 0.999923
        y_00 = ygk * 0.999923 + strefa * 1000000 + 500000
        return x_00, y_00
    
    def x_y_1992(self, X, Y, Z):
       '''
        

        Parameters
        ----------
        x : TYPE
            DESCRIPTION.
        y : TYPE
            DESCRIPTION.
        z : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        a = self.a
        e2 = self.e2
        fi= self.hirvonen(X,Y,Z)[0]
        lam = self.hirvoen(X,Y,Z)[1]
        m = 0.9993
        
        N = self.a/sqrt(1-self.e2*sin(fi)**2)
        e2p = self.e2/(1-self.e2)
        t = tan(fi)
        n2 = e2p * cos(fi)**2
        l0 = radians(19)
        lam1 = lam - l0
        
        A0 = 1 - (self.e2/4) - ((3*(self.e2**2))/64) - ((5*(self.e2**3))/256)
        A2 = (3/8) * (self.e2 + ((self.e2**2)/4) + ((15 * (self.e2**3))/128))
        A4 = (15/256) * (self.e2**2 + ((3*(self.e2**3))/4))
        A6 = (35 * (self.e2**3))/3072
        
        sig = self.a * ((A0*fi) - (A2*sin(2*fi)) + (A4*sin(4*fi)) - (A6*sin(6*fi)))
        x = sig + ((lam1**2)/2) * (N*sin(fi)*cos(fi)) * (1 + ((lam1**2)/12) * ((cos(fi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((lam1**4)/360) * ((cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = lam1 * (N * cos(fi)) * (1 + ((((lam1**2)/6) * (cos(fi))**2) * (1-(t**2) + n2)) +  (((lam1**4)/(120)) * (cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
            
        x1992 = x * m - 5300000
        y1992 = y * m + 500000
            
        return(x1992, y1992)
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

