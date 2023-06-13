import numpy as np
from math import *
from argparse import ArgumentParser

o = object()

class Transformation:
    def __init__(self, model: str = "grs80"):
        '''       
        Wykorzystywane Parametry elipsoid
        a = dłuższa półos (promień rownikowy)
        b = krótsza półos (promień południkowy)
        flat = spłaszczenie
        e2 = mimosród^2
 
        Inicjuje obiekt elipsoidy na podstawie wybranego modelu.
        Dostępne modele to: WGS84, GRS80 i Mars.
        
        Argumenty:
        model (str): Łańcuch znaków określający model elipsoidy. 
                  Domyślnie ustawione na 'wgs84'.
        '''
        if model == "wgs84":
            self.a = 6378137 
            self.b = 6356752.31424518
        elif model == "grs80":
            self.a = 6378137
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a        
        self.e = sqrt(2 * self.flat - self.flat ** 2) 
        self.e2 = (2 * self.flat - self.flat ** 2)      
        

    def hirvonen(self,X,Y,Z):
        '''
        Zmiana współrzędnych geocentrycznych na współrzędne geodezyjne (Algorytm Hirvonena)

        Parameters
        ----------
        X : współrzędna geocentryczna [m]
        Y : współrzędna geocentryczna [m]
        Z : współrzędna geocentryczna [m]

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
        przeliczenie wszpolrzednych geodezyjnych 
        na wspolrzedne ukladu 2000

        Parameters
        ----------
        x : współrzędna geocentryczna [m]
        y : współrzędna geocentryczna [m]
        z : współrzędna geocentryczna [m]

        Returns
        -------
        x_00 : współrzędna geocentryczna w układzie PL - 2000
        y_00 : współrzędna geocentryczna w układzie PL - 2000

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
        x2000 = xgk * 0.999923
        y2000 = ygk * 0.999923 + strefa * 1000000 + 500000
        return (x2000, y2000)
    
    def x_y_1992(self, x, y, z):
        """   
        przeliczenie wszpolrzednych geodezyjnych 
        na wspolrzedne ukladu 1992
            
        Parameters
        -------
        x : wspolrzedna geocentryczna [m]
        y : wspolrzedna geocentryczna [m]
        z : wspolrzedna geocentryczna [m]
        a : dluzsza polos elipsoidy [m]
        e2: mimosrod elipsoidy [niemianowana]
            
        Returns
        -------
        x92 : wspolrzedna w ukladzie PL-1992 [m]
        y92 : wspolrzedna w ukladzie PL-1992 [m]  
        
        """ 
        a = self.a
        e2 = self.e2
        fi = self.hirvonen(x,y,z)[0]
        l = self.hirvonen(x,y,z)[1]
        m = 0.9993
        
        N = self.a/sqrt(1-self.e2*sin(fi)**2)
        e2p = self.e2/(1-self.e2)
        t = tan(fi)
        n2 = e2p * cos(fi)**2
        l0 = radians(19)
        lam = l - l0
        
        A0 = 1 - (self.e2/4) - ((3*(self.e2**2))/64) - ((5*(self.e2**3))/256)
        A2 = (3/8) * (self.e2 + ((self.e2**2)/4) + ((15 * (self.e2**3))/128))
        A4 = (15/256) * (self.e2**2 + ((3*(self.e2**3))/4))
        A6 = (35 * (self.e2**3))/3072
        
        sig = self.a * ((A0*fi) - (A2*sin(2*fi)) + (A4*sin(4*fi)) - (A6*sin(6*fi)))
        x = sig + ((lam**2)/2) * (N*sin(fi)*cos(fi)) * (1 + ((lam**2)/12) * ((cos(fi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((lam**4)/360) * ((cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = lam * (N * cos(fi)) * (1 + ((((lam**2)/6) * (cos(fi))**2) * (1-(t**2) + n2)) +  (((lam**4)/(120)) * (cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
            
        x1992 = x * m - 5300000
        y1992 = y * m + 500000
            
        return(x1992, y1992)

    def flh2XYZ(self, f, l, h):
        """
        Odwrotny algorytm Hirvonena - algorytm transformacji współrzędnych geodezyjnych 
        długość szerokość i wysokośc elipsoidalna(fi, lambda, h) na współrzędne ortokartezjańskie  (X, Y, Z).

        Parametery
        ----------
        f, l, h : FLOAT
            [dec_degree] współrzędne geodezyjne, 

        Returns
        -------
        X, Y, Z : FLOAT
            [metry] współrzędne ortokartezjańskie
        """
        f = radians(f)
        l = radians(l)
        N = self.a / sqrt(1 - self.e2 * sin(f)**2)
        X = (N + h) * cos(f) * cos(l)
        Y = (N + h) * cos(f) * sin(l)
        Z = (N * (1 - self.e2) + h) * sin(f)
        return X, Y, Z
    
    def Rneu(self, f, l):
        """
        Macierz R w transformacji współrzędnych XYZ na NEU jest macierzą rotacji, która pozwala przeliczyć 
        współrzędne z układu kartezjańskiego na współrzędne związanego z Ziemią układu współrzędnych geodezyjnych NEU.
        Wykorzystujemy bibliotekę numpy.
        
        Parametery
        ----------
        fi, lam: FLOAT
            [dec_degree] współrzędne fi, lambda w układzie geodezyjnym, 
       
        Returns
        -------
        R : array
            [niemianowane] macierz rotacji
        """
        R = np.array([[-np.sin(f) * np.cos(l), -np.sin(l), np.cos(f) * np.cos(l)],
                      [-np.sin(f) * np.sin(l), np.cos(l), np.cos(f) * np.sin(l)],
                      [np.cos(f), 0, np.sin(f)]])
        return(R)
    
    def XYZ_neu(self, xa, ya, za, xb, yb, zb):
        """
        Transformacja XYZ -> NEU - algorytm transformacji współrzędnych wektora pomiędzy dwoma punktami w układzie współrzędnych 
        ortokartezjańskich (X, Y, Z) na współrzędne wektora pomiędzy dwoma punktami w układzie NEU: North, East, Up (N, E, U). 
        Wykorzystujemy bibliotekę numpy.

        Parametery
        ----------
        Xa, Ya, Za : FLOAT
            [metry] współrzędne w układzie orto-kartezjańskim, 
        
        Xb, Yb, Zb : FLOAT
            [metry] współrzędne punktu referencyjnego w układzie orto-kartezjańskim, 

        Returns
        -------
        N, E, U : FLOAT
            [metry] współrzędne w układzie NEU
        """
        
        a = self.a
        e2 = self.e2
        dxyz = np.array([xb, yb, zb]) - np.array([xa, ya, za])
        f = self.hirvonen(xa, ya, za)[0]
        l = self.hirvonen(xa, ya, za)[1]
        R = self.Rneu(f, l)
        dneu = -np.linalg.solve(R, dxyz)
        N = dneu[0]
        E = dneu[1]
        U = dneu[2]

        return N, E, U

    


    def file_open(self, name):
        """
        Wczytanie pliku .txt i wyodrębnienie podanych w nim współrzędnych 
        za pomocą pętli for. Odczytane dane dodajemy do list.
        """
        
        X = []
        Y = []
        Z = []
        with open(name, 'r') as plik:
            lines = plik.readlines()
            t = 0
            for i in lines:
                x = i.split(',')
                X.append(float(x[0]))
                Y.append(float(x[1]))
                Z.append(float(x[2]))
            
        return X, Y, Z

    def file_save92(self, X, Y, Z, file_out):
        """
        Funkcja ta pozwala zapisać współrzędne otrzymane po transformacji XYZ -> PL-1992 
        do pliku wyjsciowego file_out w formacie .txt
        """
        
        X92 = []
        Y92 = []
        for a, b, c in zip(X, Y, Z):
            x92, y92 = geo.x_y_1992(a, b, c)
            X92.append(x92)
            Y92.append(y92)
            
        plik=open(file_out,"w")
        plik.write(f'# PL-1992---------------------------------------------- \n')
        plik.write(f'  X[m]         Y[m] \n')
        plik.write(f'# ----------------------------------------------------- \n')
        for a,b in zip(X92,Y92):
            a = f'{a:7.3f}'
            b = f'{b:7.3f}'
            plik.write(f'{a},   {b} \n')
        plik.close()
            
    def file_save00(self, X, Y, Z, file_out):
        """
        Funkcja ta pozwala zapisać współrzędne otrzymane po transformacji XYZ -> PL-2000 
        do pliku wyjsciowego file_out w formacie .txt
        """
        
        X00 = []
        Y00 = []
        for a, b, c in zip(X, Y, Z):
            x00, y00 = geo.x_y_2000(a, b, c)
            X00.append(x00)
            Y00.append(y00)
            
        plik=open(file_out,"w")
        plik.write(f'# PL-2000---------------------------------------------- \n')
        plik.write(f'  X[m]         Y[m] \n')
        plik.write(f'# ----------------------------------------------------- \n')
        for a,b in zip(X00,Y00):
            a = f'{a:7.3f}'
            b = f'{b:7.3f}'
            plik.write(f'{a},   {b} \n')
        plik.close()

            
    def file_saveFLH(self, X, Y, Z, file_out):
        """
        Funkcja ta pozwala zapisać współrzędne otrzymane po transformacji XYZ -> BLH 
        do pliku wyjsciowego file_out w formacie .txt
        """
        
        F = []
        L = []
        H = []
        for a, b, c in zip(X, Y, Z):
            f, l, h = geo.hirvonen(a, b, c)
            F.append(degrees(f))
            L.append(degrees(l))
            H.append(h)
            
        plik=open(file_out,"w")
        plik.write(f'  B[d]         L[d]         H[m] \n')
        plik.write(f'# ----------------------------------------------------- \n')
        for a,b,c in zip(F,L,H):
            a = f'{a:7.4f}'
            b = f'{b:7.4f}'
            c = f'{c:7.3f}'
            plik.write(f'{a},      {b},      {c} \n')
        plik.close()
            
    def file_saveNEU(self, X, Y, Z, xref, yref, zref, file_out):
        """
        Funkcja ta pozwala zapisać współrzędne otrzymane po transformacji XYZ -> NEU 
        do pliku wyjsciowego file_out w formacie .txt
        """
        
        N = []
        E = []
        U = []
        for a, b, c in zip(X, Y, Z):
            n, e, u = geo.XYZ_neu(a, b, c, xref, yref, zref)
            N.append(n)
            E.append(e)
            U.append(u)

        plik=open(file_out,"w")
        plik.write(f'  N[m]         E[m]         U[m] \n')
        plik.write(f'# ----------------------------------------------------- \n')

        for a,b,c in zip(N,E,U):
            a = f'{a:7.3f}'
            b = f'{b:7.3f}'
            c = f'{c:7.3f}'
            plik.write(f'{a},   {b},      {c} \n')
        plik.close()
        
    def file_saveXYZ(self, B, L, H, file_out):
        """
        Funkcja ta pozwala zapisać współrzędne otrzymane po transformacji BLH -> XYZ 
        do pliku wyjsciowego file_out w formacie .txt
        """
        
        X = []
        Y = []
        Z = []
        for a, b, c in zip(B, L, H):
            x, y, z = geo.flh2XYZ(a, b, c)
            X.append(x)
            Y.append(y)
            Z.append(z)
            
        plik=open(file_out,"w")
        plik.write(f'  X[m]         Y[m]         Z[m] \n')
        plik.write(f'# ----------------------------------------------------- \n')
        for a,b,c in zip(X,Y,Z):
            a = f'{a:7.3f}'
            b = f'{b:7.3f}'
            c = f'{c:7.3f}'
            plik.write(f'{a},      {b},      {c} \n')
        plik.close()

'''
    
X = []
Y = []
Z = []
F = []
L = []
H = []
X_92 = []
Y_92 = []
X_00 = []
Y_00 = []
N = []
E = []
U = []


            
with open('wsp_inp.txt', 'r') as plik:
    lines = plik.readlines()
    t = 0
    for i in lines:
        t = t + 1
        if t > 1:
            x = i.split(',')
            X.append(float(x[0]))
            Y.append(float(x[1]))
            Z.append(float(x[2]))
            
            #print(X)
   '''     
        
if __name__ == "__main__":
    #tworze obiekt
    #geo = Transformation(model = "wgs84")
    
    parser = ArgumentParser() 
    parser.add_argument('--dane', type = str, choices = ['XYZ', 'BLH'], default = 'XYZ', help = 'Typ wprowadzanych współrzędnych (BLH lub XYZ), domyslnie: XYZ' )
    parser.add_argument('-x_ref', type = float, help = 'Współrzędna X punktu referencyjnego', default = 100.00)
    parser.add_argument('-y_ref', type = float, help = 'Współrzędna Y punktu referencyjnego', default = 100.00)
    parser.add_argument('-z_ref', type = float, help = 'Współrzędna Z punktu referencyjnego', default = 100.00)
    parser.add_argument('--model', type = str, choices = ['wgs84', 'grs80', 'krasowski'], default = 'wgs84', help = 'Model elipsoidy (wgs84, grs80 lub krasowski), domyslnie: wgs84')
    parser.add_argument('--uklad', type = str, choices = ['PL-1992', 'PL-2000', 'BLH', 'NEU'], default = 'BLH', help= 'System współrzędnych (PL-1992, PL-2000, BLH, NEU), domyslnie: BLH')
    parser.add_argument('-txt', type = str, help = 'Nazwa pliku wejsciowego z rozszerzeniem .txt, zawierającego współrzędne do transformacji')
    parser.add_argument('-txt_out', type = str, help = 'Nazwa pliku wyjsciowego z rozszerzeniem .txt zawierającego przetransformowane współrzędne')
    args = parser.parse_args()
    geo = Transformation(model = args.model)
    if args.dane == 'XYZ':
            if args.uklad == 'PL-1992':
                X, Y, Z = geo.file_open(args.txt)
                geo.file_save92(X, Y, Z, args.txt_out)
            elif args.uklad == 'PL-2000':
                X, Y, Z = geo.file_open(args.txt)
                geo.file_save00(X, Y, Z, args.txt_out)
            elif args.uklad == 'BLH':
                X, Y, Z = geo.file_open(args.txt)
                geo.file_saveFLH(X, Y, Z, args.txt_out)
            elif args.uklad == 'NEU':
                X, Y, Z = geo.file_open(args.txt)
                geo.file_saveNEU(X, Y, Z, args.x_ref, args.y_ref, args.z_ref, args.txt_out)
    elif args.dane == 'BLH':
            B, L, H = geo.file_open(args.txt)
            geo.file_saveXYZ(B, L, H, args.txt_out)
            
    #wsp geocentryczne
#    x = 100; y = 120; z = 0
#   x2000,y2000 = geo.x_y_2000(x,y,z)
#    x1992,y1992 = geo.x_y_1992(x, y, z)
#    f1, l1, h = geo.hirvonen(x,y,z)
#   # N, E, U = geo.XYZ_neu(x, y, z )
#    f = f1 * 180 / pi
#    l = l1 * 180 / pi
#    print(x2000,y2000)
#    print('')
#    print(x1992,y1992)
#    print('')
#    print(f, l, h)
#    print('')
   #print(N, E, U)
#    print('')
    

    

















    '''
    X = []
    Y = []
    Z = []
    F = []
    L = []
    H = []
    X_92 = []
    Y_92 = []
    X_00 = []
    Y_00 = []
    N = []
    E = []
    U = []


                
    with open('wsp_inp.txt', 'r') as plik:
        lines = plik.readlines()
        t = 0
        for i in lines:
            t = t + 1
            if t > 1:
                x = i.split(',')
                X.append(float(x[0]))
                Y.append(float(x[1]))
                Z.append(float(x[2]))
                
                
    for A,B,C in zip(X,Y,Z):
        f, l, h = geo.hirvonen(A, B, C)
        F.append(degrees(f))
        L.append(degrees(l))
        H.append(h)
        x1992, y1992 = geo.x_y_1992(A, B, C)
        X_92.append(x1992)
        Y_92.append(y1992)
        x2000, y2000 = geo.x_y_2000(A, B, C)
        X_00.append(x2000)
        Y_00.append(y2000)
        n, e, u = geo.XYZ_neu(A, B, C)
        N.append(n)
        E.append(e)
        U.append(u)
    





    plik=open(args.t,"w")
    plik.write(f'Współrzędne flh, PL_1992, PL_2000, NEU stacji permanentnej GNSS \n')
    plik.write(f'Obserwatorium Astronomiczno-Geodezyjne w Józefosławiu \n')
    plik.write(f'# ************************************* \n')
    plik.write(f'# fLh **********************************\n')
    plik.write(f'  f[d]         l[d]         h[m] \n')
    plik.write(f'# ************************************* \n')
    for A,B,C in zip(F,L,H):
        A = f'{A:7.4f}'
        B = f'{B:7.4f}'
        C = f'{C:7.4f}'
        plik.write(f'{A},      {B},      {C} \n')
        
      
    plik.write(f'# ************************************* \n')
    plik.write(f'# PL_2000 ************************************* \n')
    plik.write(f'  X[m]         Y[m] \n')
    plik.write(f'# ************************************* \n')
    for A,B in zip(X_00,Y_00):
        A = f'{A:7.3f}'
        B = f'{B:7.3f}'
        plik.write(f'{A},   {B} \n')
        
    plik.write(f'# ************************************* \n')
    plik.write(f'# PL_1992 ************************************* \n')
    plik.write(f'  X[m]         Y[m] \n')
    plik.write(f'# ************************************* \n')
    for A,B in zip(X_92,Y_92):
        A = f'{A:7.3f}'
        B = f'{B:7.3f}'
        plik.write(f'{A},   {B} \n')
    
    plik.write(f'# ************************************* \n')
    plik.write(f'# NEU ************************************* \n')
    plik.write(f'  N[m]         E[m]         U[m] \n')
    plik.write(f'# ************************************* \n')
    
    for A,B,C in zip(N,E,U):
        A = f'{A:7.3f}'
        B = f'{B:7.3f}'
        C = f'{C:7.3f}'
        plik.write(f'{A},   {B},      {C} \n')
    plik.close()
     


'''



    

    
    
    
    
    
    
    
    
    

