# Projekt Transformacje
Skryp właściwy znajduje się w pliku script_12.py.
Skrypt potrzebny do uzywania biblioteki argparse znajduje sie w pliku proj2.py

# Q&A

# Do czego służy program i jaką funkcjonalność oferuje?  
Transformacje:
- XYZ -> BLH, 
- BLH -> XYZ 
- XYZ -> NEUp
- BL(GRS80, WGS84, ew. Krasowski) -> 2000 
- BL(GRS80, WGS84, ew. Krasowski) -> 1992

# Jakie elipsidy są obsługiwane?
- GRS80 
- WGS84

# Jakie wymagania trzeba spełnić, by program działał na danym komputerze?
- trzeba mieć pythona w wersji 10.9 lub nowszej 
- zainstalowaną bibliotekę argparse (która miewa problemy z odczytem)

# Dla jakiego systemu operacyjnego został napisany program?
- program działa zarówno na Windows 10 jak i macOS Ventura 13.0

# Instrukcja jak używać program na 2 sposoby: 

- Wczytanie pliku.txt w skrypcie python
plik.txt - plik ze współrzędnymi X, Y, Z 
Należy je wpisać w odpowiedniej kolejności X,Y,Z oddzieljąc przecinkami
NP:
88888.000,99999.000,100000,000
Plik ze współrzędnymi należy mieć w danym folderze w którym znajduje się skrypt programu
Żeby program wywołał wyniki nazwę pliku należy wpisac do skryptu w miejsce open('nazwa_pliku.txt','r')
Nastepnie pliku o nazwie wyniki.txt utworzy się w folderze w którym znajduję się skrypt 
Kolejność wyników:wsp. 1- geodezyjne, 2- PL-2000, 3- Pl-1992, 4- NEU
Plik wejściowy musi mieć zawsze nazwę "nazwa_pliku.txt" , aby użytkownik miał możliwość łatwej edycji nazwy pliku.

- Wykorzystanie biblioteki argparse w celu przeliczenia współrzędnych w wybranym odwzorowaniu:
Należy otworzyć command window ( w folderze w którym znajduje się skrypt programu należy zmienić scieżkę tego folderu na cmd), 
w otwrtym command window należy wpisac nazwę pliku proj2.py -m "nazwę modelu" -x 120 -y 200 -z 300 (nazwa modelu np. wgs84 lub grs80)
Wartości i argumenty oddzielamy spacją. 
Kolejnosć wyników: 1 - wsp. geo. 2 - Pl-2000 3 - PL-1992 4 - NEU (błedne wyniki)
- Za pomocą komendy wyglądajcej następująco  (python script_12.py -g dane.txt -t wyniki.txt) można w oknie cmd wprowadzic nazwe pliku (np.dane.txt) z danymi w formacie .txt oraz dodać nazwę pliku(np. wyniki.txt)->(plik utworzy sie w folderze w ktorym jest program) bedzie zawierac  obliczone wartości w formacie .txt. Przykładowy plik z danymi to wsp.inp.txt 

# Znane błędy i nietypowe zachowania programu, które nie zostały jeszcze naprawione:
- problemy z przeliczaniem współrzędnych xyz do NEU
