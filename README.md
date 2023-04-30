# Projekt Transformacje
Skryp właściwy znajduje się w pliku script_12.py, a w pliku script_12-14.py zapisane sa wszystkie commity  (zostalismy zmuszeni do tego 
ponieważ python czyta - jako minus i nie mozna zaimportowac funkcji do 2 skryptu.
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

- Wykorzystanie biblioteki argparse
Najpierw należy otworzyć command window ( w folderze w którym znajduje się skrypt programu należy zmienić scieżkę tego folderu na cnd), 
w otwrtym command window należy wpisac nazwę pliku proj2.py -m"nazwę modelu" -x 120 -y200 -z300 (nazwa modelu np. WGS84 lub GRS80)
Wartości i argumenty oddzielamy spacją 
Kolejnosć wyników: 1 - wsp. geo. 2 - Pl-2000 3 - PL-1992 4 - NEU


    
# Znane błędy i nietypowe zachowania programu, które nie zostały jeszcze naprawione:
- problemy z wywoływaniem w bibliotece argparse 
- problemy z przeliczaniem współrzędnych xyz do NEU
