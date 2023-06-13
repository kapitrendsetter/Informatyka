# Projekt Transformacje
Skrypt potrzebny do uzywania biblioteki argparse znajduje sie w pliku proj2.py i script_12.py

# Q&A

# Do czego służy program i jaką funkcjonalność oferuje?  
Transformacje:
- XYZ -> BLH, 
- BLH -> XYZ 
- XYZ -> NEUp
- BL(GRS80, WGS84,) -> 2000 
- BL(GRS80, WGS84,) -> 1992

# Jakie elipsidy są obsługiwane?
- GRS80 
- WGS84

# Jakie wymagania trzeba spełnić, by program działał na danym komputerze?
- trzeba mieć pythona w wersji 10.9 lub nowszej 
- zainstalowaną bibliotekę argparse (która miewa problemy z odczytem)

# Dla jakiego systemu operacyjnego został napisany program?
- program działa zarówno na Windows 10 jak i macOS Ventura 13.0

# Instrukcja jak używać program na 2 sposoby: 

- Wykorzystanie biblioteki argparse w celu przeliczenia współrzędnych w wybranym odwzorowaniu: \
Należy otworzyć command window ( w folderze w którym znajduje się skrypt programu należy zmienić scieżkę tego folderu na cmd), 
w otwrtym command window należy wpisac nazwę pliku proj2.py -m "nazwę modelu" -x 120 -y 200 -z 300 (nazwa modelu np. wgs84 lub grs80)
Wartości i argumenty oddzielamy spacją. \
Kolejnosć wyników: 1 - wsp. geo. 2 - Pl-2000 3 - PL-1992 4 - NEU (błedne wyniki)

- Wykorzystanie biblioteki argparse w celu wprowadzania danych za pomoca pliku tekstowego w wybranym odwzorowaniu: \
przykładowy kod: \
*python script_12.py --dane XYZ --model grs80 --uklad BLH -txt dane.txt -txt_out wyniki.txt* \
kolejno: \
--dane XYZ (w metrach) lub BLH (BL w stopniach dziesietnych i H w metrach) \
--model grs80 lub wgs84 \
Wyjatkowo dla NEU przykladowy kod wyglada: \
*python script_12.py --dane XYZ --model grs80 --uklad NEU -x_ref 3665000 -y_ref 1410000 -z_ref 5010000 -txt dane.txt -txt_out wyniki2.txt*

dodatkowo nalezy podac:
nazwe pliku (np.dane.txt) z danymi w formacie .txt oraz dodać nazwę pliku(np. wyniki.txt)->(plik utworzy sie w folderze w ktorym jest program) bedzie zawierac  obliczone wartości w formacie .txt

**Przykładowy plik z wprowadzonymi danymi to wsp.inp.txt** 

# Znane błędy i nietypowe zachowania programu, które nie zostały jeszcze naprawione:
- problemy z przeliczaniem współrzędnych xyz do NEU

