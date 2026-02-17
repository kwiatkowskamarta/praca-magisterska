# praca magisterska
## Analiza porównawcza skuteczności algorytmów uczenia maszynowego do klasyfikacji podtypów białaczki.

### Krok po kroku

1.  **Sklonuj repozytorium:**
    ```bash
    git clone https://github.com/kwiatkowskamarta/AML_Classification.git
    ```
    
2.  **Otwórz projekt:** Kliknij plik `AML_Classification.Rproj

3.  **Zainstaluj biblioteki:**
    Uruchom komendę, która odtworzy środowisko z wymaganymi wersjami pakietów:
    ```r
    renv::restore()
    ```
    
4.  **Uruchom analizę:**
    ```r
    source("scripts/01_data_download.R")
    source("scripts/02_preprocessing.R")
    source("scripts/03_feature_selection.R")
    ```
