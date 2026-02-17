# SKRYPT INSTALACYJNY (SETUP)
# Uruchom ten skrypt tylko raz na nowym środowisku (np. po sklonowaniu repo).
# ------------------------------------------------------------------------

# 1. Instalacja menedżera środowisk 'renv' (jeśli brak)
if (!require('renv', quietly = TRUE)) install.packages('renv')

# 2. Aktywacja projektu
# (To sprawi, że biblioteki zainstalują się w folderze projektu, a nie globalnie)
renv::init(bare = TRUE) 

# 3. Instalacja menedżera pakietów bioinformatycznych
if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')

# 4. Definicja wymaganych pakietów

# A. Bioinformatyka (Bioconductor)
bio_packages <- c(
  'TCGAbiolinks',       # Pobieranie danych z TCGA
  'SummarizedExperiment', # Struktury danych genomowych
  'edgeR',              # Normalizacja (CPM)
  'ComplexHeatmap'      # Zaawansowane mapy ciepła
)

# B. Machine Learning & Statystyka
ml_packages <- c(
  'caret',              # Framework ML
  'randomForest',       # Algorytm Lasów Losowych
  'kernlab',            # Algorytm SVM
  'Boruta',             # Selekcja cech
  'e1071',              # Metryki statystyczne
  'pheatmap',           # Proste mapy ciepła
  'ggplot2',            # Wykresy statyczne
  'dplyr'               # Manipulacja danymi
)

# C. Aplikacja GUI (Shiny Dashboard) - NOWE!
gui_packages <- c(
  'shiny',              # Silnik aplikacji webowej
  'shinydashboard',     # Profesjonalny interfejs (dashboard)
  'plotly',             # Wykresy interaktywne
  'DT',                 # Interaktywne tabele
  'shinycssloaders'     # Animacje ładowania
)

# 5. Instalacja

message('--- KROK 1/3: Instalacja pakietów Bioinformatycznych ---')
BiocManager::install(bio_packages, update = FALSE, ask = FALSE)

message('--- KROK 2/3: Instalacja pakietów ML i GUI ---')
install.packages(c(ml_packages, gui_packages))

# 6. Zamrożenie stanu (Snapshot)
message('--- KROK 3/3: Aktualizacja pliku renv.lock ---')
# To kluczowe: zapisuje DOKŁADNE wersje, które masz teraz zainstalowane
renv::snapshot()

message('>>> Środowisko gotowe do pracy. <<<')