# SKRYPT INSTALACYJNY 
# instalacja menedżera środowisk 'renv' (jeśli brak)
if (!require('renv', quietly = TRUE)) install.packages('renv')

# aktywacja projektu
renv::init(bare = TRUE) 

# instalacja menedżera pakietów bioinformatycznych
if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')


# Bioconductor
bio_packages <- c(
  'TCGAbiolinks',       # Pobieranie danych z TCGA
  'SummarizedExperiment', # Struktury danych genomowych
  'edgeR',              # Normalizacja (CPM)
  'ComplexHeatmap'      # Zaawansowane mapy ciepła
)

# Machine Learning & Statystyka
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

# aplikacja GUI 
gui_packages <- c(
  'shiny',              # Silnik aplikacji webowej
  'shinydashboard',     # Profesjonalny interfejs (dashboard)
  'plotly',             # Wykresy interaktywne
  'DT',                 # Interaktywne tabele
  'shinycssloaders'     # Animacje ładowania
)

# instalacja
message('instalacja pakietów bioinformatycznych')
BiocManager::install(bio_packages, update = FALSE, ask = FALSE)

message('instalacja pakietów ML i GUI')
install.packages(c(ml_packages, gui_packages))

message('aktualizacja pliku renv.lock')
renv::snapshot()

message('Środowisko gotowe do pracy.')