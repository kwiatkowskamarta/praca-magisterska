# SKRYPT 01: Data Download (TCGA-LAML)
# pobranie danych klinicznych i transkryptomicznych
# źródło: GDC (Genomic Data Commons) TCGAbiolinks

library(TCGAbiolinks)
library(SummarizedExperiment)

# utworzenie folderu na pliki tymczasowe
if(!dir.exists("GDCdata")) dir.create("GDCdata")

# dane kliniczne (etykiety)
message("Pobieranie danych klinicznych.")

query_clin <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

# pobieranie do folderu GDCdata
GDCdownload(query_clin, directory = "GDCdata")
clinical_data <- GDCprepare(query_clin, directory = "GDCdata")

# główna tabela pacjentów
clinical_patient <- clinical_data$clinical_patient_laml

message(paste("   Pobrano dane kliniczne dla", nrow(clinical_patient), "pacjentów."))

#dane RNA-seq (cechy)
message("Pobieranie danych RNA-seq.")

query_rna <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts" 
)

GDCdownload(query_rna, 
            method = "api", 
            files.per.chunk = 10,
            directory = "GDCdata")

data_rna <- GDCprepare(query_rna, directory = "GDCdata")

message(paste("   Pobrano macierz ekspresji:", nrow(data_rna), "genów x", ncol(data_rna), "próbek."))

# zapis suronych danych
message("Zapisywanie danych do folderu data/raw/.")

if(!dir.exists("data/raw")) dir.create("data/raw", recursive = TRUE)

save(data_rna, clinical_patient, file = "data/raw/TCGA_LAML_Raw.RData")

message(">>> Dane zostały zapisane w pliku 'data/raw/TCGA_LAML_Raw.RData'. <<<")