# SKRYPT 02: Preprocessing (Czyszczenie i Normalizacja)
# przygotowanie tabeligotowych do analizy z surowych danych.
# input: data/raw/TCGA_LAML_Raw.RData
# output: data/processed/TCGA_LAML_Cleaned.RData

library(SummarizedExperiment)
library(dplyr)
library(edgeR) 

# wczytanie surowych danych
load("data/raw/TCGA_LAML_Raw.RData")

# macierz (geny x próbki)
counts_matrix <- assay(data_rna)
message("Wczytano macierz: ", nrow(counts_matrix), " genów x ", ncol(counts_matrix), " próbek.")

# przygotowanie danych kllinicznych
target_col <- "fab_category"

# sprawdzenie czy kolumna istnieje
if(!target_col %in% colnames(clinical_patient)) {
  stop(paste("BŁĄD: Nie znaleziono kolumny", target_col))
}

message("Czyszczenie danych klinicznych.")

# ID pacjenta i target, czyszczenie śmieci
clin_clean <- clinical_patient %>%
  select(bcr_patient_barcode, target = all_of(target_col)) %>%
  # usunięcie metadanych (wiersze nagłówkowe z pliku Biotab)
  filter(!grepl("CDE_ID", target)) %>%            # usuniecie wiersza z ID
  filter(!grepl("morphology_code", target)) %>%   # usuniecie wiersza z długim napisem
  # usunięcie braków diagnozy
  filter(target != "Not Classified") %>%
  filter(!is.na(target)) %>%
  # ujednolicenie nazw (M0 Undifferentiated -> M0)
  mutate(target = ifelse(target == "M0 Undifferentiated", "M0", target))

message("Liczba pacjentów po wyczyszczeniu diagnozy: ", nrow(clin_clean))
print(table(clin_clean$target)) # podgląd

# parowanie próbek (geny <-> pacjenci) 
# ID w RNA (np. TCGA-AB-2805-03A...) są dłuższe niż w klinice (TCGA-AB-2805).
rna_barcodes <- colnames(counts_matrix)
rna_patient_ids <- substr(rna_barcodes, 1, 12) 

# wspólne ID
common_patients <- intersect(clin_clean$bcr_patient_barcode, rna_patient_ids)
common_patients <- sort(common_patients) #sortowanie zeby kolejnosc sie zgadzala wszedzie
message("Liczba pacjentów do sparowania: ", length(common_patients))

if(length(common_patients) < 50) {
  stop("UWAGA: Zbyt mało wspólnych pacjentów! Sprawdź formaty ID.")
}

# filtracja obu tabel
# a) klinicznie - tylko ci, co mają RNA
clin_final <- clin_clean %>% 
  filter(bcr_patient_barcode %in% common_patients) %>%
  arrange(match(bcr_patient_barcode, common_patients)) #alfabetycznie

# b) genetycznie - tylko ci, co mają diagnozę
match_indices <- match(common_patients, rna_patient_ids)
counts_final <- counts_matrix[, match_indices]

colnames(counts_final) <- common_patients

# ostateczne sprawdzenie kolejności
if(!all(colnames(counts_final) == clin_final$bcr_patient_barcode)) {
  stop("BŁĄD KRYTYCZNY: Kolejność pacjentów się nie zgadza!")
}

# filtracja genów (szum)
# usuniecie genow, które mają bardzo niską ekspresję
# (zostaja geny, które mają > 5 zliczeń u przynajmniej 20% pacjentów)
keep_genes <- rowSums(counts_final > 5) >= (0.2 * ncol(counts_final))
counts_filtered <- counts_final[keep_genes, ]

message("Filtracja genów: zredukowano z ", nrow(counts_final), " do ", nrow(counts_filtered))

# normalizacja (Log2 CPM)
message("Normalizacja danych (Log2 CPM).")
dge <- DGEList(counts = counts_filtered)
dge <- calcNormFactors(dge) # TMM normalization
cpm_log <- cpm(dge, log = TRUE, prior.count = 1)

# transpozycja  (do dalszej analizy)
ml_matrix <- t(cpm_log) 
ml_df <- as.data.frame(ml_matrix)
ml_df$Target <- factor(clin_final$target)

# usuniecie rzadki klas (M6 i M7)
print(table(ml_df$Target))
ml_df <- ml_df %>%
  filter(!Target %in% c("M6", "M7"))
ml_df$Target <- droplevels(ml_df$Target)

# czyszczenie nazwy kolumn (geny)
colnames(ml_df) <- make.names(colnames(ml_df))

message("Gotowa tabela do analizy: ", nrow(ml_df), " wierszy x ", ncol(ml_df), " kolumn.")

# zapis
if(!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
save(ml_df, file = "data/processed/TCGA_LAML_Cleaned.RData")

message(">>> Dane gotowe do analizy. <<<")