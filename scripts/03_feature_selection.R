# SKRYPT 03: Feature Selection (Boruta)
# identyfikacja biomarkerów genetycznych bez wstępnej redukcji
# algorytm Boruta

library(Boruta)
library(dplyr)
library(caret)

set.seed(42)

# wczytanie danych
load("data/processed/TCGA_LAML_Cleaned.RData")

message("Wczytano dane.")
message("   Liczba pacjentów: ", nrow(ml_df))
message("   Liczba genów do sprawdzenia: ", ncol(ml_df) - 1)

# przygotowanie danych (x i y)
# rozdzielenie macierzy na features i target
target <- ml_df$Target
features <- ml_df %>% select(-Target)

# boruta
message("\nRozpoczynam obliczenia algorytmem Boruta.")

boruta_output <- Boruta(
  x = features,
  y = target,
  doTrace = 2, #szczegolowy podglad postepu
  maxRuns = 200 
)

message("\nKoniec obliczen.")

# wyniki
print(boruta_output)

# naprawa "tentative" -- jeśli po 200 rundach Boruta nadal nie jest pewna co do niektórych genów, zmuszam ją do podjęcia decyzji (TentativeRoughFix).
if(any(boruta_output$finalDecision == "Tentative")) {
  message("Znaleziono atrybuty niepewne (Tentative) - ostateczna decyzja...")
  boruta_output <- TentativeRoughFix(boruta_output)
  print(boruta_output)
}

# ekstrakcja kluczowych genów - tylko "confirmed"
final_genes <- getSelectedAttributes(boruta_output, withTentative = FALSE)
message("\nLiczba zidentyfikowanych biomarkerów: ", length(final_genes))

# top 10 najważniejszych (wg "importance")
imps <- attStats(boruta_output)
imps_confirmed <- imps[imps$decision == "Confirmed", ]
top_genes <- row.names(imps_confirmed)[order(imps_confirmed$meanImp, decreasing = TRUE)]

message("Top 10 najważniejszych genów:")
print(head(top_genes, 10))

# zapis
# (pełny obiekt Boruta (do wykresów) oraz dataset do modeli)
ml_data_final <- ml_df %>% 
  select(all_of(final_genes), Target)

if(!dir.exists("results/models")) dir.create("results/models", recursive = TRUE)

save(boruta_output, ml_data_final, final_genes, file = "results/models/Boruta_Results.RData")

message(">>> Wyniki zapisane w results/models/Boruta_Results.RData <<<")