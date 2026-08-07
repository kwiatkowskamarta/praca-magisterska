# SKRYPT 04: Visualizations
# Generowanie wykresów do analizy danych i modelu

library(ggplot2)
library(dplyr)
library(tidyr) 
library(Boruta)


load("data/processed/TCGA_LAML_Cleaned.RData")
load("results/models/Boruta_Results.RData")

message("Generowanie Wykresu 1: Rozkład klas FAB.")

# wykres słupkowy
plot_fab <- ggplot(ml_df, aes(x = Target, fill = Target)) +
  geom_bar(color = "black", alpha = 0.8) +
  # dodanie etykiet z dokładną liczbą nad słupkami
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5, size = 5) +
  theme_minimal() +
  labs(title = "Rozkład podtypów białaczki (Klasyfikacja FAB)",
       subtitle = "Zbiór danych TCGA-LAML po filtracji",
       x = "Podtyp FAB ",
       y = "Liczba pacjentów") +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", size = 14))

print(plot_fab)

ggsave("results/plot_1_fab_distribution.png", plot = plot_fab, width = 8, height = 6, dpi = 300)

message("Generowanie Wykresu 2: Zoptymalizowany wykres Boruty (Top 10).")

# 10 najważniejszych genów
imps <- attStats(boruta_output)
imps_confirmed <- imps[imps$decision == "Confirmed", ]
top10_genes <- rownames(imps_confirmed)[order(imps_confirmed$meanImp, decreasing = TRUE)[1:10]]

# historia ważności z modelu (tylko dla Top 10 i "cieni")
columns_to_keep <- c(top10_genes, "shadowMax", "shadowMean", "shadowMin")
imp_history <- as.data.frame(boruta_output$ImpHistory)
imp_subset <- imp_history[, colnames(imp_history) %in% columns_to_keep]

#transformacja tabeli
imp_df_long <- pivot_longer(imp_subset, cols = everything(), 
                            names_to = "Feature", values_to = "Importance")

# usunięcie braków danych/nieskonczonosci (częsty artefakt z pierwszej iteracji Boruty)
imp_df_long <- imp_df_long[is.finite(imp_df_long$Importance), ]

# oznaczenie genow i szumu
imp_df_long$Typ <- ifelse(grepl("shadow", imp_df_long$Feature), "Cienie (Szum losowy)", "Geny (Biomarkery)")

plot_boruta <- ggplot(imp_df_long, aes(x = reorder(Feature, Importance, FUN = median), 
                                       y = Importance, fill = Typ)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  coord_flip() + # odwrócenie osi dla lepszej czytelności nazw genów
  scale_fill_manual(values = c("Cienie (Szum losowy)" = "dodgerblue", "Geny (Biomarkery)" = "palegreen3")) +
  theme_minimal() +
  labs(title = "Ważność Top 10 genów vs Cechy losowe (Boruta)",
       subtitle = "Rozkład punktacji 'Importance' na przestrzeni iteracji algorytmu",
       x = "Gen / Cień",
       y = "Ważność (Z-score)") +
  theme(legend.title = element_blank(), 
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

print(plot_boruta)
ggsave("results/plot_2_boruta_top10.png", plot = plot_boruta, width = 8, height = 6, dpi = 300)

message("Generowanie Wykresu 3: Analiza PCA.")

# odizolowanie samych wartości liczbowych (macierz X) - bez kolumny Target
pca_data <- ml_df %>% select(-Target)

# wykonanie algorytmu PCA
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

# wyciągnięcie dwóch pierwszych osi (PC1 i PC2) dla każdego pacjenta
pca_df <- as.data.frame(pca_result$x[, 1:2])
# dodanie z powrotem diagnozy FAB
pca_df$Target <- ml_df$Target

# obliczenie ile informacji (% wariancji) zachowały nowe osie
variance <- summary(pca_result)$importance[2, 1:2] * 100

# rysowanie wykresu
plot_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Target)) +
  geom_point(size = 3.5, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Analiza Głównych Składowych (PCA) - Ekspresja RNA-seq",
       subtitle = "Projekcja danych transkryptomicznych na 2 wymiary",
       x = paste0("Główna Składowa 1 (", round(variance[1], 1), "% wariancji)"),
       y = paste0("Główna Składowa 2 (", round(variance[2], 1), "% wariancji)"),
       color = "Podtyp FAB") +
  theme(plot.title = element_text(face = "bold"))

print(plot_pca)
ggsave("results/plot_3_pca.png", plot = plot_pca, width = 8, height = 6, dpi = 300)


# SKRYPT 04: Visualizations
# Generowanie wykresów do analizy danych i modelu

library(ggplot2)
library(dplyr)
library(tidyr) 
library(Boruta)

load("data/processed/TCGA_LAML_Cleaned.RData")
load("results/models/Boruta_Results.RData")

message("Generowanie Wykresu 1: Rozkład klas FAB.")

# wykres słupkowy
plot_fab <- ggplot(ml_df, aes(x = Target, fill = Target)) +
  geom_bar(color = "black", alpha = 0.8) +
  # dodanie etykiet z dokładną liczbą nad słupkami
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5, size = 5) +
  theme_minimal() +
  labs(title = "Distribution of AML Subtypes (FAB Classification)",
       subtitle = "Filtered TCGA-LAML Dataset",
       x = "FAB Subtype",
       y = "Number of Patients") +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", size = 14))

print(plot_fab)

ggsave("results/plot_1_fab_distribution.png", plot = plot_fab, width = 8, height = 6, dpi = 300)

message("Generowanie Wykresu 2: Zoptymalizowany wykres Boruty (Top 10).")

# 10 najważniejszych genów
imps <- attStats(boruta_output)
imps_confirmed <- imps[imps$decision == "Confirmed", ]
top10_genes <- rownames(imps_confirmed)[order(imps_confirmed$meanImp, decreasing = TRUE)[1:10]]

# historia ważności z modelu (tylko dla Top 10 i "cieni")
columns_to_keep <- c(top10_genes, "shadowMax", "shadowMean", "shadowMin")
imp_history <- as.data.frame(boruta_output$ImpHistory)
imp_subset <- imp_history[, colnames(imp_history) %in% columns_to_keep]

# transformacja tabeli
imp_df_long <- pivot_longer(imp_subset, cols = everything(), 
                            names_to = "Feature", values_to = "Importance")

# usunięcie braków danych/nieskonczonosci (częsty artefakt z pierwszej iteracji Boruty)
imp_df_long <- imp_df_long[is.finite(imp_df_long$Importance), ]

# oznaczenie genow i szumu (ANGIELSKIE NAZWY)
imp_df_long$Typ <- ifelse(grepl("shadow", imp_df_long$Feature), "Shadow Features (Noise)", "Confirmed Genes (Biomarkers)")

plot_boruta <- ggplot(imp_df_long, aes(x = reorder(Feature, Importance, FUN = median), 
                                       y = Importance, fill = Typ)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  coord_flip() + # odwrócenie osi dla lepszej czytelności nazw genów
  scale_fill_manual(values = c("Shadow Features (Noise)" = "dodgerblue", "Confirmed Genes (Biomarkers)" = "palegreen3")) +
  theme_minimal() +
  labs(title = "Top 10 Feature Importance vs. Shadow Features",
       subtitle = "Z-score distribution across Boruta algorithm iterations",
       x = "Feature",
       y = "Importance (Z-score)") +
  theme(legend.title = element_blank(), 
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

print(plot_boruta)
ggsave("results/plot_2_boruta_top10.png", plot = plot_boruta, width = 8, height = 6, dpi = 300)

message("Generowanie Wykresu 3: Analiza PCA.")

# odizolowanie samych wartości liczbowych (macierz X) - bez kolumny Target
pca_data <- ml_df %>% select(-Target)

# wykonanie algorytmu PCA
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

# wyciągnięcie dwóch pierwszych osi (PC1 i PC2) dla każdego pacjenta
pca_df <- as.data.frame(pca_result$x[, 1:2])
# dodanie z powrotem diagnozy FAB
pca_df$Target <- ml_df$Target

# obliczenie ile informacji (% wariancji) zachowały nowe osie
variance <- summary(pca_result)$importance[2, 1:2] * 100

# rysowanie wykresu (ANGIELSKIE NAZWY)
plot_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Target)) +
  geom_point(size = 3.5, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Principal Component Analysis (PCA) of RNA-seq Data",
       subtitle = "2D Projection of Patient Transcriptomic Profiles",
       x = paste0("Principal Component 1 (", round(variance[1], 1), "% variance)"),
       y = paste0("Principal Component 2 (", round(variance[2], 1), "% variance)"),
       color = "FAB Subtype") +
  theme(plot.title = element_text(face = "bold"))

print(plot_pca)
ggsave("results/plot_3_pca.png", plot = plot_pca, width = 8, height = 6, dpi = 300)
library(ggplot2)

library(ggplot2)

# ==========================================================
# 1. WYKRES: Confusion Matrix (English version)
# ==========================

# Przygotowanie danych
cm_data <- data.frame(
  Reference = rep(c("M0", "M1", "M2", "M3", "M4", "M5"), times = 6),
  Prediction = rep(c("M0", "M1", "M2", "M3", "M4", "M5"), each = 6),
  Value = c(
    2, 0, 0, 0, 0, 0,
    0, 2, 2, 0, 1, 0,
    1, 4, 7, 0, 0, 0,
    0, 0, 0, 3, 0, 0,
    0, 2, 0, 0, 6, 1,
    0, 0, 0, 0, 0, 2
  )
)

cm_data$Prediction <- factor(cm_data$Prediction, levels = rev(c("M0", "M1", "M2", "M3", "M4", "M5")))
cm_data$Reference <- factor(cm_data$Reference, levels = c("M0", "M1", "M2", "M3", "M4", "M5"))

# Tworzenie wykresu
p1 <- ggplot(cm_data, aes(x = Reference, y = Prediction, fill = Value)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = Value), fontface = "bold", size = 5) +
  scale_fill_gradient(low = "white", high = "#4388C5", name = "Number of\nsamples") +
  labs(
    title = "Confusion Matrix - Random Forest",
    subtitle = "Comparison of model predictions with actual FAB subtypes on the test set",
    x = "Actual physician diagnosis (Reference)",
    y = "Predicted diagnosis by the model (Prediction)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle = element_text(size = 12, hjust = 0, margin = margin(b = 15)),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 11, color = "gray40"),
    panel.grid = element_blank()
  )

# Wyświetlenie wykresu
print(p1)

# AUTOMATYCZNY ZAPIS DO PLIKU
ggsave("confusion_matrix_en.png", plot = p1, width = 10, height = 8, dpi = 300)


# ==========================================================
# 2. WYKRES: FAB Distribution (English version)
# ==========================================================

# Przygotowanie danych
bar_data <- data.frame(
  Subtype = factor(c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7"), 
                   levels = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7")),
  Count = c(15, 35, 38, 14, 29, 15, 3, 2)
)

# Tworzenie wykresu
p2 <- ggplot(bar_data, aes(x = Subtype, y = Count)) +
  geom_bar(stat = "identity", fill = "#4B8CDB", width = 0.65) +
  geom_text(aes(label = Count), vjust = -0.6, fontface = "bold", size = 4.5) +
  labs(
    title = "Distribution of FAB Subtype Classes",
    x = "FAB Subtype",
    y = "Number of patients"
  ) +
  scale_y_continuous(limits = c(0, 42), breaks = seq(0, 30, by = 10)) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 20)),
    axis.title.x = element_text(face = "bold", margin = margin(t = 15)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 15)),
    axis.text = element_text(size = 11, color = "gray30"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5)
  )

# Wyświetlenie wykresu
print(p2)

# AUTOMATYCZNY ZAPIS DO PLIKU
ggsave("fab_distribution_en.png", plot = p2, width = 12, height = 7, dpi = 300)
