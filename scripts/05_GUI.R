library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(dplyr)
library(caret)
library(randomForest)
library(ggplot2)

load_data_smart <- function() {
  path_root <- "results/models/"
  path_shiny <- "../results/models/"
  
  if (file.exists(paste0(path_root, "Boruta_Results.RData"))) {
    base_path <- path_root
  } else if (file.exists(paste0(path_shiny, "Boruta_Results.RData"))) {
    base_path <- path_shiny
  } else {
    return(NULL)
  }
  
  message(paste("Wczytuję dane z folderu:", base_path))
  
  load(paste0(base_path, "Boruta_Results.RData"), envir = .GlobalEnv)
  load(paste0(base_path, "Final_Models.RData"), envir = .GlobalEnv)
  return(base_path)
}

base_path <- load_data_smart()
if(is.null(base_path)) {
  stop("BŁĄD: Nie znaleziono plików .RData. Uruchom skrypty 03 i 04.")
}

if (!exists("model_rf") && exists("models_list")) {
  model_rf <- models_list$Random_Forest
}

# przygotowanie danych do wykresów
df <- ml_data_final
all_genes <- sort(as.character(final_genes)) # Pełna lista genów
top_genes <- as.character(final_genes[1:20])

# UI
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Klasyfikacja AML"),
  
  # pasek boczny
  dashboardSidebar(
    sidebarMenu(
      menuItem("strona główna", tabName = "dashboard"),
      menuItem("przegląd danych", tabName = "eda"),
      menuItem("wyniki modeli", tabName = "models"),
      menuItem("analiza biomarkerów", tabName = "genes"),
      menuItem("symulator diagnostyczny", tabName = "simulator")
    )
  ),
  
  dashboardBody(
    tabItems(
      
      # strona glowna
      tabItem(tabName = "dashboard",
              h2("System wspomagania diagnostyki Ostrej Białaczki Szpikowej"),
              br(),
              fluidRow(
                valueBoxOutput("box_patients")
              ),
              fluidRow(
                box(title = "Opis Projektu", status = "primary", solidHeader = TRUE, width = 12,
                    "Aplikacja stanowi część pracy magisterskiej. Wykorzystano dane transkryptomiczne (RNA-seq) z projektu TCGA-LAML https://portal.gdc.cancer.gov/projects/TCGA-LAML. " ,
                    "Kluczowe etapy analizy obejmują: czyszczenie danych, selekcję cech oraz trening modeli służących do klasyfikacji podtypów Ostrej Białaczki Szpikowej. ")
              )
      ),
      
      # pca i rozklad podtypoe
      tabItem(tabName = "eda",
              h2("Eksploracyjna Analiza Danych"),
              fluidRow(
                box(title = "Rozkład Podtypów Białaczki", status = "primary", solidHeader = TRUE,
                    plotlyOutput("plot_class_dist")),
                box(title = "Wizualizacja PCA (2D)", status = "primary", solidHeader = TRUE,
                    plotlyOutput("plot_pca"),
                    helpText("Każdy punkt to pacjent. Bliskość punktów oznacza podobieństwo genetyczne."))
              )
      ),
      
      # modele
      tabItem(tabName = "models",
              h2("Ewaluacja Modeli ML"),
              fluidRow(
                box(title = "Porównanie Dokładności (Accuracy)", width = 6,
                    plotlyOutput("plot_model_compare")),
                box(title = "Macierz Pomyłek (Random Forest)", width = 6,
                    plotlyOutput("plot_confusion_matrix"),
                    helpText("Wiersze: Prawdziwa diagnoza. Kolumny: Predykcja modelu."))
              )
      ),
      
      # biomarkery
      tabItem(tabName = "genes",
              h2("Analiza Biomarkerów Genetycznych"),
              fluidRow(
                column(width = 4,
                       box(title = "Wybierz Gen", width = 12, status = "warning",
                           # Zostawiamy puste choices - wypełnimy je z serwera (server side)
                           selectInput("selected_gene", "Gen do analizy:", choices = NULL),
                           helpText("Wybrano 20 najważniejszych genów wg algorytmu Boruta.")
                       )
                ),
                column(width = 8,
                       box(title = "Ekspresja Genu w Podtypach AML", width = 12, solidHeader = TRUE,
                           plotlyOutput("plot_gene_boxplot"))
                )
              )
      ),
      
      # symulator
      tabItem(tabName = "simulator",
              h2("Symulator Diagnostyczny"),
              fluidRow(
                box(title = "Panel Pacjenta", status = "danger", solidHeader = TRUE, width = 4,
                    helpText("Wybierz pacjenta ze zbioru testowego, aby sprawdzić działanie modelu."),
                    selectInput("sim_patient_id", "ID Pacjenta:", choices = rownames(test_set)),
                    actionButton("btn_predict", "ZDIAGNOZUJ", 
                                 class = "btn-danger btn-lg", style = "width:100%; margin-top: 20px;")
                ),
                box(title = "Wynik Diagnostyki AI", status = "success", solidHeader = TRUE, width = 8,
                    h3(textOutput("pred_class")),
                    hr(),
                    plotlyOutput("plot_pred_prob")
                )
              )
      )
    )
  )
)

server <- function(input, output, session) {
  
  observe({
    updateSelectInput(session, "selected_gene",
                      choices = all_genes,
                      selected = all_genes[1])
  })
  
  # 1 logika (str glowna)
  output$box_patients <- renderValueBox({
    valueBox(nrow(df), "Liczba Pacjentów", color = "light-blue")
  })
  
  # 2 logika (pca i rozklad)
  output$plot_class_dist <- renderPlotly({
    p <- ggplot(df, aes(x = Target, fill = Target)) +
      geom_bar() +
      theme_minimal() +
      labs(x = "Podtyp", y = "Liczba pacjentów") +
      theme(legend.position = "none")
    ggplotly(p)
  })
  
  output$plot_pca <- renderPlotly({
    pca_res <- prcomp(df[, -ncol(df)], scale. = TRUE)
    pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Target = df$Target, ID = rownames(df))
    
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Target, text = ID)) +
      geom_point(size = 3, alpha = 0.8) +
      theme_minimal() +
      labs(title = "Analiza Głównych Składowych")
    ggplotly(p)
  })
  
  # 3 logika (modele)
  output$plot_model_compare <- renderPlotly({
    # Wyciągamy dane z obiektu results (bezpiecznie - sprawdzenie czy kNN istnieje)
    knn_acc <- if("kNN" %in% names(models_list)) mean(results$values$`kNN~Accuracy`) else 0
    
    res_df <- data.frame(
      Model = c("Random Forest", "k-NN"),
      Accuracy = c(mean(results$values$`Random_Forest~Accuracy`), knn_acc)
    )
    
    p <- ggplot(res_df, aes(x = Model, y = Accuracy, fill = Model)) +
      geom_col(width = 0.5) +
      ylim(0, 1) +
      theme_minimal() +
      labs(y = "Średnia Dokładność (CV)")
    ggplotly(p)
  })
  
  output$plot_confusion_matrix <- renderPlotly({
    cm_table <- as.data.frame(cm_rf$table)
    p <- ggplot(cm_table, aes(x = Reference, y = Prediction, fill = Freq)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "steelblue") +
      geom_text(aes(label = Freq)) +
      theme_minimal() +
      labs(x = "Prawdziwa Klasa", y = "Predykcja Modelu")
    ggplotly(p)
  })
  
  # 4 logika (geny)
  output$plot_gene_boxplot <- renderPlotly({
    req(input$selected_gene)
    
    # Bezpieczne pobieranie danych (zabezpieczenie przed kropkami w nazwach genów)
    plot_data <- df[, c("Target", input$selected_gene)]
    colnames(plot_data) <- c("Target", "Expression") # Ujednolicamy nazwę
    
    p <- ggplot(plot_data, aes(x = Target, y = Expression, fill = Target)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.3) +
      theme_minimal() +
      labs(title = paste("Gen:", input$selected_gene), y = "Poziom Log2 CPM") +
      theme(legend.position = "none")
    ggplotly(p)
  })
  
  # 5 logika (symualtor)
  observeEvent(input$btn_predict, {
    req(input$sim_patient_id)
    
    # dane pacjenta
    patient_data <- df[input$sim_patient_id, , drop = FALSE]
    
    #predykcja klasy
    pred_class <- predict(model_rf, newdata = patient_data)
    
    # predykcja prawdopodobieństw
    pred_prob <- predict(model_rf, newdata = patient_data, type = "prob")
    prob_df <- data.frame(Klasa = colnames(pred_prob), Prawdopodobienstwo = as.numeric(pred_prob[1,]))
    
    output$pred_class <- renderText({
      paste("Zdiagnozowano podtyp:", as.character(pred_class))
    })
    
    # wykres słupkowy prawdopodobieństw
    output$plot_pred_prob <- renderPlotly({
      p <- ggplot(prob_df, aes(x = Klasa, y = Prawdopodobienstwo, fill = Klasa)) +
        geom_col() +
        ylim(0, 1) +
        geom_text(aes(label = round(Prawdopodobienstwo, 2)), vjust = -0.5) +
        theme_minimal() +
        labs(title = "Pewność modelu dla poszczególnych klas", y = "Prawdopodobieństwo")
      ggplotly(p)
    })
  })
}


shinyApp(ui = ui, server = server)