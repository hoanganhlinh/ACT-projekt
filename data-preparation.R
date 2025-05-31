library(readr)
library(dplyr)    
library(stringr)  
library(tidyr)
library(ggplot2)

# 2. Wczytanie danych
df_metabric <- read_csv("./data/Breast Cancer METABRIC.csv")
df_metabric

# 3. Unikalne wartości w kolumnie "Cancer Type"
unique(df_metabric$`Cancer Type`)

# 4. Nazwy kolumn
colnames(df_metabric)

######
#Przygotowanie danych 
# 1. Zmiana nazw kolumn

przeksztalc_nazwe <- function(nazwa) {
  nazwa <- str_replace_all(nazwa, "[()]", "")
  nazwa <- str_replace_all(nazwa, "[^a-zA-Z]", "_")
  nazwa <- str_replace_all(nazwa, "_+", "_")
  nazwa <- str_replace(nazwa, "^_+|_+$", "")
  nazwa <- tolower(nazwa)
  return(nazwa)
}

# Przekształć wszystkie nazwy kolumn w ramce
colnames(df_metabric) <- sapply(colnames(df_metabric), przeksztalc_nazwe)

colnames(df_metabric)

# 2. Usunięcie pozostałych typów nowotworów 
df_metabric <- df_metabric %>%
  filter(cancer_type == "Breast Cancer")

# Usuń kolumnę cancer_type
df_metabric <- df_metabric %>%
  select(-cancer_type)

df_metabric

# 3. Usuwanie wartości zbędnych wartości Patient_s_Vital_Status + definicja event (c)

unique(df_metabric$patient_s_vital_status)

df_metabric <- df_metabric %>%
  filter(patient_s_vital_status != "Died of Other Causes")


df_metabric <- df_metabric %>%
  mutate(
    c = case_when(
      patient_s_vital_status == "Died of Disease" ~ 1L,
      patient_s_vital_status == "Living"          ~ 0L,
      TRUE                                         ~ NA_integer_
    )
  )


df_metabric <- df_metabric %>%
  select(-patient_s_vital_status)


df_metabric

# 4. Definicja zmiennej t
df_metabric <- df_metabric %>%
  rename(
    t = overall_survival_months
  )

# 5. Sprawdzanie braków danych
df_metabric <- df_metabric %>%
  filter(
    !is.na(t),
    !is.na(c)
  )
braki <- colSums(is.na(df_metabric))
procent_brakow <- (braki / nrow(df_metabric)) * 100

braki_info <- data.frame(
  Liczba_brakow = braki,
  Procent_brakow = procent_brakow
)

braki_info$Kolumna <- rownames(braki_info)
rownames(braki_info) <- NULL

braki_info_filtered <- braki_info[braki_info$Liczba_brakow > 0, ]
print(braki_info_filtered)

df_metabric <- na.omit(df_metabric)

df_before <- df_metabric
# 6. Kategoryzacja zmiennych ciągłych
df_metabric <- df_metabric %>%
  mutate(
    age_c = case_when(
      age_at_diagnosis <  49                                 ~ 1L,
      age_at_diagnosis >= 49 & age_at_diagnosis <  59        ~ 2L,
      age_at_diagnosis >= 59 & age_at_diagnosis <  68        ~ 3L,
      age_at_diagnosis >= 68                                 ~ 4L,
      TRUE                                                   ~ NA_integer_
    ),
    ts_c = case_when(
      tumor_size == 0                                        ~ 0L,
      tumor_size > 0  & tumor_size <  5                      ~ 1L,
      tumor_size >= 5 & tumor_size < 10                      ~ 2L,
      tumor_size >= 10 & tumor_size < 20                     ~ 3L,
      tumor_size >= 20 & tumor_size < 50                     ~ 4L,
      tumor_size >= 50 & tumor_size < 80                     ~ 5L,
      tumor_size >= 80                                        ~ 6L,
      TRUE                                                   ~ NA_integer_
    ),
    npi_c = case_when(
      nottingham_prognostic_index <  1                       ~ 1L,
      nottingham_prognostic_index >= 1 & nottingham_prognostic_index <  2 ~ 2L,
      nottingham_prognostic_index >= 2 & nottingham_prognostic_index <  3 ~ 3L,
      nottingham_prognostic_index >= 3 & nottingham_prognostic_index <  4 ~ 4L,
      nottingham_prognostic_index >= 4 & nottingham_prognostic_index <  5 ~ 5L,
      nottingham_prognostic_index >= 5 & nottingham_prognostic_index <  6 ~ 6L,
      nottingham_prognostic_index >= 6                        ~ 7L,
      TRUE                                                   ~ NA_integer_
    ),
    lnep_c = case_when(
      lymph_nodes_examined_positive == 0                                ~ 0L,
      lymph_nodes_examined_positive > 0  & lymph_nodes_examined_positive <= 2  ~ 1L,
      lymph_nodes_examined_positive > 2  & lymph_nodes_examined_positive <= 6  ~ 2L,
      lymph_nodes_examined_positive > 6  & lymph_nodes_examined_positive <= 15 ~ 3L,
      lymph_nodes_examined_positive > 15                                    ~ 4L,
      TRUE                                                                 ~ NA_integer_
    )
  ) %>%
  select(
    -age_at_diagnosis,
    -tumor_size,
    -nottingham_prognostic_index,
    -lymph_nodes_examined_positive
  )

glimpse(df_metabric)

# 7. Usunięcie identyfikatora oraz kolumn powiązanych
df_metabric <- df_metabric %>%
  select(
    -patient_id,
    -er_status_measured_by_ihc,
    -her_status_measured_by_snp,
    -overall_survival_status,
    -relapse_free_status_months,
    -relapse_free_status,
    -sex
  )
View(df_metabric)

# 8. Dummy classifier na wszystkich kategorycznych kolumnach
df_metabric_num <- df_metabric

df_metabric_num <- df_metabric_num %>%
  mutate(
    across(
      .cols = where(is.character),
      .fns  = ~ as.integer(factor(., levels = unique(.))) - 1
    )
  )

glimpse(df_metabric_num)

write.csv(df_metabric_num, "metabric_cleaned.csv", row.names = FALSE)
View(df_metabric_num)

# TABLICE KONTYNGENCJI
count_censoring_separately <- function(df, varlist, censor_col = "c") {
  tables <- list()
  for (var in varlist) {
    ct <- table(df[[var]], df[[censor_col]])
    if (!("0" %in% colnames(ct))) {
      ct <- cbind(ct, "0" = 0)
    }
    if (!("1" %in% colnames(ct))) {
      ct <- cbind(ct, "1" = 0)
    }
    colnames(ct) <- c("not_censored (c=0)", "censored (c=1)")
    tables[[var]] <- ct
  }
  tables
}

varlist <- setdiff(names(df_metabric), c("t", "c"))
tables_dict <- count_censoring_separately(df_metabric, varlist)

for (col in varlist) {
  cat("\n", col, "\n")
  print(tables_dict[[col]])
}

### Histogramy 

primary_vars <- setdiff(
  colnames(df_metabric),
  c("c", "age_c", "ts_c", "npi_c", "lnep_c")
)

for (var in primary_vars) {
  if (is.numeric(df_metabric[[var]])) {
    rng <- range(df_metabric[[var]], na.rm = TRUE)
    bw <- (rng[2] - rng[1]) / 30
    print(
      ggplot(df_metabric, aes_string(x = var)) +
        geom_histogram(binwidth = bw) +
        ggtitle(var)
    )
  } else {
    print(
      ggplot(df_metabric, aes_string(x = var)) +
        geom_bar() +
        ggtitle(var)
    )
  }
}

pairs <- list(
  age_at_diagnosis = "age_c",
  tumor_size = "ts_c",
  nottingham_prognostic_index = "npi_c",
  lymph_nodes_examined_positive = "lnep_c"
)

for (cont in names(pairs)) {
  cat_var <- pairs[[cont]]
  rng_before <- range(df_before[[cont]], na.rm = TRUE)
  bw_before <- (rng_before[2] - rng_before[1]) / 30
  print(
    ggplot(df_before, aes_string(x = cont)) +
      geom_histogram(binwidth = bw_before) +
      ggtitle(paste0(cont, " przed kategoryzacją"))
  )
  print(
    ggplot(df_metabric, aes_string(x = cat_var)) +
      geom_bar() +
      ggtitle(paste0(cat_var, " po kategoryzacji"))
  )
}
