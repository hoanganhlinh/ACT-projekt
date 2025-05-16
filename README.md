# METABRIC Breast Cancer dataset Survival Analysis Project
Repository dedicated to an academic project on Survival Analysis (Subject: *Analiza Czasu Trwania*).

1. **Dataset name and source**
   **Breast Cancer (METABRIC) Clinical Data**
   Source: [https://www.kaggle.com/datasets/gunesevitan/breast-cancer-metabric](https://www.kaggle.com/datasets/gunesevitan/breast-cancer-metabric)

2. **Number of observations and variables**
   2,509 observations (after filtering out cases where *Cancer Type* = "Breast Sarcoma" and removing rows with missing data, 854 observations remain), 34 variables

3. **Start event**
   Inclusion of the patient in the observation/study.

4. **Event of interest (eliminating event)**
   Death due to breast cancer (*"Died of Disease"* in the column *Patient's Vital Status*).

5. **Primary variable**
   Cancer type – only breast cancer cases are included (*Cancer Type* = "Breast Cancer"); other types such as "Breast Sarcoma" are excluded.

6. **Secondary variables (examples)**

* Age at Diagnosis
* Type of Breast Surgery
* Cancer Type Detailed
* Cellularity
* Chemotherapy
* HER2 Status
* Hormone Therapy
* Inferred Menopausal State
* Integrative Cluster
* Primary Tumor Laterality
* Lymph Nodes Examined Positive
* Mutation Count
* Radio Therapy
* 3-Gene Classifier Subtype
* Tumor Size
* Tumor Stage

7. **Definition of the time-to-event variable**
   The variable *Overall Survival (Months)* – number of months from diagnosis to death or last follow-up.

8. **Definition of the censoring variable**
   The censoring indicator takes the value 1 for uncensored cases (i.e., patient died): when *Patient's Vital Status* = "Died of Disease" (370 cases).
   It takes the value 0 for censored cases (i.e., patient did not die): when *Patient's Vital Status* = "Living" (484 cases).
