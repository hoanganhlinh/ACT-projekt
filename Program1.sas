options encoding='utf-8';
/* KROK 1: Filtrowanie – zostawiamy tylko przypadki raka piersi         */
/************************************************************************/
data work.metabric_m;
  set work."BREAST CANCER METABRIC"n;
  if "Cancer Type"n = 'Breast Cancer';
run;
/* KROK 2: Usunięcie kolumny 'Cancer Type', bo już nie będzie potrzebna */
/************************************************************************/
data work.metabric_m;
  set work.metabric_m(drop="Cancer Type"n);
run;
/* KROK 3: Sprawdzenie dostępnych kategorii w kolumnie Vital Status     */
/************************************************************************/
proc freq data=work.metabric_m;
  tables "Patient's Vital Status"n / missing;
  title "Kategorie w Patient's Vital Status";
run;

/* KROK 4: Zostawiamy tylko Living oraz Died of Disease                 */
/************************************************************************/
data work.metabric_m;
  set work.metabric_m;
  if "Patient's Vital Status"n in ('Living', 'Died of Disease');
run;

/* KROK 5: Tworzymy zmienną event                                       */
/************************************************************************/
data work.metabric_m;
  set work.metabric_m;
  if "Patient's Vital Status"n = 'Died of Disease' then event = 1;
  else if "Patient's Vital Status"n = 'Living'      then event = 0;
run;

proc freq data=work.metabric_m;
  tables 
    Cohort
    _character_
    _numeric_
    / missing noprint out=freqs(drop=percent cumfreq cumpercent);
run;

/* Podstawowe statystyki i liczba brak�w dla wszystkich numerycznych */
/************************************************************************/
proc means data=work.metabric_m 
           n nmiss mean std min p25 median p75 max;
  var _numeric_;
  title 'Statystyki opisowe i liczba brak�w dla zmiennych numerycznych';
run;

proc sql;
  select name, type
    from dictionary.columns
   where libname='WORK'
     and memname='METABRIC_M'
   order by type, name;
quit;



proc freq data=work.metabric_m;
  tables
 
    Cellularity
    Chemotherapy
	Sex
  
    "3-Gene classifier subtype"n
    "Cancer Type Detailed"n
    "ER status measured by IHC"n
    "Pam50 + Claudin-low subtype"n
    "Patient's Vital Status"n
	"ER Status"n
	"HER2 Status"n
	"HER2 status measured by SNP6"n
	"Hormone Therapy"n
	"Inferred Menopausal State"n
	"Integrative Cluster"n
	"Oncotree Code"n
	"Overall Survival Status"n
	"PR Status"n
	"Primary Tumor Laterality"n
	"Radio Therapy"n
	"Relapse Free Status"n
	"Tumor Other Histologic Subtype"n
	"Type of Breast Surgery"n

   
  / missing;  
  title 'Rozk�ady zmiennych jako�ciowych';
run;
/*Badanie ograniczamy do os�b w wieku od 1 roku do 100 lat:*/
proc means data=work.metabric_m;
var "Age at Diagnosis"n;
run; /*max 96.29 wi�c nie trzeba dalej filtrowa�*/
/*Wyznaczamy pomocnicz� zmienn� jako�ciow� age_c :*/
data work.metabric_m;
set work.metabric_m;
if "Age at Diagnosis"n<34 then age_c=1;
if 34<="Age at Diagnosis"n<45 then age_c=2;
if 45<="Age at Diagnosis"n<57 then age_c=3;
if "Age at Diagnosis"n>=57 then age_c=4;
run;
/* Podgl�d efektu */
proc freq data=work.metabric_m;
  tables age_c / missing;
  title 'Rozk�ad nowej zmiennej age_c (grupy wiekowe)';
run;