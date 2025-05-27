proc import datafile="C:\Users\laura\Documents\GitHub\ACT-projekt\Breast Cancer METABRIC.csv"
out=work.metabric_m
dbms=csv
replace;

getnames=yes;
guessingrows=max;
run;


options encoding='utf-8';
/* KROK 1: Filtrowanie – zostawiamy tylko przypadki raka piersi         */
/************************************************************************/
data work.metabric_m;
/*  set work."BREAST CANCER METABRIC"n;*/
  set work.metabric_m;
  if "Cancer Type"n = 'Breast Cancer';
run;
/* KROK 1.5: Zmieniamy wszystkie spacje w nazwach kolumn na znak pod³ogi "_" */
/********************************************************************************/
proc sql noprint;
    select 
        cats('"', name, '"n = ',
             substr(
               case 
                 when prxmatch('/^[0-9]/', name) = 1 
                 then cats('var_',
                    prxchange('s/_+$//', 1, 
                      prxchange('s/_+/_/o', -1,
                        prxchange('s/[^a-zA-Z0-9]/_/o', -1, name))))
                 else 
                    prxchange('s/_+$//', 1,
                      prxchange('s/_+/_/o', -1,
                        prxchange('s/[^a-zA-Z0-9]/_/o', -1, name)))
               end,
             1, 32)
        )
    into :renamelist separated by ' '
    from dictionary.columns
    where upcase(libname)='WORK' and upcase(memname)='METABRIC_M';
quit;

data work.metabric_m;
    set work.metabric_m(rename=(&renamelist));
run;



/* KROK 2: Usuniêcie kolumny 'Cancer Type', bo ju¿ nie bêdzie potrzebna */
/************************************************************************/
data work.metabric_m;
  set work.metabric_m(drop=Cancer_Type);
run;
/* KROK 3: Sprawdzenie dostêpnych kategorii w kolumnie Vital Status     */
/************************************************************************/
proc freq data=work.metabric_m;
  tables Patient_s_Vital_Status / missing;
  title "Kategorie w Patient's Vital Status";
run;

/* KROK 4: Zostawiamy tylko Living oraz Died of Disease                 */
/************************************************************************/
data work.metabric_m;
  set work.metabric_m;
  if Patient_s_Vital_Status in ('Living', 'Died of Disease');
run;

/* KROK 5: Tworzymy zmienn¹ event                                       */
/************************************************************************/
data work.metabric_m;
  set work.metabric_m;
  if Patient_s_Vital_Status = 'Died of Disease' then event = 1;
  else if Patient_s_Vital_Status = 'Living'      then event = 0;
run;

proc freq data=work.metabric_m;
  tables 
    Cohort
    _character_
    _numeric_
    / missing noprint out=freqs(drop=percent cumfreq cumpercent);
run;

/* Podstawowe statystyki i liczba brak w dla wszystkich numerycznych */
/************************************************************************/
proc means data=work.metabric_m 
           n nmiss mean std min p25 median p75 max;
  var _numeric_;
  title 'Statystyki opisowe i liczba brak w dla zmiennych numerycznych';
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
    var_3_Gene_classifier_subtype
    Cancer_Type_Detailed
    ER_status_measured_by_IHC
    Pam50_Claudin_low_subtype
    Patient_s_Vital_Status
	ER_Status
	HER2_Status
	HER2_status_measured_by_SNP6
	Hormone_Therapy
	Inferred_Menopausal_State
	Integrative_Cluster
	Oncotree_Code
	Overall_Survival_Status
	PR_Status
	Primary_Tumor_Laterality
	Radio_Therapy
	Relapse_Free_Status
	Tumor_Other_Histologic_Subtype
	Type_of_Breast_Surgery
  / missing;  
  title 'Rozklady zmiennych jakosciowych';
run;
/*Badanie ograniczamy do os b w wieku od 1 roku do 100 lat:*/
proc means data=work.metabric_m;
var Age_at_Diagnosis;
run; /*max 96.29 wi c nie trzeba dalej filtrowa */


/*Age at diagnosis*/
proc means data=work.metabric_m Q1 median Q3;
    var Age_at_Diagnosis;
run;

proc sgplot data=work.metabric_m;
    histogram Age_at_Diagnosis;
    title "Histogram wieku w momencie diagnozy";
run;

proc sgplot data=work.metabric_m;
    histogram Tumor_Size;
run;


/*----------------------------------------------------------*/
/*----------------------------------------------------------*/
/*Wyznaczamy pomocnicza zmienna jakosciowa age_c :*/
data work.metabric_m;
set work.metabric_m;
if Age_at_Diagnosis<49 then age_c=1;
if 49<=Age_at_Diagnosis<59 then age_c=2;
if 59<=Age_at_Diagnosis<68 then age_c=3;
if Age_at_Diagnosis>=68 then age_c=4;
run;


data work.metabric_m;
set work.metabric_m;
if Tumor_Size=0 then ts_c=0;
if 0<Tumor_Size<5 then ts_c=1;
if 5<=Tumor_Size<10 then ts_c=2;
if 10<=Tumor_Size<20 then ts_c=3;
if 20<=Tumor_Size<50 then ts_c=4;
if 50<=Tumor_Size<80 then ts_c=5;
if 80<=Tumor_Size then ts_c=6;
run;


data work.metabric_m;
set work.metabric_m;
if Nottingham_prognostic_index <1 then npi_c=1;
if 1<=Nottingham_prognostic_index<2 then npi_c=2;
if 2<=Nottingham_prognostic_index<3 then npi_c=3;
if 3<=Nottingham_prognostic_index<4 then npi_c=4;
if 4<=Nottingham_prognostic_index<5 then npi_c=5;
if 5<=Nottingham_prognostic_index<6 then npi_c=6;
if 6<=Nottingham_prognostic_index then npi_c=7;
run;



data work.metabric_m;
set work.metabric_m;
if Lymph_nodes_examined_positive=0 then lnep_c=0;
if 0<Lymph_nodes_examined_positive<=2 then lnep_c=1;
if 2<Lymph_nodes_examined_positive<=6 then lnep_c=2;
if 6<Lymph_nodes_examined_positive<=15 then lnep_c=3;
if 15<Lymph_nodes_examined_positive then lnep_c=4;
run;


/*data work.metabric_m;*/
/*set work.metabric_m;*/
/*if Mutation_Count=0 then lnep_c=0;*/
/*if 0<Mutation_Count<=2 then lnep_c=1;*/
/*if 2<Mutation_Count<=6 then lnep_c=2;*/
/*if 6<Mutation_Count<=15 then lnep_c=3;*/
/*if 15<Mutation_Count then lnep_c=4;*/
/*run;*/



/*Age at diagnosis*/
proc means data=work.metabric_m Q1 median Q3;
    var Age_at_Diagnosis;
run;

proc sgplot data=work.metabric_m;
    histogram Age_at_Diagnosis;
    title "Histogram wieku w momencie diagnozy";
run;






/* Podglad efektu */
proc freq data=work.metabric_m;
  tables age_c / missing;
  title 'Rozklad nowej zmiennej age_c (grupy wiekowe)';
run;

data work.metabric_m;
  set work.metabric_m(rename=(
    Overall_Survival_Months = t
    event = c
  ));
run;


/*proc lifetest data=work.metabric_m plots=survival(cl);*/
/*   time t*c(0);*/
/*   strata age_c;*/
/*   title "Krzywe prze¿ycia Kaplana-Meiera wed³ug grup wiekowych";*/
/*run;*/




/************************************************************************/
/************************************************************************/
/*MODEL */

proc phreg data=work.metabric_m;
	class Cancer_Type_Detailed
	Cellularity
	Chemotherapy
	ER_Status
	HER2_Status
	Hormone_Therapy
	Inferred_Menopausal_State
	Integrative_Cluster
	Oncotree_Code
	PR_Status
	Pam50_Claudin_low_subtype
	Primary_Tumor_Laterality
	Radio_Therapy
	Tumor_Other_Histologic_Subtype
	Type_of_Breast_Surgery
	var_3_Gene_classifier_subtype
	age_c
	ts_c
	npi_c
	lnep_c
	Tumor_Stage
	Cohort
	Neoplasm_Histologic_Grade;
	model t*c(0)= Cancer_Type_Detailed
	Cellularity
	Chemotherapy
	ER_Status
	HER2_Status
	Hormone_Therapy
	Inferred_Menopausal_State
	Integrative_Cluster
	Oncotree_Code
	PR_Status
	Pam50_Claudin_low_subtype
	Primary_Tumor_Laterality
	Radio_Therapy
	Tumor_Other_Histologic_Subtype
	Type_of_Breast_Surgery
	var_3_Gene_classifier_subtype
	age_c
	Tumor_Stage
	Cohort
	Mutation_Count
	Neoplasm_Histologic_Grade
	ts_c
	npi_c
	lnep_c
	/ties=efron selection=stepwise;
run;


/*data work.metabric_selected;*/
/*    set work.metabric_m(keep=*/
/*	Nottingham_prognostic_index*/
/*	age_c*/
/*	var_3_Gene_classifier_subtype*/
/*	Radio_Therapy*/
/*	Tumor_Size*/
/*	Hormone_Therapy*/
/*	HER2_Status*/
/*	Neoplasm_Histologic_Grade*/
/*	t c);*/
/*run;*/

data work.metabric_selected;
    set work.metabric_m(keep=
	lnep_c
	var_3_Gene_classifier_subtype
	age_c
	Radio_Therapy
	ts_c
	Hormone_Therapy
	HER2_Status
	Primary_Tumor_Laterality
	t c);
run;


/*Weryfikacja za³o¿enia proporcjonalnych hazardów – metoda graficzna*/

proc phreg data=work.metabric_selected;
	model t*c(0)= lnep_c / ties=efron;
	strata lnep_c;
	baseline out=zb_lls_lnep loglogs=lls / method=pl;
	run;
proc gplot data=zb_lls_lnep;
	plot lls*t=lnep_c;
	symbol1 I=join color=blue line=1 value=none;
	symbol2 I=join color=red line=1 value=none;
	symbol3 I=join color=green line=1 value=none;
	symbol4 I=join color=orange line=1 value=none;
	run;



proc phreg data=work.metabric_selected;
	model t*c(0)= age_c / ties=efron;
	strata age_c;
	baseline out=zb_lls_age loglogs=lls / method=pl;
	run;
proc gplot data=zb_lls_age;
	plot lls*t=age_c;
	symbol1 I=join color=blue line=1 value=none;
	symbol2 I=join color=red line=1 value=none;
	symbol3 I=join color=green line=1 value=none;
	symbol4 I=join color=orange line=1 value=none;
	run;



proc phreg data=work.metabric_selected;
	class var_3_Gene_classifier_subtype;
	model t*c(0)= var_3_Gene_classifier_subtype / ties=efron;
	strata var_3_Gene_classifier_subtype;
	baseline out=zb_lls_var3 loglogs=lls / method=pl;
	run;
proc gplot data=zb_lls_var3 ;
	plot lls*t=var_3_Gene_classifier_subtype;
	symbol1 I=join color=blue line=1 value=none;
	symbol2 I=join color=red line=1 value=none;
	symbol3 I=join color=green line=1 value=none;
	symbol4 I=join color=orange line=1 value=none;
	run;



proc phreg data=work.metabric_selected;
	class Radio_Therapy;
	model t*c(0)= Radio_Therapy / ties=efron;
	strata Radio_Therapy;
	baseline out=zb_lls_rt loglogs=lls / method=pl;
	run;
proc gplot data=zb_lls_rt ;
	plot lls*t=Radio_Therapy;
	symbol1 I=join color=blue line=1 value=none;
	symbol2 I=join color=red line=1 value=none;
	symbol3 I=join color=green line=1 value=none;
	symbol4 I=join color=orange line=1 value=none;
	run;



proc phreg data=work.metabric_selected;
	model t*c(0)= ts_c / ties=efron;
	strata ts_c;
	baseline out=zb_lls_ts loglogs=lls / method=pl;
	run;
proc gplot data=zb_lls_ts;
	plot lls*t=ts_c;
	symbol1 I=join color=blue line=1 value=none;
	symbol2 I=join color=red line=1 value=none;
	symbol3 I=join color=green line=1 value=none;
	symbol4 I=join color=orange line=1 value=none;
	symbol5 I=join color=green line=1 value=none;
	symbol6 I=join color=orange line=1 value=none;
	run;



proc phreg data=work.metabric_selected;
	class Hormone_Therapy;
	model t*c(0)= Hormone_Therapy / ties=efron;
	strata Hormone_Therapy;
	baseline out=zb_lls_ht loglogs=lls / method=pl;
	run;
proc gplot data=zb_lls_ht;
	plot lls*t=Hormone_Therapy;
	symbol1 I=join color=blue line=1 value=none;
	symbol2 I=join color=red line=1 value=none;
	run;



proc phreg data=work.metabric_selected;
/*	class HER2_Status;*/
	model t*c(0)= HER2_Status / ties=efron;
	strata HER2_Status;
	baseline out=zb_lls_hs loglogs=lls / method=pl;
	run;
proc gplot data=zb_lls_hs;
	plot lls*t=HER2_Status;
	symbol1 I=join color=blue line=1 value=none;
	symbol2 I=join color=red line=1 value=none;
	run;



/*Model tylko ze zmienn¹ age*/
proc phreg data=work.metabric_selected;
	class age_c;
	model t*c(0)= age_c / ties=efron;
	output out=R_Sch_age
	(keep = t  age_c  age_RS) 
	ressch= age_RS;
	run;

goptions reset=all;
goptions htext=1.5 ;
option nodate nonumber;
axis1 order=(-1 0 1)
label= ( angle=90 'Schoenfeld residuals');
axis2 order=(0 73 146 219 292 365)
label= ('t');
legend1 label=none;
proc gplot data=R_Sch_age;
plot age_RS*t / vaxis=axis1  haxis=axis2;
symbol1 v=point i=sm90s width=1 c=blue;
run;



/*Pelny model*/
proc phreg data=work.metabric_selected;
	class 
	HER2_Status
	Hormone_Therapy
	Primary_Tumor_Laterality
	Radio_Therapy
	var_3_Gene_classifier_subtype
	age_c
	ts_c
	lnep_c;
	model t*c(0)=
	HER2_Status
	Hormone_Therapy
	Primary_Tumor_Laterality
	Radio_Therapy
	var_3_Gene_classifier_subtype
	age_c
	ts_c
	lnep_c
	/ties=efron;
	output out = R_Sch_full
	ressch= full_RS;
run;


goptions reset=all;
goptions htext=1.5 ;
option nodate nonumber;
axis1 order=(-1 0 1)
label= ( angle=90 'Schoenfeld residuals');
axis2 order=(0 73 146 219 292 365)
label= ('t');
legend1 label=none;
proc gplot data=R_Sch_full;
plot full_RS*t / vaxis=axis1  haxis=axis2;
symbol1 v=point i=sm90s width=1 c=blue;
run;



data R_Sch_full_Ft;
set R_Sch_full;
lt=log(t);
t2=t**2;
run;
proc corr data=R_Sch_full_Ft;
var t lt t2 full_RS;
run;





proc phreg data=dane.covid_model;
 model t*c(0)=Sex Intubed Pneumonia  Diabetes Copd Inmsupr Hypertension Other_disease Obesity 
Renal_chronic Tobacco Contact_other_covid age_c
 Sex_t Intubed_t Pneumonia_t Diabetes_t Copd_t Inmsupr_t Hypertension_t Other_disease_t Obesity_t
 Renal_chronic_t Tobacco_t Contact_other_covid_t age_c_t
 /ties=efron;
 Sex_t=Sex*t;
 Intubed_t=Intubed*t; 
Pneumonia_t=Pneumonia*t;
 Diabetes_t=Diabetes*t;
 Copd_t=Copd*t;
 Inmsupr_t=Inmsupr*t;
Hypertension_t=Hypertension*t;
 Other_disease_t= Other_disease*t;
 Obesity_t=Obesity*t;
 Renal_chronic_t= Renal_chronic*t;
 Tobacco_t=Tobacco*t;
 Contact_other_covid_t= Contact_other_covid*t;
 age_c_t= age_c*t;
 run;



proc phreg data=work.metabric_selected;
	model t*c(0)=
	HER2_Status
	Hormone_Therapy
	Primary_Tumor_Laterality
	Radio_Therapy
	var_3_Gene_classifier_subtype
	age_c
	ts_c
	lnep_c

	HER2_Status_t
	Hormone_Therapy_t
	Primary_Tumor_Laterality_t
	Radio_Therapy_t
	var_3_Gene_classifier_subtype_t
	age_c_t
	ts_c_t
	lnep_c_t

	/ties=efron;
	HER2_Status_t=HER2_Status*t
	Hormone_Therapy_t=Hormone_Therapy*t
	Primary_Tumor_Laterality_t=Primary_Tumor_Laterality*t
	Radio_Therapy_t=Radio_Therapy*t
	var_3_Gene_classifier_subtype_t=var_3_Gene_classifier_subtype*t
	age_c_t=age_c*t
	ts_c_t=ts_c*t
	lnep_c_t=lnep_c*t;
run;

























/*proc phreg data=work.metabric_m;*/
/*	class age_c;*/
/*	model t*c(0)= age_c/ ties=efron;*/
/*run;*/
/**/
/**/
/*proc phreg data=work.metabric_m;*/
/*	model t*c(0)= age_c / ties=efron;*/
/*	strata age_c;*/
/*	baseline out=zb_lls_sex loglogs=lls / method=pl;*/
/*run;*/
/**/
/*proc gplot data=zb_lls_sex ;*/
/*	plot lls*t=sex;*/
/*	symbol1 I=join color=blue line=1 value=none;*/
/*	symbol2 I=join color=red line=1 value=none;*/
/*run;*/

















