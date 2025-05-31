proc datasets library=work kill nolist;
quit;


/*Import przygotowanych zmiennych*/
proc import datafile="C:\Users\laura\Documents\GitHub\ACT-projekt\data\metabric_cleaned.csv"
out=work.metabric_cleaned
dbms=csv
replace;

getnames=yes;
guessingrows=max;
run;

/*options encoding='utf-8';*/

/*MODEL SEMIPARAMETRYCZNY: W UJECIU KLASYCZNYM*/
/*Selekcja zmiennych*/
proc phreg data=work.metabric_cleaned;
	class type_of_breast_surgery
	cancer_type_detailed
	cellularity
	chemotherapy
	pam_claudin_low_subtype
	cohort
	er_status
	neoplasm_histologic_grade
	her_status
	tumor_other_histologic_subtype
	hormone_therapy
	inferred_menopausal_state
	integrative_cluster
	primary_tumor_laterality
	oncotree_code
	pr_status
	radio_therapy
	sex
	gene_classifier_subtype
	tumor_stage
	age_c
	ts_c
	npi_c
	lnep_c;
	model t*c(0)= type_of_breast_surgery
	cancer_type_detailed
	cellularity
	chemotherapy
	pam_claudin_low_subtype
	cohort
	er_status
	neoplasm_histologic_grade
	her_status
	tumor_other_histologic_subtype
	hormone_therapy
	inferred_menopausal_state
	integrative_cluster
	primary_tumor_laterality
	mutation_count
	oncotree_code
	pr_status
	radio_therapy
	sex
	gene_classifier_subtype
	tumor_stage
	age_c
	ts_c
	npi_c
	lnep_c
	/ties=efron selection=stepwise;
run;

%let vars = lnep_c gene_classifier_subtype age_c radio_therapy ts_c hormone_therapy her_status;

/*Weryfikacja zalozenia proporcjonalnych hazardow: metoda graficzna*/
%macro check_ph_simple(varlist=, data=WORK.METABRIC_CLEANED);

  %let n=%sysfunc(countw(&varlist));

  %do i=1 %to &n;
    %let var=%scan(&varlist, &i);

    proc phreg data=&data;
      model t*c(0) = &var / ties=efron;
      strata &var;
      baseline out=zb_lls_&var loglogs=lls / method=pl;
    run;

    proc gplot data=zb_lls_&var;
      plot lls*t = &var;
      symbol1 i=join color=blue    line=1 value=none;
      symbol2 i=join color=red     line=1 value=none;
      symbol3 i=join color=green   line=1 value=none;
      symbol4 i=join color=purple  line=1 value=none;
      symbol5 i=join color=orange  line=1 value=none;
      symbol6 i=join color=cyan    line=1 value=none;
      symbol7 i=join color=black   line=1 value=none;
      symbol8 i=join color=gray    line=1 value=none;
      symbol9 i=join color=brown   line=1 value=none;
      symbol10 i=join color=magenta line=1 value=none;
    run;

  %end;

%mend check_ph_simple;

/*%check_ph_simple(varlist=lnep_c gene_classifier_subtype age_c radio_therapy ts_c hormone_therapy her_status);*/
%check_ph_simple(varlist=&vars);
/*Zalozenie proporcjonalnych hazardow spelnione dla: ts_c i her_status*/

/*Weryfikacja zalozenia proporcjonalnych hazardow: reszty Schoenfelda*/
/* Makro generujace liste reszt Schoenfelda */
%macro schoenfeld_list(vars);
  %local i var n result;
  %let n = %sysfunc(countw(&vars));
  %let result=;
  %do i=1 %to &n;
    %let var = %scan(&vars, &i);
    %let result = &result &var._RS;
  %end;
  &result
%mend;

/* Uzycie makra do generowania listy */
proc phreg data=work.metabric_cleaned;
  class &vars;
  model t*c(0) = &vars / ties=efron;
  output out=R_Sch_s ressch=%schoenfeld_list(&vars);
run;


goptions reset=all;
goptions htext=1.5;
option nodate nonumber;

axis1 order=(-1 0 1)
  label=(angle=90 'Schoenfeld residuals');
/*axis2 order=(0 20 40 60 80 100)*/
axis2 order=(0 67 135 202 270 337)
  label=('t');
legend1 label=none;

/* Makro generujace wykresy dla kazdej zmiennej */
%macro plot_schoenfeld(vars);
  %local i var n;
  %let n = %sysfunc(countw(&vars));

  %do i=1 %to &n;
    %let var = %scan(&vars, &i);

    proc gplot data=R_Sch_s;
      plot &var._RS*t / vaxis=axis1 haxis=axis2;
      symbol1 v=point i=sm90s width=1 c=blue;
      title "Schoenfeld Residuals for &var";
    run;
    quit;
  %end;
%mend;

%plot_schoenfeld(&vars);
/*Zalozenie proporcjonalnych hazardow spelnione dla: gene_classifier_subtype i radio_therapy*/


data R_Sch_s_Ft;
  set R_Sch_s;
  lt = log(t);     /* logarytm czasu */
  t2 = t**2;       /* czas do kwadratu */
run;

%macro corr_schoenfeld(vars);
  %local i var n;
  %let n = %sysfunc(countw(&vars));

  %do i=1 %to &n;
    %let var = %scan(&vars, &i);

    title "Correlation of &var._RS with time functions";
    proc corr data=R_Sch_s_Ft nosimple;
      var t lt t2;
      with &var._RS;
    run;
  %end;
%mend;

%corr_schoenfeld(&vars);

/* Model Coxa z interakcj¹ z czasem */
proc phreg data=work.metabric_cleaned;
  model t*c(0) = lnep_c gene_classifier_subtype age_c radio_therapy ts_c hormone_therapy her_status
	lnep_c_t gene_classifier_subtype_t age_c_t radio_therapy_t ts_c_t hormone_therapy_t her_status_t
	/ ties=efron;
	lnep_c_t=lnep_c*t;
	gene_classifier_subtype_t=gene_classifier_subtype*t;
	age_c_t=age_c*t;
	radio_therapy_t=radio_therapy*t;
	ts_c_t=ts_c*t;
	hormone_therapy_t=hormone_therapy*t;
	her_status_t=her_status*t;
run;

proc phreg data=work.metabric_cleaned;
	model t*c(0) = lnep_c age_c ts_c hormone_therapy her_status
	/ ties=efron;
	output out=work.diagnostic_out
	xbeta = Xb         /* predyktor liniowy (log-hazard) */
	resmart = Mart     /* reszty martynga³owe */
	resdev = Dev;      /* reszty odchylenia */
run;


proc sgplot data=work.diagnostic_out;
yaxis grid;
refline 0 / axis=y;
scatter y=Mart x=Xb;
title "Reszty Martyngalowe";
run;


proc sgplot data=work.diagnostic_out;
yaxis grid;
refline 0 / axis=y;
scatter y=Dev x=Xb;
title "Reszty Odchylenia";
run;


/*Usuwanie obserwacji odstajacych*/
data work.diagnostic_out_;
set work.diagnostic_out;
if Mart<-2 then delete;
if-2>Dev or 3<Dev  then delete;
run;


/*Wyniki modelu przed usunieciem obserwacji odstajacych*/
proc phreg data=work.metabric_cleaned;
	model t*c(0) = lnep_c gene_classifier_subtype age_c radio_therapy ts_c hormone_therapy her_status
	/ ties=efron;
run;


/*Wyniki modelu po usunieciu obserwacji odstajacych*/
proc phreg data=work.diagnostic_out_;
	model t*c(0) = lnep_c gene_classifier_subtype age_c radio_therapy ts_c hormone_therapy her_status
	/ ties=efron;
run;
/*Wyniki niewiele lepsze po usunieciu skrajnych obserwacji, zostawiamy oryginalne dane*/



/*Analiza wynikow z modelu*/

/*Krzywa przezycia*/
/*Wersja z 7 zmiennymi*/
proc phreg data=work.metabric_cleaned plots=survival;
	class lnep_c gene_classifier_subtype age_c radio_therapy ts_c hormone_therapy her_status;
	model t*c(0) = lnep_c gene_classifier_subtype age_c radio_therapy ts_c hormone_therapy her_status
	/ ties=efron;
	baseline covariates=work.metabric_cleaned out=work.metabric_pred_sur survival=_all_ / diradj;
run;

/*Wersja z 5 zmiennymi*/
proc phreg data=work.metabric_cleaned plots=survival;
	class lnep_c age_c ts_c hormone_therapy her_status;
	model t*c(0) = lnep_c age_c ts_c hormone_therapy her_status
	/ ties=efron;
	baseline covariates=work.metabric_cleaned out=work.metabric_pred_sur survival=_all_ / diradj;
run;



/*Krzywa AUC*/
/*model na 7 zmiennych*/
proc phreg data=work.metabric_cleaned plots(overlay=individual)=roc rocoptions(at= 67 135 202 270 337);
	class lnep_c gene_classifier_subtype age_c radio_therapy ts_c hormone_therapy her_status;
	model t*c(0) = lnep_c gene_classifier_subtype age_c radio_therapy ts_c hormone_therapy her_status
	/ ties=efron;
run;
/*gorszy*/

/*model na 5 zmiennych*/
proc phreg data=work.metabric_cleaned plots(overlay=individual)=roc rocoptions(at= 67 135 202 270 337);
	class lnep_c age_c radio_therapy ts_c her_status;
	model t*c(0) = lnep_c age_c radio_therapy ts_c her_status / ties=efron;
run;


/*Wartosc AUC w zaleznosci od czasu*/
/*model na 7 zmiennych*/
proc phreg data=work.metabric_cleaned plots=auc rocoptions(method=ipcw(cl seed=1234) iauc);
	class lnep_c gene_classifier_subtype age_c radio_therapy ts_c hormone_therapy her_status;
	model t*c(0) = lnep_c gene_classifier_subtype age_c radio_therapy ts_c hormone_therapy her_status
	/ ties=efron;
run;

/*model na 5 zmiennych*/
proc phreg data=work.metabric_cleaned plots=auc rocoptions(method=ipcw(cl seed=1234) iauc);
	class lnep_c age_c radio_therapy ts_c her_status;
	model t*c(0) = lnep_c age_c radio_therapy ts_c her_status / ties=efron;
run;




