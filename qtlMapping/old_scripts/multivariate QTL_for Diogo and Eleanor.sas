
*reads in phenotypes;
PROC IMPORT OUT= mouse
            DATAFILE= "D:\Work\Imprinting project\multivariate analysis\Multivariate\F3Phenotypes_further corrected family data_corrected litter sizes.xlsx"
            DBMS=xlsx REPLACE; 
     GETNAMES=YES; 
quit;


*This just makes sure that we are analyzing only the F3 since some files I use will contain both generations;
data mouse; set mouse;
where gen = 'F3';
run;

*just renaming the variable to make it less confusing!;
data mouse; set mouse;
rename sexAN = sex;
run;

*this creates growth traits from the weights - makes sure that they are all present in the phenotypes since we don't always use them all;
data mouse; set mouse;
grow12 = week2 - week1;
grow23 = week3 - week2;
grow34 = week4 - week3;
grow45 = week5 - week4;
grow56 = week6 - week5;
grow67 = week7 - week6;
grow78 = week8 - week7;
run;

*This takes the wide data set and makes it tall - I redid this using my own approach, not the Chenoweth version I showed you;
data mouse_tall;
set mouse;
value = grow12; trait = 1; output;
value = grow23; trait = 2; output;
value = grow34; trait = 3; output;
value = grow45; trait = 4; output;
value = grow56; trait = 5; output;
value = grow67; trait = 6; output;
value = grow78; trait = 7; output;
run;


PROC SORT; BY ID; quit;

*Reads in the genotypes file - here it will only be used to determine which animals were actually genotyped
so that they can be used in a reduced model;
PROC IMPORT OUT= mouse2
            DATAFILE= "D:\Work\Imprinting project\multivariate analysis\Multivariate\chrom1.XLSx"
			DBMS=xlsx REPLACE;
     GETNAMES=YES;
quit;

PROC SORT; BY ID; quit;

*merges genotypes to phenotypes;
DATA mouse2a; MERGE mouse_tall mouse2; BY ID; run;

*this deletes any missing observations that have no genotype - so that the reduced model 
only contains individuals that have a genotype and contribute to the full model;
data mouse2a; set mouse2a;
if LocusA = "" then delete;
run;

proc sort data = mouse2a; by ID; quit;

*Runs the reduced model that does not contain any genotype information;
proc mixed data = mouse2a covtest lognote IC Method = ML;
class gen sire dam family id trait sex; 
model value = trait trait*sex / noint solution ddfm=satterthwaite;
random trait /sub=family type=un;
repeated trait/sub=id(family) type=un;
ODS output InfoCrit = ICR;
quit;


*extracts the model fit parameters - i.e., -2LogLikelihood and number of parameters;
data ICRa; set ICR;
rename Neg2LogLike = Reduced;
rename Parms = ParmsReduced ;
run;

data ICRa; set ICRa;
keep Reduced ParmsReduced;
run;

*enters the macro to run across the genome;
%MACRO TEST;
%LET Y = 0;

*iterates by chromosome;
%DO L = 1 %TO 19;
%LET S = %EVAL(&L);


    %IF &S = 1 %THEN %DO; %LET N = 31;  %END;
    %IF &S = 2 %THEN %DO; %LET N = 28;  %END;
    %IF &S = 3 %THEN %DO; %LET N = 26;  %END;
    %IF &S = 4 %THEN %DO; %LET N = 24;  %END;
    %IF &S = 5 %THEN %DO; %LET N = 19;  %END;
    %IF &S = 6 %THEN %DO; %LET N = 22;  %END;
    %IF &S = 7 %THEN %DO; %LET N = 18;  %END;
    %IF &S = 8 %THEN %DO; %LET N = 13;  %END;
    %IF &S = 9 %THEN %DO; %LET N = 19;  %END;
    %IF &S = 10 %THEN %DO; %LET N = 18; %END;
    %IF &S = 11 %THEN %DO; %LET N = 21; %END;
    %IF &S = 12 %THEN %DO; %LET N = 19; %END;
    %IF &S = 13 %THEN %DO; %LET N = 18; %END;
    %IF &S = 14 %THEN %DO; %LET N = 14; %END;
    %IF &S = 15 %THEN %DO; %LET N = 16; %END;
    %IF &S = 16 %THEN %DO; %LET N = 13; %END;
    %IF &S = 17 %THEN %DO; %LET N = 13; %END;
    %IF &S = 18 %THEN %DO; %LET N = 13; %END;
    %IF &S = 19 %THEN %DO; %LET N = 8;  %END;

*reads in the chromosome file needed for this run through the loop;
PROC IMPORT OUT= mouse2b
            DATAFILE= "D:\Work\Imprinting project\multivariate analysis\Multivariate\chrom&S..XLS"
			DBMS=EXCEL97 REPLACE;
     GETNAMES=YES;
quit;

PROC SORT data = mouse2b; BY ID; quit;

*merges genotype and phenotype;
DATA mouse3; MERGE mouse_tall mouse2b; BY ID; run;


*now starts the inner loop that runs through markers on the chromosome
N is the number of markers on the chromosome and X will keep track of the marker number each run through, while Y
just counts up the genomic position from 1 to the total number of markers in the genome;
  %LET X = 0;
  %DO I = 1 %TO &N;
  %LET X = %EVAL(&X+1);
  %LET Y = %EVAL(&Y+1);

*assigns the additive, dominance, and if needed, imprinting, values to run for a locus. One could use
  the macro variables here (i.e., A&X, D&X, I&X), but SAS doesn't seem to like to include them in an interaction term, so I just assign the
  macro variables to fixed names;
data mouse3; set mouse3;
LocusA = A&X;
LocusD = D&X;
LocusI = I&X;
run;

*This runs the multivariate model for the locus in question, using the additive and dominance effects;
proc mixed data = mouse3 covtest method = ML IC lognote; 
class sire dam family id trait sex; 
model value = trait trait*sex LocusA*trait LocusD*trait / noint solution ddfm=satterthwaite;
random trait /sub=family type=un;
repeated trait/sub=id(family) type=un;
ODS output solutionF = solutions;
ODS output covparms = covars;
ods output  SolutionF = vectors;
ods output Tests3 = Tests;
ODS output InfoCrit = ICF ;
quit;

*this extracts the model fit for the full model;
data ICFa; set ICF;
rename Neg2LogLike = Full;
rename Parms = ParmsFull ;
run;

data ICFa; set ICFa;
keep Full ParmsFull;
run;

*this merges the model fit info for the full and reduced models;
data ICRF; merge ICFa ICRa;
run;

*this uses the full and reduced model fit info for the likelihood ratio test;
data icrf2; set icrf;
chivalue = Reduced - Full;
degfree = ParmsFull - ParmsReduced ;
mPROB = 1-PROBchi(chivalue, degfree);  
mLOD = LOG10(1/mPROB); 
LOCUS = &X;
CHR = &S;
Count = &Y;
run;

*this appends the LRT info to a data set so it will accumulate the tests as it walks through the genome;
Proc Append Base=Significance data=icrf2 force; quit; 

*This next section goes through the tests generated by the full model and pulls out the significance tests 
for the vectors of effects;

data testsa; set tests;
If Effect = 'LocusA*trait' then merger = 1;
where Effect = 'LocusA*trait';
run;

data testsd; set tests;
If Effect = 'LocusD*trait' then merger = 1;
where Effect = 'LocusD*trait';
run;

data testsa; set testsa; 
rename Fvalue = FvalueA;
rename numDF = numDFA;
rename denDF = denDFA;
run;

data testsd; set testsd; 
rename Fvalue = FvalueD;
rename numDF = numDFD;
rename denDF = denDFD;
run;


Data testsADI; 
Merge  testsa testsd; by merger; run; 

*this uses the F values to generate p values and the p values to calculate the negative log of p (same as log(1/p), 
which we use as the significance test;
data testsFINALADI; set testsADI;
aPROB = 1 - PROBF(FvalueA, numDFA, dendfA);
dPROB = 1 - PROBF(FvalueD, numDFD, dendfD);
aLOD = LOG10(1/aPROB); 
dLOD = LOG10(1/dPROB);
LOCUS = &X;
CHR = &S;
Count = &Y;

Keep  LOCUS CHR Count aLOD dLOD  FvalueA numDFA denDFA FvalueD numDFD denDFD ;

run;

*this appends these vectors to a data set to accumulate tests as it walks through the genome;
Proc Append Base=LocusTESTS data=testsFINALADI; run; 

********************;
data ALLtests; merge significance locustests; by count; run;

*This next section goes through the results from the full model and pulls out the estimates
of the effects of the vectors;

data vectorA; set solutions;
If Effect = 'LocusA*trait' then merger = 1;
where Effect = 'LocusA*trait';
run;

data vectorD; set solutions;
If Effect = 'LocusD*trait' then merger = 1;
where Effect = 'LocusD*trait';
run;

data vectorA; set vectorA; 
rename Estimate = Additive;
rename StdErr = AdditiveSE;
rename DF = AdditiveDF;
rename tvalue = AdditiveT;
run;

data vectorD; set vectorD; 
rename Estimate = Dominance;
rename StdErr = DominanceSE;
rename DF = DominanceDF;
rename tvalue = DominanceT;
run;
= ImprintT;
run;

Data Vectors; 
Merge  VectorA VectorD; by merger; run; 


*this generates the significance tests for the vectors;
data vectors; set vectors;
LOCUS = &X;
CHR = &S;
Count = &Y;
Aprob = 1 - probt(additiveT, additiveDF);
Dprob = 1 - probt(dominanceT, dominanceDF);
aLOD = LOG10(1/Aprob); 
dLOD = LOG10(1/Dprob); 
Keep  LOCUS CHR Count trait 
additive additiveSE additiveDF additiveT ALOD
dominance dominanceSE dominanceDF dominanceT DLOD;
run;

*****************************************************;
*this now appends the vector estimates to a data set to accumulate them as it walks through the genome;
Proc Append Base=vectorsALL data=vectors; run; 


*this next section goes through and keeps the estimates of the residual G matrices -
but I have no idea why they might be useful;
data genetic; set covars;
If Subject = 'FAMILY' then merger = 1;
where Subject = 'FAMILY';
run;

data Envi; set covars;
If Subject = 'ID(FAMILY)' then merger = 1;
where Subject = 'ID(FAMILY)';
run;


data genetic; set genetic; 
rename Estimate = Genetic;
rename StdErr = GeneticSE;
rename ZValue = GeneticZ;
run;

data Envi; set Envi; 
rename Estimate = Envi;
rename StdErr = EnviSE;
rename ZValue = EnviZ;
run;

data matrices; merge genetic envi; by merger; run;

data matrices; set matrices;
LOCUS = &X;
CHR = &S;
Count = &Y;
Keep  LOCUS CHR Count CovParm
Genetic GeneticSE GeneticZ
Envi EnviSE EnviZ;
run;

Proc Append Base=matricesALL data=matrices; run; 

%END;

%END;


*This is the end of the macro;
%MEND TEST;
%TEST
run;

*these steps just output the results by taking the data sets that have been accumulating the analyses and 
saving them as Excel files;
Proc export data=ALLtests outfile="D:\Work\Imprinting project\multivariate analysis\growth\TESTS6 growths_uncorr_no i.xls"
DBMS = excel97 replace;
quit;

Proc export data=vectorsALL outfile="D:\Work\Imprinting project\multivariate analysis\growth\vectors6 growths_uncorr_no i.csv"
DBMS = csv replace;
quit;

Proc export data=matricesALL outfile="D:\Work\Imprinting project\multivariate analysis\growth\matrices6 growths_uncorr_no i.csv"
DBMS = csv replace;
quit;

