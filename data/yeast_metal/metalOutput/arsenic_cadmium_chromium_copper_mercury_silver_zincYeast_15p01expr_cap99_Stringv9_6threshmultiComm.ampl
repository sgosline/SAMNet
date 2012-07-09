model metalOutput/arsenic_cadmium_chromium_copper_mercury_silver_zincYeast_15p01expr_cap99_Stringv9_6threshmultiCommchangeflow.mod;
data metalOutput/arsenic_cadmium_chromium_copper_mercury_silver_zincYeast_15p01expr_cap99_Stringv9_6threshmultiComm.dat;
option solver cplexamp; 
solve; 
printf {(i,j) in all_interactions,k in commodities}: "%s\t%s\t%s\t%f\n", i,j,k,X[i,j,k]>metalOutput/arsenic_cadmium_chromium_copper_mercury_silver_zincYeast_15p01expr_cap99_Stringv9_6threshmultiComm.txt; 
