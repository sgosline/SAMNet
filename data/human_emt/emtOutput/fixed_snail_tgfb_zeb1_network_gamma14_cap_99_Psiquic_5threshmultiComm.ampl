model emtOutput/fixed_snail_tgfb_zeb1_network_gamma14_cap_99_Psiquic_5threshmultiCommchangeflow.mod;
data emtOutput/fixed_snail_tgfb_zeb1_network_gamma14_cap_99_Psiquic_5threshmultiComm.dat;
option solver cplexamp; 
solve; 
printf {(i,j) in all_interactions,k in commodities}: "%s\t%s\t%s\t%f\n", i,j,k,X[i,j,k]>emtOutput/fixed_snail_tgfb_zeb1_network_gamma14_cap_99_Psiquic_5threshmultiComm.txt; 
