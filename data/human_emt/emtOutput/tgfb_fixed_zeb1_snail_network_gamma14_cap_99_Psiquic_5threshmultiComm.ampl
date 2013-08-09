model emtOutput/tgfb_fixed_zeb1_snail_network_gamma14_cap_99_Psiquic_5threshmultiCommchangeflow.mod;
data emtOutput/tgfb_fixed_zeb1_snail_network_gamma14_cap_99_Psiquic_5threshmultiComm.dat;
option solver cplexamp; 
solve; 
printf {(i,j) in all_interactions,k in commodities}: "%s\t%s\t%s\t%f\n", i,j,k,X[i,j,k]>emtOutput/tgfb_fixed_zeb1_snail_network_gamma14_cap_99_Psiquic_5threshmultiComm.txt; 
