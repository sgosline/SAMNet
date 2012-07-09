set proteins;
set source;
set sink;
set commodities;
set initnodes = source union proteins;
set endnodes = sink union proteins;
set sourcesink = source union sink;
set nodes = sourcesink union proteins;
set interactions within {proteins cross proteins};
set all_interactions within {nodes cross nodes};
set sink_interactions{commodities} within {i in proteins,j in sink};
set source_interactions{commodities} within {i in source,j in proteins};
param cost {all_interactions,commodities} >=0;
param capacity {all_interactions} >=0;
var X {all_interactions,commodities} >=0;
minimize Total_Cost:
sum{k in commodities}sum{(i,j) in interactions}-log(cost[i,j,k])*X[i,j,k] + sum{k in commodities}sum{(i,j) in source_interactions[k]}-log(cost[i,j,k])*X[i,j,k] + sum{k in commodities}sum{(i,j) in sink_interactions[k]}-log(cost[i,j,k])*X[i,j,k] - sum{k in commodities}sum{(i,j) in source_interactions[k]}15*X[i,j,k];
subject to Kirkoff {k in proteins,c in commodities}: sum {(i,k) in all_interactions} X[i,k,c]=sum{(k,j) in all_interactions} X[k,j,c];
subject to sourcesinkcond{k in commodities}: sum{(i,j) in source_interactions[k]} X[i,j,k]=sum{(i,j) in sink_interactions[k]} X[i,j,k];
subject to multi {(i,j) in all_interactions}: sum {k in commodities}X[i,j,k]<=capacity[i,j];