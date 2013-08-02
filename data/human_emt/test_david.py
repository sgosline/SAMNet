import subprocess

subprocess.Popen("python ../../src/go_enrichment/DAVID.py --listName=siffile --listType=CYTOSIF --threshold=Pvalue --thresholdValue=0.1 --outputDir=./goOutput ./emtOutput/fixed_snail_tgfb_zeb1_network_gamma14_cap_99_Psiquic_5threshmultiComm_mcfs.sif", shell=True)
