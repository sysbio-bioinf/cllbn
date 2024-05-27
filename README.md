# CLL to Richter Syndrome: integrating network-based strategies with experimental readouts elucidating disease-drivers and personalized therapies

# Boolean network simulation and experimental data evaluation

Simulation, analysis and code for reproducing the results of the draft can be found in CLL_FigureScripts.R. The script is ordered according to the figure numbers in the draft. All additional custom functions required to run the script are provided within the Function folder, and loaded automatically within the main script. All additional file sources are provided in the repository and listed below: 
- The modeled Boolean network is given in the text-based BoolNet-format (CLL_Model.txt, 		CLL_Model_Microenvironment.txt)
- Experimental data on the anonymised patient cohort can be found in 								BMI1_Expression_Anonymous.csv and CD5_Expression_Anonymous.csv.
- The input output analysis readouts are provided as CLL_IO-Notch_TP53.RData 
 

# Validation of CLL/Richter Boolean network using Single Cell RNA Seq data

The network predictions were validated using a single cell RNA Seq dataset of human CLL patients by Nadeu et al., Nature Medicine, 2022. For the analysis, we used the publicly available Seurat-objects (https://zenodo.org/records/6631966). 
Code for the validation of the Boolean network prediction based on single cell RNA sequencing data can be found in the RMarkdown-file SingleCellAnalysis.Rmd.

Additional files : 
 - Intermediary results are given as RData-objects (binarizedDataset.RData, toBinarizeRevised.RData). 
 - To compare the scRNA seq data with the network predictions, we used the attractors given in the files RS_AKT_KI_attrmatrix.Rdata, RS_CDKN2AB_TP53_KO_attrmatrix.Rdata, RS_NFAT_KO_attrmatrix.Rdata, RS_NOTCH_KI_attrmatrix.Rdata, CLL_attrmatrix.Rdata)
 - A mapping of network gene names to ENSEMBL-IDs can be found in Attrgene_plus_EnsID.Rdata
 - SingleCellAnalysis.nb.html is a compiled html-file showing the results as generated within  SingleCellAnalysis.Rmd
