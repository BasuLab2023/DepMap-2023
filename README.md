Manuscript: 

# A novel pipeline for prioritizing cancer type-specific therapeutic vulnerabilities using DepMap identifies PAK2 as a target in head and neck squamous cell carcinomas

Malay K. Sannigrahi* 1; Austin C. Cao* 1; Pavithra Rajagopalan1; Lova Sun2, Robert M. Brody1; Lovely Raghav1; Phyllis A. Gimotty3; and Devraj Basu1,4

*authors contributed equally

1 Department of Otorhinolaryngology-Head and Neck Surgery, University of Pennsylvania, Philadelphia, PA
2 Department of Medicine, University of Pennsylvania, Philadelphia, PA
3 Department of Biostatistics, Epidemiology and Informatics, University of Pennsylvania, Philadelphia, PA
4 Ellen and Ronald Caplan Cancer Center, The Wistar Institute, Philadelphia, PA

Running Title: DepMap reveals PAK2 is a target in HPV(-) HNSCC

Methodology for Cell line data analysis

Gene dependency data processed by the Chronos algorithm was extracted initially from the 22Q1 public data release from the DepMap at the Broad Institute for all HPV(-) HNSCC (n=63) and esophageal squamous cell carcinoma (ESCC) (n=24) cell lines in this resource. A re-analysis prior to publication was performed with the 23Q2 release. The Drug Gene Interaction Database (DGIdb) was used to filter for the genes predicted to encode for druggable proteins. The Open Targets Platform was used to identify approved or investigational drugs known to target druggable proteins and describe their application in clinical trials to date along with their FDA approval status. Drug response information was sourced from the GDSC2 dataset of the Genomics of Drug Sensitivity in Cancer database (GDSC, Release 8.3). Additional drug response information and RNAi screen results were obtained from DepMap PRISM Repurposing Public 23Q2 dataset and Gene Effect RNAi (DEMETER2) data, respectively, via the DepMap data portal. 
The R code used to prioritize dependencies and visualize them as a dot plot showing the number of cell lines with each dependency vs. median gene effect score in those cell lines. Chronos Gene effect score and gene probability score of cell-lines can be downloaded from https://depmap.org/portal/.
