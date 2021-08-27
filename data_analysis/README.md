# Data analyses
R codes used for data analyses. 

All numbers and figures (Table 3 and Figures 2, 3, 4) presented in the data analyses can be reproduced (with minor differences, because a synthetic dataset is used instead of the real NHID-2010 data) using codes in this folder.

Run `data_analysis_master.R` to run the complete flow of data analyses. Input and output files are noted in comments.

Intermediate results (`synthetic_harmriskfactor.rda`, `synthetic_landmarktime.rda`, and `synthetic_gsprofile.rda`), and final results (`figure2a.jpeg`, `figure2b.jpeg`, `figure3.jpeg`, `figure4a.jpeg`, `figure4b.jpeg`, `figure4c.jpeg`, `figure4d.jpeg`, `figure4e.jpeg`, `figure4f.jpeg`, and `table3.csv`) are stored in subfolder `data_analysis_result`. Move relevant rda files to `data_analysis` folder to run specific parts of the data analysis without running the prerequisites.