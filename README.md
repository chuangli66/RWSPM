# Folder code_RWSPM: contains all core functions for imposing “RegionWise Statistical Parametric Mapping”.
   - Function RWSPM in “EWSPM.R” performs regional partitioning, distribution aggregation within each sub-region, and region-wise independence test.
   - Function RW.MTCCT in “EWSPM.R” identifies spatially contiguous regions where multiple adjacent subregions exhibit moderate but collectively meaningful associations.
   - Function golden_section_search in “slide_width.R” is used to select a sliding-region width by minimizing loss function defined in section 2.1 of manuscript.

## Folder compared_method: 
-	The folder contains five subdirectories (code_GBJ, code_AdaMant, code_Rdcov, code_sKPCR, code_fvGWAS), each corresponding to one of the compared methods. 
- 	Code “Y_slide.R” is used for image partitioning and serves for the five compared methods in simulation.

## Folder Simu_setting: 
- This folder stores the simulation setting codes: “sim_s1.py” generates the data under the null hypothesis; “sim_s2.py”, “sim_s3.py”, and “sim_s4.py” generate data for Model 1 of manuscript under different parameter settings. “sim_s5_x2.py”, “sim_s6_x2.py”, and “sim_s7_x2.py” generate data for Model 2 of manuscript under different parameter settings. 

## Folder code_simu_result:
-	Code “DOR.R” generates box plots for six compared methods under settings 2-7, all results are summarized in Figure4 of manuscript.
-	Code “load_sim_result.R” loads and stores the results of the six compared methods under Settings 1-7.
-	Code “sum_sim_result.R” summarizes Type I error rate and power results under settings 1-7, all results are summarized in Table2 of manuscript.

## Folder code_realdata:  
-  Due to data privacy and storage limitations, we did not put the
original genotype data from ADNI online. Instead, a general illustrative set of scripts is included in the code_realdata folder, providing the code necessary to reproduce the real-data analysis 	presented in Section 5 of manuscript.
-	Before the performing the data application, the following input files are required:
1.The genetic data file (in our case, named as ADNI_1_3.RData),
2.The phenotype data files (in our case, named as left_hipp_ori.Rdata and right_hipp_ori.Rdata),
3.The covarites information file (in our case, named as covarites_inf.Rdata),
4.The genetic data information file including chr, SNP, position, Alleles (in our case, named as snp_inf1_3.Rdata).

### Step 1. 
- 	Run “processing_realdata.R” to perform preprocessing of the phenotype data including regional partitioning, adjustment for covariate effects, and other related preparatory procedures.

### Step 2.
-	Run “rwspm rwspm_realdata.R” to impose our proposed method RWSPM on processed genetic-phenotype data.
- 	Run “realdata_sub_region.R” to segment phenotype data and serves for five compared methods. 
- 	Run “realdata_gbj.R”, “realdata_adaMant.R”, “realdata_Rdcov.R”, “realdata_ sKPCR.py”, “realdata_fvGWAS.m” to impose five compared methods on processed genetic-phenotype data. 

### Step 3. 
- 	Run “summarize_result.R” to summarize all reasults in real data analysis. 

## The RWSPM is organized into an R package ( /.../RWSPM_0.1.0.tar.gz).
