# GERMLINE ANALYSIS PIPELINE - Rare Genetic diseases and Familial cancer predisposition

**Author:** Nagarajan Paramasivam, DKFZ

Pipeline checks for the rare, protein function affecting variants in the germline. It also generates a cohort summary for all the samples analyzed in a run.

- **INPUT DATA:**
   - SNVs and Indel varaiants in VCF format produced by the germline platypus pipeline with gene annotation and population frequencies already annotated.
   - Basic information file. Create a text file with some basic information PID, FAMILY_ID, CASE/CONTROL are mandatory. If it is a TRIO analysis, program gets the mother and father info from this file.
    Use the below format

PID | FAMILY_ID | CASE_CONTROL | GENDER | MOTHER_ID | FATHER_ID | PROJECT
--- | --- | --- | --- | --- | --- | ---
PID_01 | FAMILy_1 |  CASE | MALE | PID_02 | PID_03 | SAMPLE
PID_02 | FAMILY_1 |  CONTROL| FEMALE |  NA | NA | SAMPLE
PID_03 | FAMILY_1 | CASE | MALE | NA | NA |  SAMPLE
PID_04 | FAMILY_2 | CASE | NA|  NA|  NA |  SAMPLE

- **Family_Structure:**
   - Trio           = patient child,mother and father
   - Single patient = only case
   - Multiple patient = If there are more than one child in the family, they are considered as separate trios a cohort summary step is performed to merge them.
  
- **Additional Annotation:**
   - Known disease association of variant and the gene.
   - Please download and add the [FunSeq2](http://funseq2.gersteinlab.org/) tools in the current dir for the non-coding annotations.
   - [CADD](https://cadd.gs.washington.edu/score), download the precomputed SNVs and Indels files, better to install the CADD tool locally for novel indels. 
   - [dbNSFP v2.9](https://sites.google.com/site/jpopgen/dbNSFP), download the precomputed files.
   - [GTEx](https://www.gtexportal.org/home/) RPKM values can be downloaded and the tissue specific expression of genes can be calculated. If 80% of the samples for a given tissues have a gene's expression above 75th percentile, then the gene is considered to have high expression in that tissue.

- **TO RUN**
   - The workflow runs on the LSF environment.
   - One sample, only use the patient in case of trio 
`sh run_germlineAnalaysisPipeline.sh -i PID`
   - Multiple patients, use the patient ids in case of families (with parent information) having multiple patients
`sh run_germlineAnalaysisPipeline.sh -i PID_1 PID_2
   - Familial pedigrees, use all the case and control names from the cohort to generate a family and cohort specific analysis.    

- **OUTPUT DATA:**
   - The pipeline generates a list of variants filtered with different genetic models and combines them into single file for each pids.
   - The files generated for each sample will be combined and place in the ANALYSIS_DIR for further analysis.
   - These filtered files were combined to craete a cohort summary.
