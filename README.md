# WDL pipeline for PCA analysis of AGD data 

## Questions & concerns 

- Docker file to use?:
    - In the VUMC examples scripts just use this docker for plink2, feel like we could use the same: hkim298/plink_1.9_2.0:20230116_20230707
    - Otherwise, BROAD docker for plink2/regenie: https://github.com/broadinstitute/warp-tools/blob/develop/3rd-party-tools/plink-regenie/Dockerfile
- Getting R to run in a WDL script - should we separate out the R script as a different WDL?: 
    - How to get R to run https://support.terra.bio/hc/en-us/articles/4404690637851--2021-July-Demo-How-to-incorporate-your-R-code-into-a-WDL 
         - Heredoc approach: less readable, some R code might not be recognized 
         - Rscript approach more readable, make sure the file is localized, just like for the input files 
            - Rscript file is provided as an input to the WDL 
                    - note that when passing arguments to the Rscript, they'll all be converted to strings, so may need to remake integer values 
    - feel like we need a different docker image for this? 
        - make sure we have the correct R libraries in the docker image? 
    

## AGD data: 35k

- Starting n: 35,024 samples
- Appears to be in HG38 according to this, see 1.2 Table agd36k_bed_all: https://sites.google.com/view/biovu/biovu-data/data-snapshot?pli=1&authuser=3#h.2p10szbhafpa

## Underlying code and resources: 

- GBMI pipeline 
    - https://docs.google.com/document/d/1DtEcfuQ0UqSDkjGd7cy-1hy8XGGAsFjHRcoEw4VGnBo/edit
    - https://github.com/globalbiobankmeta/pca_projection


## Inputs 

- GRCh 38 precomputed PCA loadings and allele frequencies 
- AGD data in plink2 format, per chromosome, in an array, as is provided to the WDL script https://app.terra.bio/#workspaces/VICTR_Cloud_Migration_Pilot/VICTR_BioVU_pilot_Davis_lab/workflows/VICTR_Cloud_Migration_Pilot/VUMCPrepareGenotypeData for example 

- Rscript for plotting PCA projection
- SD/BioVU table with FID, IID (GRID) and EHR race for plotting 


## WDL tasks

- Prepare AGD data: This should basically just be this WDL: https://app.terra.bio/#workspaces/VICTR_Cloud_Migration_Pilot/VICTR_BioVU_pilot_Davis_lab/workflows/VICTR_Cloud_Migration_Pilot/VUMCPrepareGenotypeData, but with the change such that plink2_filter_option is --extract variants.extract 
    - make variant list 
    - extract using plink 2 {scatter per chromosome}
    - merge using plink 2 {can use biostats merge script (I think)}

- Run projection 
    - https://github.com/globalbiobankmeta/pca_projection/blob/master/project_pc.sh
        - maybe just localize this to the WDL script as just one plink2 command 

- Plot projection 
    - since this is a R script (a little different, might require a different docker image), maybe this should be 
