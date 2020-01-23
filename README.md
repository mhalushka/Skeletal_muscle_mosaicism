# Skeletal_muscle_mosaicism
R codes related to Seurat and SimpleSingleCell analyses of skeletal muscle scRNA-seq data

This repository is coding for the manuscript: Proteogenomic single cell analysis of skeletal muscle myocytes by Fomchenko and Verma et al.

Full Jupyter notebooks are also available showing output from these codes.  These are available on request (but large!).

The required files for the code are below:

Mann Cell Data - This file holds the collection of genes from the Mann proteomics paper
https://livejohnshopkins-my.sharepoint.com/:x:/g/personal/mhalush1_jh_edu/ESW25yQwBrlNr6BlW8Ey2zwBnF6mmG_X1mBQgeqvwVboQg?e=uHZlx8

top3K (note, this will download as a xlsx, but save it as a csv --- This is also jsut the list of genes found in rownames(Sobj$gene[["SCT"]]@scale.data)
https://livejohnshopkins-my.sharepoint.com/:x:/r/personal/mhalush1_jh_edu/Documents/Shared%20with%20Everyone/top3kchosengenesSCEpkg.csv?d=w5ae307ccc4db4eca9071054170b8113e&csf=1&e=JhoHqp

009skeljup - This is the main Rdata file you need to load that has the counts
https://livejohnshopkins-my.sharepoint.com/:u:/r/personal/mhalush1_jh_edu/Documents/Shared%20with%20Everyone/009skeljupyter.RData?csf=1&e=erXxXV
