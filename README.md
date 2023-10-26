# B_ALL_genome_wide_crispr

## Count Data

### Broad PoolQ count data for the primary screen are located in these directories:
- 200406
- 200627
- RH_resequenced_counts_files

### Broad PoolQ count data for the validation screen are located in these directories:

- crispr_041822_validationScreen/fromBroad

## Primary screen 

### data were imported and organized using the script:
- mageck/mageckPrimary_Screen.Rmd

### Pool/Sample assembled count matricies are in mageck directory:
The files are:
- invitro.pool[ABCDEF].txt
- invivo.pool[ABCDEF].txt

### mageck (version 0.5.9.4) test was run on Linux system with these scripts:
- mageck/mageck_invitro_pool[ABCDEF].sh
- mageck/mageck_invivo_pool[ABCDEF].sh

all associated mageck output for these runs are included in the mageck folder

## Validation screen

### Broad PoolQ count data for the validation screen are located in this directory:

- crispr_041822_validationScreen/fromBroad

### data were imported and organized using the script:

- crispr_041822_validationScreen/mageck/mageckFlute_RamosValidationScreen.html

### mageck (version 0.5.9.4) test was run on Linux system with this script:

- crispr_041822_validationScreen/mageck/mageck_sr.sh

mageck output for this run is included in the mageck folder
