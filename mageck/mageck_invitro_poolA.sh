##start interactive node and activate conda mageckenv

#module add miniconda3/v4
#source /home/software/conda/miniconda3/bin/condainit
#conda info --envs   
#source activate /home/charliew/.conda/envs/mageckenv
#module add samtools/1.10

## NOTE - Index numbers are 0-based after gene columsn

## Omitting #4 from run

mageck test -k invitro.poolA.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 1,2,3,5,6,7,8,9,10,11,12,13,14 -c 0 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invitro.poolA_all_v_input

mageck test -k invitro.poolA.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 2,3,5 -c 0 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invitro.poolA_hEGFRv3.CAR_ET.1.10_v_input

mageck test -k invitro.poolA.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 9,10,11 -c 2,3,5 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invitro.poolA_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10

mageck test -k invitro.poolA.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 12,13,14 -c 6,7,8 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invitro.poolA_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2
