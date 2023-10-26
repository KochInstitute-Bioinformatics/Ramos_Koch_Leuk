##start interactive node and activate conda mageckenv

#module add miniconda3/v4
#source /home/software/conda/miniconda3/bin/condainit
#conda info --envs   
#source activate /home/charliew/.conda/envs/mageckenv
#module add samtools/1.10

## NOTE - Index numbers are 0-based after gene columsn

mageck test -k invitro.poolE.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 1,2,3,4,5,6,7,8,9,10,11,12 -c 0 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invitro.poolE_all_v_input

mageck test -k invitro.poolE.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 2,3 -c 0 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invitro.poolE_hEGFRv3.CAR_ET.1.10_v_input

mageck test -k invitro.poolE.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 7,8,9 -c 2,3 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invitro.poolE_mCD19.CAR_ET.1.10_v_hEGFRv3.CAR_ET.1.10

mageck test -k invitro.poolE.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 10,11,12 -c 4,5,6 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invitro.poolE_mCD19.CAR_ET.1.2_v_hEGFRv3.CAR_ET.1.2
