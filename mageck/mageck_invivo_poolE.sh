##start interactive node and activate conda mageckenv

#module add miniconda3/v4
#source /home/software/conda/miniconda3/bin/condainit
#conda info --envs   
#source activate /home/charliew/.conda/envs/mageckenv
#module add samtools/1.10

## NOTE - Index numbers are 0-based after gene columsn

mageck test -k invivo.poolE.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27 -c 0 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invivo.poolE_all_v_input

mageck test -k invivo.poolE.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 1,2,3,4,5,6 -c 0 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invivo.poolE_BM_hEGFRv3_15m_v_input

mageck test -k invivo.poolE.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 7,8,9,10,11 -c 1,2,3,4,5,6 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invivo.poolE_BM_mCD19_10m_v_BM_hEGFRv3_15m

mageck test -k invivo.poolE.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 12,13,14,15,16 -c 1,2,3,4,5,6 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invivo.poolE_BM_mCD19_15m_v_BM_hEGFRv3_15m

mageck test -k invivo.poolE.txt --norm-method control --control-sgrna sgrna_controlsBarcodes.txt \
-t 23,24,25,26,27 -c 17,18,19,20,21,22 \
--gene-test-fdr-threshold 0.05 \
--adjust-method fdr --gene-lfc-method median --normcounts-to-file -n invivo.poolE_SP_mCD19_10m_v_SP_hEGFRv3_15m
