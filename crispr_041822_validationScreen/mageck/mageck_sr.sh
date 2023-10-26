#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3 \
#-c ahEGFRv3_1to10rep1,ahEGFRv3_1to10rep2,ahEGFRv3_1to10rep3 \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n cd19_1to10_v_egfr_1to10

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t amCD19_1to2rep1,amCD19_1to2rep2,amCD19_1to2rep3 \
#-c ahEGFRv3_1to2rep1,ahEGFRv3_1to2rep2,ahEGFRv3_1to2rep3 \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n cd19_1to2_v_egfr_1to2

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3 \
#-c No.CAR.T_rep1,No.CAR.T_rep2,No.CAR.T_rep3 \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n cd19_1to10_v_NOCAR_T

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3 \
#-c Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n cd19_1to10_v_Input1_DOI

#######

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t AR2466_BM,AR2467_BM,AR2468_BM,AR2469_BM \
#-c Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n BM_egfr.15m_v_Input1_DOI

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t AR2466_SP,AR2467_SP,AR2468_SP,AR2469_SP \
#-c Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n SP_egfr.15m_v_Input1_DOI

######

mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
-t AR2473_BM,AR2474_BM,AR2475_BM,AR2476_BM \
-c Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 \
--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n BM_cd19.10m_v_Input1_DOI

mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
-t AR2473_SP,AR2474_SP,AR2475_SP,AR2476_SP \
-c Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 \
--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n SP_cd19.10m_v_Input1_DOI

mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
-t AR2480_BM,AR2481_BM,AR2482_BM,AR2483_BM \
-c Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 \
--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n BM_cd19.15m_v_Input1_DOI

mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
-t AR2480_SP,AR2481_SP,AR2482_SP,AR2483_SP \
-c Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 \
--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n SP_cd19.15m_v_Input1_DOI

#########

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t ahEGFRv3_1to10rep1,ahEGFRv3_1to10rep2,ahEGFRv3_1to10rep3 \
#-c Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n egfr_1to10_v_Input1_DOI

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t ahEGFRv3_1to2rep1,ahEGFRv3_1to2rep2,ahEGFRv3_1to2rep3 \
#-c Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n egfr_1to2_v_Input1_DOI

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3 \
#-c Input2_ACT_rep1,Input2_ACT_rep2,Input2_ACT_rep3 \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n cd19_1to10_v_Input1_ACT

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t amCD19_1to2rep1,amCD19_1to2rep2,amCD19_1to2rep3 \
#-c Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3 \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n cd19_1to2_v_Input1_DOI

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t amCD19_1to2rep1,amCD19_1to2rep2,amCD19_1to2rep3 \
#-c Input2_ACT_rep1,Input2_ACT_rep2,Input2_ACT_rep3 \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n cd19_1to2_v_Input1_ACT

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3 \
#-c amCD19_1to2rep1,amCD19_1to2rep2,amCD19_1to2rep3 \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n cd19_1to10_v_cd19_1to2

###########

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t AR2480_BM,AR2481_BM,AR2482_BM,AR2483_BM \
#-c AR2466_BM,AR2467_BM,AR2468_BM,AR2469_BM \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n BM_cd19.15m_v_egfr.15m

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t AR2473_BM,AR2474_BM,AR2475_BM,AR2476_BM \
#-c AR2466_BM,AR2467_BM,AR2468_BM,AR2469_BM \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n BM_cd19.10m_v_egfr.15m

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t AR2480_SP,AR2481_SP,AR2482_SP,AR2483_SP \
#-c AR2466_SP,AR2467_SP,AR2468_SP,AR2469_SP \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n SP_cd19.15m_v_egfr.15m

#mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt \
#-t AR2473_SP,AR2474_SP,AR2475_SP,AR2476_SP \
#-c AR2466_SP,AR2467_SP,AR2468_SP,AR2469_SP \
#--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n SP_cd19.10m_v_egfr.15m

#AR2466_BM,AR2466_SP
#AR2467_BM,AR2467_SP
#AR2468_BM,AR2468_SP
#AR2469_BM,AR2469_SP
#AR2473_BM,AR2473_SP
#AR2474_BM,AR2474_SP
#AR2475_BM,AR2475_SP
#AR2476_BM,AR2476_SP
#AR2480_BM,AR2480_SP
#AR2481_BM,AR2481_SP
#AR2482_BM,AR2482_SP
#AR2483_BM,AR2483_SP

#ahEGFRv3_1to10rep1,ahEGFRv3_1to10rep2,ahEGFRv3_1to10rep3
#ahEGFRv3_1to2rep1,ahEGFRv3_1to2rep2,ahEGFRv3_1to2rep3
#amCD19_1to10rep1,amCD19_1to10rep2,amCD19_1to10rep3
#No.CAR.T_rep1,No.CAR.T_rep2,No.CAR.T_rep3
#amCD19_1to2rep1,amCD19_1to2rep2,amCD19_1to2rep3
#Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3
#Input2_ACT_rep1,Input2_ACT_rep2,Input2_ACT_rep3

#pDNA_CP1704
