mageck test -k stockCounts.txt --norm-method control --control-sgrna negCont_guides.txt -t KA1.RPMI_5,KA2.RPMI_5,KA3.RPMI_5,\
KA4.MPM_5,KA5.MPM_5,KA6.MPM_5,KA7.RPMI_21,KA8.RPMI_21,KA9.RPMI_21,KA10.MPM_21,KA11.MPM_21,\
KA12.MPM_21,AR2466_BM,AR2466_SP,AR2467_BM,AR2467_SP,AR2468_BM,AR2468_SP,AR2469_BM,AR2469_SP,AR2473_BM,AR2473_SP,AR2474_BM,\
AR2474_SP,AR2475_BM,AR2475_SP,AR2476_BM,AR2476_SP,AR2480_BM,AR2480_SP,AR2481_BM,AR2481_SP,AR2482_BM,AR2482_SP,AR2483_BM,AR2483_SP,\
ahEGFRv3_1to10rep1,ahEGFRv3_1to10rep2,ahEGFRv3_1to10rep3,ahEGFRv3_1to2rep1,ahEGFRv3_1to2rep2,ahEGFRv3_1to2rep3,amCD19_1to10rep1,\
amCD19_1to10rep2,amCD19_1to10rep3,NoCAR_T_rep1,NoCAR_T_rep2,NoCAR_T_rep3,amCD19_1to2rep1,amCD19_1to2rep2,amCD19_1to2rep3,\
Input1_DOI_rep1,Input1_DOI_rep2,Input1_DOI_rep3,Input2_ACT_rep1,Input2_ACT_rep2,Input2_ACT_rep3 \
-c pDNA_CP1704 \
--gene-test-fdr-threshold 0.05 --adjust-method fdr --gene-lfc-method median --normcounts-to-file -n batch

