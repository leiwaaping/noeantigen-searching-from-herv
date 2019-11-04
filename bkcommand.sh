grep -E "_3M_|_B7_|_2N_|_BR_|_CD_|_CG_|_D7_|_EQ_|_F1_|_FP_|_H1_|_HA_|_HF_|_HU_|_IN_|_IP_|_KB_|_KZ_|_MX_|_QW_|_R5_|_RD_|_SW_|_TL_|_V6_|_VA_|_VQ_|_VX_|_YX_|_ZA_|_ZQ_" data/JCI121476.sdt12.txt |sort -k1 |less -S
sed -i -e '/_BR_/d' -e '/_CG_/d' -e '/_AB_/d' data/JCI121476.sdt12.txt
jupyter notebook --ip=172.16.36.14 --NotebookApp.iopub_data_rate_limit=2147483647
awk -v FS="," '{if ($5=="f") print $2,$4,$4+$3*3,$5;else print $2,$4-$3*3,$4,$5}'115_unique_pep_by_IEDBcancer.csv
