#!/usr/bin/bash

#./check_correct_TAD.sh

all_ds=(
"ENCSR079VIJ_G401"
"ENCSR312KHQ_SK-MEL-5"
"ENCSR346DCU_LNCaP"
"ENCSR401TBQ_Caki2"
"ENCSR444WCZ_A549"
"ENCSR489OCU_NCI-H460"
"ENCSR549MGQ_T47D"
"ENCSR862OGI_RPMI-7951"
"GSE105194_cerebellum"
"GSE105194_spinal_cord"
"GSE105318_DLD1"
"GSE105381_HepG2"
"GSE118514_RWPE1"
"GSE58752_liver"
"GSE73782_PC3"
"GSE75070_MCF-7_shNS"
"K562"
"Panc1_rep12"
)

all_chrs=( {1..22} )

for ds in ${all_ds[@]};do
    for chr in ${all_chrs[@]};do
        diff ${ds}_40kb/FINAL_DOMAINS/${ds}_chr${chr}_YL_40kb_final_domains.txt ${ds}_40kb/FINAL_DOMAINS_CORRECTED/${ds}_chr${chr}_YL_40kb_final_domains.txt 
    done
done
