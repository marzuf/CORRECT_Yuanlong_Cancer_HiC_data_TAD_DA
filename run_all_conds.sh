#!/usr/bin/bash

 # LIVER
 ./run_pipeline.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc   # 04.07.2019 - 06.07.2019
 ./run_pipeline.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1  #  04.07.2019 - 06.07.2019

 ./run_pipeline.sh GSE58752_liver_40kb TCGAlihc_norm_lihc #  04.07.2019 - 06.07.2019
 ./run_pipeline.sh GSE58752_liver_40kb TCGAlihc_wt_mutCTNNB1  #04.07.2019 - 06.07.2019


#  BREAST

 ./run_pipeline.sh ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas   # 04.07.2019 - 06.07.2019
 ./run_pipeline.sh GSE75070_MCF-7_shNS_40kb TCGAbrca_lum_bas  # 04.07.2019 - 06.07.2019


#  KIDNEY
 ./run_pipeline.sh ENCSR079VIJ_G401_40kb TCGAkich_norm_kich  # 04.07.2019 - 06.07.2019
 ./run_pipeline.sh ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich # 04.07.2019 - 06.07.2019


#  SKIN
 ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf # 04.07.2019 - 06.07.2019
 ./run_pipeline.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf  # 04.07.2019 - 06.07.2019

 ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF # 04.07.2019 - 06.07.2019
 ./run_pipeline.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF # 04.07.2019 - 06.07.2019

 ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1	# 04.07.2019 - 06.07.2019
 ./run_pipeline.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1	#  04.07.2019 - 06.07.2019


#  LUNG
 ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_norm_luad   # 05.07.2019 - 07.07.2019
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad   # 05.07.2019 - 07.07.2019

 ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR # 05.07.2019 - 07.07.2019
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR  #  05.07.2019 - 07.07.2019

 ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker # 05.07.2019 - 07.07.2019
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker  #  05.07.2019 - 07.07.2019

 ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS # 05.07.2019 - 07.07.2019
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS #  05.07.2019 - 07.07.2019

 ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc  # 05.07.2019 - 07.07.2019
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc #  05.07.2019 - 07.07.2019

	
#  PANCREAS
 ./run_pipeline.sh Panc1_rep12_40kb TCGApaad_wt_mutKRAS # 04.07.2019 - 07.07.2019
	

#  PROSTATE	
 ./run_pipeline.sh ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad    # 04.07.2019 - 07.07.2019
 ./run_pipeline.sh GSE73782_PC3_40kb TCGAprad_norm_prad         # 04.07.2019 - 07.07.2019
 ./run_pipeline.sh GSE118514_RWPE1_40kb TCGAprad_norm_prad       #  04.07.2019 - 07.07.2019

#  GBM
 ./run_pipeline.sh GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal # 05.07.2019 - 07.07.2019
 ./run_pipeline.sh GSE105194_cerebellum_40kb TCGAgbm_classical_neural # 05.07.2019 - 07.07.2019
 ./run_pipeline.sh GSE105194_cerebellum_40kb TCGAgbm_classical_proneural # 05.07.2019 - 07.07.2019
 ./run_pipeline.sh GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc # 05.07.2019 - 07.07.2019

 ./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal # 05.07.2019  - 07.07.2019
 ./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_neural # 05.07.2019  - 07.07.2019
 ./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural # 05.07.2019 - 07.07.2019
 ./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc # 05.07.2019  - 07.07.2019


 # COLORECTAL
 ./run_pipeline.sh GSE105318_DLD1_40kb TCGAcoad_msi_mss # 04.07.2019 - 07.07.2019

 # LYMPHOBLAST
 ./run_pipeline.sh K562_40kb TCGAlaml_wt_mutFLT3 # 04.07.2019 - 07.07.2019

