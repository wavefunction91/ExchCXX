#!/bin/sh


for xc in \
  SlaterExchange \
  VWN3 \
  VWN_RPA \
  VWN\
  PW91_LDA \
  PW91_LDA_MOD \
  PW91_LDA_RPA \
  PZ81 \
  PZ81_MOD \
  EPC17_1 \
  EPC17_2 \
  EPC18_1 \
  EPC18_2 
do
  cp cuda_impl_template.impl cuda_impl_${xc}.cu
  sed -i "s/XCTYPE/LDA/g" cuda_impl_${xc}.cu
  sed -i "s/XCNAME/${xc}/g" cuda_impl_${xc}.cu
done

for xc in \
  B88   \
  LYP   \
  PBE_X \
  RevPBE_X \
  PBE_C \
  B97_D \
  ITYH_X \
  ITYH_X_033 \
  ITYH_X_015 \
  P86_C \
  P86VWN_FT_C \
  PW91_C \
  PBE_SOL_C \
  BMK_C \
  N12_C \
  N12_SX_C \
  SOGGA11_X_C \
  PW91_X \
  MPW91_X \
  OPTX_X \
  RPBE_X \
  SOGGA11_X_X \
  PW86_X \
  WB97_XC \
  WB97X_XC \
  WB97X_V_XC \
  WB97X_D_XC \
  WB97X_D3_XC \
  HJS_PBE_X \
  LCwPBE_wPBEh_X \
  LRCwPBE_HJS_PBE_X \
  LRCwPBEh_HJS_PBE_X \
  WPBEh_X_default0 \
  HSE03_wPBEh_X \
  HSE06_wPBEh_X 
do
  cp cuda_impl_template.impl cuda_impl_${xc}.cu
  sed -i "s/XCTYPE/GGA/g" cuda_impl_${xc}.cu
  sed -i "s/XCNAME/${xc}/g" cuda_impl_${xc}.cu
done



for xc in \
  SCAN_X \
  SCAN_C \
  R2SCAN_X \
  R2SCAN_C \
  FT98_X \
  M062X_X \
  M062X_C \
  PKZB_X \
  PKZB_C \
  TPSS_X \
  RevTPSS_X \
  M06_L_X \
  M06_X \
  M06_HF_X \
  RevM06_L_X \
  M06_SX_X \
  M06_L_C \
  M06_C \
  M06_HF_C \
  RevM06_L_C \
  M06_SX_C \
  M05_2X_C \
  M05_C \
  M08_HX_C \
  M08_SO_C \
  CF22D_C \
  M11_C \
  MN12_L_C \
  MN12_SX_C \
  MN15_C \
  MN15_L_C \
  TPSS_C \
  RevTPSS_C \
  RSCAN_C \
  BC95_C \
  MBEEF_X \
  RSCAN_X \
  BMK_X \
  M08_HX_X \
  M08_SO_X \
  MN12_L_X \
  MN15_L_X \
  MN15_X \
  CF22D_X \
  MN12_SX_X \
  M11_X \
  M05_X \
  M05_2X_X \
  PC07_K \
  PC07OPT_K \
  SCANL_C \
  SCANL_X \
  R2SCANL_C \
  R2SCANL_X 
do
  cp cuda_impl_template.impl cuda_impl_${xc}.cu
  sed -i "s/XCTYPE/MGGA/g" cuda_impl_${xc}.cu
  sed -i "s/XCNAME/${xc}/g" cuda_impl_${xc}.cu
done
