#!/bin/zsh
set -euo pipefail

echo "==> Comprehensive cleanup to get under 2GB limit"

# Set environment variable to suppress git-filter-branch warnings
export FILTER_BRANCH_SQUELCH_WARNING=1

# Remove ALL large data files from Git history
echo "==> Removing all large data files from Git history"
git filter-branch --force --index-filter \
  'git rm --cached --ignore-unmatch \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/processed/GSE222520_combined_expression_matrix.rds" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/processed/GSE222520_combined_expression_matrix_test.csv" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925398_IWR4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925395_IWR1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925386_IMR2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925393_IWP3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925392_IWP2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925385_IMR1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925396_IWR2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925383_IMP3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925378_NGB1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925379_NGB2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925390_IMR6.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925389_IMR5.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925388_IMR4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925387_IMR3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925384_IMP4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925382_IMP2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925381_IMP1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925380_IMG2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925377_IMG1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925376_IGP4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925375_IGP3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925374_IGP2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925373_IGP1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925372_IGR4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925371_IGR3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925370_IGR2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925369_IGR1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925368_IGM4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925367_IGM3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925366_IGM2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925365_IGM1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925364_IGL4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925363_IGL3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925362_IGL2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925361_IGL1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925360_IGK4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925359_IGK3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925358_IGK2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925357_IGK1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925356_IGJ4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925355_IGJ3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925354_IGJ2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925353_IGJ1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925352_IGI4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925351_IGI3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925350_IGI2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925349_IGI1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925348_IGH4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925347_IGH3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925346_IGH2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925345_IGH1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925344_IGF4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925343_IGF3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925342_IGF2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925341_IGF1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925340_IGE4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925339_IGE3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925338_IGE2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925337_IGE1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925336_IGD4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925335_IGD3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925334_IGD2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925333_IGD1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925332_IGC4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925331_IGC3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925330_IGC2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925329_IGC1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925328_IGB4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925327_IGB3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925326_IGB2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925325_IGB1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925324_IGA4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925323_IGA3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925322_IGA2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925321_IGA1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925320_IG9.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925319_IG8.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925318_IG7.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925317_IG6.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925316_IG5.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925315_IG4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925314_IG3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925313_IG2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925312_IG1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925311_IF4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925310_IF3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925309_IF2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925308_IF1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925307_IE4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925306_IE3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925305_IE2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925304_IE1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925303_ID4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925302_ID3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925301_ID2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925300_ID1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925299_IC4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925298_IC3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925297_IC2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925296_IC1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925295_IB4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925294_IB3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925293_IB2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925292_IB1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925291_IA4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925290_IA3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925289_IA2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925288_IA1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925287_I9.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925286_I8.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925285_I7.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925284_I6.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925283_I5.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925282_I4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925281_I3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925280_I2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925279_I1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925278_H4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925277_H3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925276_H2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925275_H1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925274_G4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925273_G3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925272_G2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925271_G1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925270_F4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925269_F3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925268_F2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925267_F1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925266_E4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925265_E3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925264_E2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925263_E1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925262_D4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925261_D3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925260_D2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925259_D1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925258_C4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925257_C3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925256_C2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925255_C1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925254_B4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925253_B3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925252_B2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925251_B1.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925250_A4.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925249_A3.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925248_A2.tar.gz" \
    "My_Research/Sc_Immunology_Analysis/Sc_Immune_Analysis/data/raw/GSM6925247_A1.tar.gz"' \
  --prune-empty --tag-name-filter cat -- --all

# Clean up the repository
echo "==> Cleaning up repository"
git for-each-ref --format="%(refname)" refs/original/ | xargs -n 1 git update-ref -d
git reflog expire --expire=now --all
git gc --prune=now --aggressive

echo "==> Repository size after cleanup:"
du -sh .git

echo "==> Ready to push (if under 2GB)"
