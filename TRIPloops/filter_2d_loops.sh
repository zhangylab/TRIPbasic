#arg1: basis loop file, bedpe format, i.e., p<0.05 & dmin>=1500 for a 250bp FitHiChIP called loop file; p<0.05 & dmin>=3000 for a 500bp FitHiChIP called loop file
#arg2: composite annotation bed file, bed format, i.e., sorted unmerged ENCODE peak file, H3K4me1 H3K4me3 H3K27ac CTCF RAD21 PolII
#arg3: qvalue threshold

# bash /mirror/Data/Gengyao_Chen/TRIP_article_figures/2D_loop_criterion/filter_2d_loops.sh /STORE3/TRIP/24-10-15_H3K27ac/results/FitHiChIP_loops/TRIP_H3K27ac_pool_967M/500/FitHiChIP_Peak2ALL_b500_L1000_U2000000/P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr/TRIP_H3K27ac_pool_967M_500.interactions_FitHiC_dmin3000_p0.05.bed /STORE1/DATABASE/Epigenetics/hg19/ChIP-seq/HeLa/composite_sorted_unmerged.bed 0.05

mkdir -p ov_analysis
awk -vOFS="\t" '{print $1,$2,$3,NR}' $1 > ov_analysis/an1.bed
awk -vOFS="\t" '{print $4,$5,$6,NR}' $1 > ov_analysis/an2.bed

awk -vOFS="\t" '{print $1,$2,$3,NR}' $1 | sort-bed -  >  ov_analysis/an1_sorted.bed
awk -vOFS="\t" '{print $4,$5,$6,NR}' $1 | sort-bed - >  ov_analysis/an2_sorted.bed

awk -vOFS="\t" '{print $1,$2,$3,NR}' <(sort-bed $2) > ov_analysis/composite_sorted_unmerged_4col.bed

bedmap --delim "\t" --bp-ovr 1 --echo --echo-map-id-uniq --indicator ov_analysis/an1_sorted.bed ov_analysis/composite_sorted_unmerged_4col.bed > ov_analysis/an1_sorted_ov.bed
bedmap --delim "\t" --bp-ovr 1 --echo --echo-map-id-uniq --indicator ov_analysis/an2_sorted.bed ov_analysis/composite_sorted_unmerged_4col.bed > ov_analysis/an2_sorted_ov.bed

SOURCE=${BASH_SOURCE[0]}
SOURCEdir=$(dirname $SOURCE)
Rscript "$SOURCEdir"/filter_2d_loops.R $1 $3