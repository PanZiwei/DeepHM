PWD=/projects/li-lab/Ziwei/DeepHM/forCompress
##process genomic features like gc percent, cpg num, dist to cpg island
mkdir $PWD/$2
cd $PWD/$2 ##This step will direct to the $pwd/$2 and use it as the working directory later
#gc percent
python3 ../src/extract_feature_from_bp3.py ../$1/gc5Base.sort.bedGraph ../$1/cpg_no_chrM ../$1/mm9_chrom_sizes cortex_gcContent_window2
#remove first 5 column (change)
cut -f5 cortex_gcContent_window2 >cortex_gcContent_window2_cut

#cpg num at each window
python3 ../src/generate_cpg_window.py ../$1/cpg_no_chrM cpg_window
bedtools coverage -a cpg_window -b ../$1/cpg_no_chrM  -counts >cpg_num_window
cut -f4 cpg_num_window >cpg_num_window_cut
Rscript ../src/format_cpgNum_data.R cpg_num_window_cut cpg_num_window_final

#dist to cpg island
bedtools intersect -a ../$1/cpg_no_chrM -b ../$1/cpgIsland.bed -wao -sorted >cpg_inter_cgi
cut -f1-3,7 cpg_inter_cgi >cpg_inter_cgi_cut
python3 ../src/cal_dist_to_cgi.py cpg_inter_cgi_cut cpg_dist_cgi

#genome encoding (add)
python3 -u ../src/create_genome_encoding2.py ../$1/mm9.fa mm9_genome_encoding

#mre site (add)
cat ../$1/ACGT_sites.bed ../$1/GCGC_sites.bed  ../$1/CCGC_sites.bed ../$1/CCGG_sites.bed >four_site.bed
sort -k1,1 -k2,2n four_site.bed >mre_four_site_sort.bed
bedtools intersect -a ../$1/cpg_no_chrM -b   mre_four_site_sort.bed -u -sorted >data_inter_mreSite
sed "s/$/\t1/" data_inter_mreSite >data_inter_mreSite_add1
bedtools intersect -a ../$1/cpg_no_chrM -b mre_four_site_sort.bed -v -sorted >data_notInter_mreSite
sed "s/$/\t0/" data_notInter_mreSite >data_notInter_mreSite_add0
cat data_inter_mreSite_add1 data_notInter_mreSite_add0 >data_mreSite
sort -k1,1 -k2,2n data_mreSite >data_mreSite_sort
cut -f4 data_mreSite_sort>data_mreSite_final
