##mlml for tab-seq and wgbs
mkdir $PWD/$3
cd $PWD/$3
bedtools intersect -a ../$2/WGBS_data -b ../$2/TAB_data -wa -wb -sorted >wgbs_inter_tab
cut -f1-5 wgbs_inter_tab >wgbs
cut -f6-10 wgbs_inter_tab >tab
python3 ../src/format_as_mlml.py wgbs wgbs_format
python3 ../src/format_as_mlml.py tab tab_format

##The program will use libgsl.so.19, you have to first find the location and assign the absoluate direction
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/li-lab/Ziwei/Anaconda3/lib  #The path to your libgsl.so.19 folder
export LD_LIBRARY_PATH
../MethPipe/bin/mlml -v -u wgbs_format -h tab_format -o data_mlml

cd $PWD/$3
cut -f5,10 wgbs_inter_tab >wgbs_tab_cov
paste data_mlml wgbs_tab_cov >data_mlml_cov
cut -f1-5,8-9 data_mlml_cov >data_mlml_mc_hmc_cov


##process medip,hmcSeal and mre features
#6w cerebellum
#hmcSeal
sort -k1,1 -k2,2n ../$2/hmC-Seal.bedGraph >hmcSeal_6w_cerebellum1.bedGraph
python3 ../src/extract_feature_from_bp3.py hmcSeal_6w_cerebellum1.bedGraph ../$1/cpg_no_chrM ../$1/mm9_chrom_sizes hmcSeal_window_6w_cerebellum1
cut -f4- hmcSeal_window_6w_cerebellum1 >hmcSeal_window_6w_cerebellum1_cut

#medip
sort -k1,1 -k2,2n ../$2/MeDIP.bedGraph >medip_6w_cerebellum1.bedGraph
python3 ../src/extract_feature_from_bp3.py medip_6w_cerebellum1.bedGraph ../$1/cpg_no_chrM ../$1/mm9_chrom_sizes medip_window_6w_cerebellum1
cut -f4- medip_window_6w_cerebellum1 >medip_window_6w_cerebellum1_cut

#mre
cut -f1-3 ../$2/MRE.bed >mre_WangT_6w-cerebellum-1_cut.filter.bed
sort -k1,1 -k2,2n mre_WangT_6w-cerebellum-1_cut.filter.bed >mre_WangT_6w-cerebellum-1_sort.filter.bed
bedtools genomecov -i mre_WangT_6w-cerebellum-1_sort.filter.bed -g ../$1/mm9_chrom_sizes -bg >mre_WangT_6w-cerebellum-1.bedGraph
python3 ../src/extract_feature_from_bp3.py mre_WangT_6w-cerebellum-1.bedGraph ../$1/cpg_no_chrM ../$1/mm9_chrom_sizes mre_window_6w_cerebellum1
cut -f4- mre_window_6w_cerebellum1 >mre_window_6w_cerebellum1_cut


##combine genomic and methylation features
#coord,dist_cgi,cpg_num, gc_content,medip,mre site, mre,hmcSeal
#add cortex_gcContent_window2_cut (change)
paste ../$1/cpg_no_chrM ../$4/cortex_gcContent_window2_cut ../$4/cpg_dist_cgi ../$4/cpg_num_window_final medip_window_6w_cerebellum1_cut ../$4/data_mreSite_final mre_window_6w_cerebellum1_cut hmcSeal_window_6w_cerebellum1_cut >data_feature
cut -f4- data_feature >data_feature_cut
Rscript ../src/convert_cpgNum_to_cpgDensity.R data_feature_cut data_final 
cut -f1-3 data_feature >data_final_coord
paste data_final_coord data_final >data_final_coord_data
