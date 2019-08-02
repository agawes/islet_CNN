##### load torch modules & dependencies
module load torch/20170221-k80-gcc5.4.0
module load python/2.7.10-gcc4.9.3
export BEDTOOLS=/apps/well/bedtools/2.24.0/
export PATH=$BEDTOOLS:$PATH
export BASSETDIR=/well/got2d/agata/Basset/
export PATH=$BASSETDIR/src:$PATH
export PYTHONPATH=$BASSETDIR/src:$PYTHONPATH
export LUA_PATH="$BASSETDIR/src/?.lua;$LUA_PATH"
export PATH=${PATH}:/well/got2d/agata/bin/weblogo/
export PATH=${PATH}:/apps/well/meme/4.11.2_2/bin


### on new compG003
module swap torch/20170221-p100-gcc5.4.0

############### run motifs predictions with best model for each of the iterations #########
### this takes long time... ####
fs=23
mkdir TF_motifs
for i in {1..100}; do
 printf "Running predictions on: $i\n"
 nohup basset_motifs.py -s 1000 -t -o TF_motifs/iter$i islets_cnn.fs$fs.iter$i\_best.th learn_islets.chr1_2.h5 > TF_motifs/iter$i.motifs.log 2>&1 </dev/null 
done

### cleanup all the unnecessary files:
find . -name "*dens.pdf" -exec rm -f {} \;
find . -name "*logo.eps" -exec rm -f {} \;
find . -name "*_logo.fa" -exec rm -f {} \;
find . -name "filter*_heat.pdf" -exec rm -f {} \;
find . -name "filter*_possum.txt" -exec rm -f {} \;
find . -name "model_out.h5" -exec rm -f {} \;
find . -name "sample.h5" -exec rm -f {} \;
find . -name "filters_meme.txt" -exec rm -f {} \;

### now summarize the motifs across the nets:
### report all the motifs with q-value <0.05 and count how many nets found them
for I in fs*/TF_motifs/iter*/tomtom/tomtom.txt; do awk '{if ($6 < 0.05){print $2}}' $I | sort -u ; done  | sort | uniq -c | sort -nr > 1000nets.TF_motifs.summary

### annotate the top motifs with:
perl ../annotate_TF_motifs.pl 1000nets.TF_motifs.summary | perl -i -pe 's/\(//g' | perl -i -pe 's/\)/ /g' | awk '{if ($1 >=50){print $1 "\t" $2 "\t" $5}}' > 1000nets.TF_motifs.top_table

### remove the redundant motifs - by checking similarity with TomTom
perl ../rm_redundant_motifs.pl 1000nets.TF_motifs.top_table > 1000nets.TF_motifs.non_redundant.top_table

### for getting figures for Table 2, find best matching filter, e.g:
cd /well/got2d/agata/Basset/Full_islet_dataset
grep M6114_1.02 fs*runs/TF_motifs/iter*/tomtom/tomtom.txt | sort -k 6 -g | head
# fs19_runs/TF_motifs/iter28/tomtom/tomtom.txt:filter312	M6114_1.02	0	7.88051e-16	5.78429e-13	1.13416e-12	10	TAAGTAAACA	AAAGCAAACA	+

########### check the influence of TF motifs on predictions #########

awk '{print NR "\t" $1}' ../Full_set_IDR_final/data/samples.txt  > targets.txt

fs=15
mkdir TF_motifs_infl
for i in {1..10}; do
 printf "Running predictions on: $i\n"
 nohup basset_motifs_infl.py -m TF_motifs/iter$i/table.txt -s 2000 -b 500 -o TF_motifs_infl/iter$i -t targets.txt --width 7 --height 40 --font 0.5 islets_cnn.fs$fs.iter$i\_best.th learn_islets.chr1_2.h5 2>&1 </dev/null
done


