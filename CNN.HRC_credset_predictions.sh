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

### control which GPU the job is running on, e.g. to run on GPU #3 :
export CUDA_VISIBLE_DEVICES=2

############## run HRC credible set predictions with best model for each of the iterations #################
fs=7
cd fs$fs\_runs/
mkdir HRC_credset_predictions
for i in {1..100}; do
 printf "Running predictions on: $i\n"
 #nohup basset_predict.py --cudnn islets_cnn.fs$fs.iter$i\_best.th ../../HRC_credset/HRC.PPA_above0.5.fasta HRC_credset_predictions/iter$i.HRC_above_0.5.snp_predict.txt
 #nohup basset_predict.py --cudnn islets_cnn.fs$fs.iter$i\_best.th ../../HRC_credset/HRC.PPA_0.1_0.5.fasta HRC_credset_predictions/iter$i.HRC_0.1_0.5.snp_predict.txt
 nohup basset_predict.py --cudnn islets_cnn.fs$fs.iter$i\_best.th ../../HRC_credset/HRC.PPA_below_0.1.fasta HRC_credset_predictions/iter$i.HRC_below_0.1.snp_predict.txt

done

## get SNP row names:
cd ../
grep ">" ../HRC_credset/HRC.PPA_above0.5.fasta | perl -i -pe 's/>//' > HRC.PPA_above0.5.snp_pred.rownames
for f in fs*/HRC_credset_predictions/*HRC_above_0.5.snp_predict.txt; do
	printf "$f\n";
	paste HRC.PPA_above0.5.snp_pred.rownames $f > tmp
	awk '{printf "\t" $1}' data/samples.txt | awk '{print "SNP_allele" $0}' > $f.ann
	cat tmp >> $f.ann
done
rm tmp

grep ">" ../HRC_credset/HRC.PPA_0.1_0.5.fasta | perl -i -pe 's/>//' > HRC.PPA_0.1_0.5.snp_pred.rownames
for f in fs*/HRC_credset_predictions/*HRC_0.1_0.5.snp_predict.txt; do
	printf "$f\n";
	paste HRC.PPA_0.1_0.5.snp_pred.rownames $f > tmp
	awk '{printf "\t" $1}' data/samples.txt | awk '{print "SNP_allele" $0}' > $f.ann
	cat tmp >> $f.ann
done
rm tmp

grep ">" ../HRC_credset/HRC.PPA_below_0.1.fasta | perl -i -pe 's/>//' > HRC.PPA_below_0.1.snp_pred.rownames
for f in fs*/HRC_credset_predictions/*HRC_below_0.1.snp_predict.txt; do
	printf "$f\n";
	paste HRC.PPA_below_0.1.snp_pred.rownames $f > tmp
	awk '{printf "\t" $1}' data/samples.txt | awk '{print "SNP_allele" $0}' > $f.ann
	cat tmp >> $f.ann
done
rm tmp
