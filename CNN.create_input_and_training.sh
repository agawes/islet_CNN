#!/bin/bash

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

#### create input files #####

preprocess_features.py -y -m 200 -s 1000 -n -o learn_islets -c ../../data/genomes/human.hg19.genome samples.txt

bedtools getfasta -fi ../../data/genomes/hg19.fa -bed learn_islets.bed -s -fo learn_islets.fa


####  let's hold out chr1 and chr2 - for test and validation sets ##################
#   43029 chr1
#   40506 chr2

CHR1=`awk '{print $1}' learn_islets.bed  | sort | uniq -c | grep -w chr1 | awk '{printf $1}'`
CHR2=`awk '{print $1}' learn_islets.bed  | sort | uniq -c | grep -w chr2 | awk '{printf $1}'`

# BED file
grep -vw chr1 learn_islets.bed | grep -vw chr2 > chr1_2_last.bed
grep -w chr2 learn_islets.bed >> chr1_2_last.bed
grep -w chr1 learn_islets.bed >> chr1_2_last.bed

### act file
grep -vw chr1 learn_islets_act.txt | grep -vw chr2 > chr1_2_last.act.txt
grep -w chr2 learn_islets_act.txt >> chr1_2_last.act.txt
grep -w chr1 learn_islets_act.txt >> chr1_2_last.act.txt

bedtools getfasta -fi ../../data/genomes/hg19.fa -bed chr1_2_last.bed -s -fo learn_islets.chr1_2.fa
seq_hdf5.permute.py -c -v $CHR2 -t $CHR1 learn_islets.chr1_2.fa chr1_2_last.act.txt learn_islets.chr1_2.h5

# 421738 training sequences
# 43029 test sequences
# 40506 validation sequences

#### parameters files for training are in params/ --> 10 different files differing in filter sizes and number of filters

## now create a folder to run 100 iterations of each param and run: nohup perl run_100_CNNs.pl &

mkdir fs7_runs
cd fs7_runs/
ln -s ../data/learn_islets.chr1_2.h5 .
cp ../../100_iterations_fs15/run_100_CNNs.pl .
## make changes in the file
nohup perl run_100_CNNs.pl &


#### cleanup check files:
rm *check.th

