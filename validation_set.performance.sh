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



########## precision-recall curves - on validation set
cd /well/got2d/agata/Basset/Full_islet_dataset/data
grep -A 1 -w chr2 learn_islets.chr1_2.fa > validation_set.chr2.fa

mkdir validation_set_predictions
for i in {1..100}; do
 printf "Running predictions on: $i\n"
 nohup basset_predict.py --cudnn islets_cnn.fs7.iter$i\_best.th /well/got2d/agata/Basset/Full_islet_dataset/data/validation_set.chr2.fa validation_set_predictions/iter$i.validation.txt
done


