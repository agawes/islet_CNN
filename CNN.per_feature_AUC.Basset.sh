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


#################### get per feature AUC's: ####################################
fs=7
mkdir AUC_per_feature
for i in {54..72}; do
 printf "AUC on: $i\n"
 basset_test.lua -cudnn islets_cnn.fs$fs.iter$i\_best.th learn_islets.chr1_2.h5 AUC_per_feature/iter$i
done
