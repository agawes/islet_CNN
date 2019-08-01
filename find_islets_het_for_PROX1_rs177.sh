cd /well/mccarthy/users/agata/allelic_imbalance/PROX1/
for I in /well/mccarthy/production/atac-seq/data/human_islets/full_merged_data/bams/HP1663*bam; do
        base=`basename $I .bam`
        /well/got2d/agata/bin/ASEQ/linux64/ASEQ vcf=../PROX1.vcf bam=$I threads=1 mbq=1 mrq=1 mdc=1 out=$base
done

#### 17/06/2019 - find all islets heterozygous for the rs177 variant:
cd /well/mccarthy/production/genotyping/data/20.imputed/islets/latest/HRC
tabix chr1.dose.vcf.gz "1:214150445-214150445" | perl -i -pe 's/\t/\n/g' | awk '{print NR "\t" $0}' | grep "1|0"
[W::hts_idx_load2] The index file is older than the data file: chr1.dose.vcf.gz.tbi
53	1|0:1.000:0.000,1.000,0.000
56	1|0:1.000:0.000,1.000,0.000
94	1|0:1.000:0.000,1.000,0.000
97	1|0:1.000:0.000,1.000,0.000
114	1|0:1.000:0.000,1.000,0.000
115	1|0:1.000:0.000,1.000,0.000
135	1|0:0.994:0.006,0.994,0.000
220	1|0:1.000:0.000,1.000,0.000
273	1|0:1.000:0.000,1.000,0.000
347	1|0:0.994:0.006,0.994,0.000
352	1|0:1.000:0.000,1.000,0.000

zcat chr1.dose.vcf.gz | head -30 | grep CHROM | cut -f 53,56,94,97,114,115,135,220,273,347,352
ISL-5	H546	H357	H407	ISL-98	ISL-99	H671	R177	R191	ISL-119	ISL-130

tabix chr1.dose.vcf.gz "1:214150445-214150445" | perl -i -pe 's/\t/\n/g' | awk '{print NR "\t" $0}' | grep "0|1"
[W::hts_idx_load2] The index file is older than the data file: chr1.dose.vcf.gz.tbi
128	0|1:1.000:0.000,1.000,0.000
141	0|1:1.000:0.000,1.000,0.000
152	0|1:1.000:0.000,1.000,0.000
153	0|1:1.000:0.000,1.000,0.000
193	0|1:0.992:0.008,0.992,0.000
201	0|1:1.000:0.000,1.000,0.000
210	0|1:0.998:0.002,0.998,0.000
318	0|1:1.000:0.000,0.999,0.001

zcat chr1.dose.vcf.gz | head -30 | grep CHROM | cut -f 128,141,152,153,193,201,210,318
ISL-86	ISL-183	R082	R094	R140	R106	ISL-167	ISL-231

### in islet base - ATAC-seq data available for:
R191 - not in Thurner 2018
ISL-183 - in Thurner 2018 (HP15-35)
ISL-167 - not in Thurner 2018 (HP15-05)
ISL-231 - not in Thurner 2018 (HP16-63)
