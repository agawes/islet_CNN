## bam files: 
#/well/mccarthy/users/agata/allelic_imbalance/rs34584161/*bam
/well/got2d/Matthias/atac_seq/Jim_pipeline/all_fully_mapped_unfiltered_up_to_date_04_Sept_2016/sorted*bam

## vcf:
cp /well/got2d/agata/Basset/HRC_credset/HRC.PPA_0.1_0.5.vcf .
cp /well/got2d/agata/Basset/HRC_credset/HRC.PPA_above0.5.vcf .
cp /well/got2d/agata/Basset/HRC_credset/HRC.PPA_below_0.1.vcf .

### add the few VCF columns
awk '{print $0 "\t.\t.\t.\tGT\t0/0"}' HRC.PPA_above0.5.vcf > HRC.PPA_above0.5.ASEQ.vcf
awk '{print $0 "\t.\t.\t.\tGT\t0/0"}' HRC.PPA_0.1_0.5.vcf > HRC.PPA_0.1_0.5.ASEQ.vcf
awk '{print $0 "\t.\t.\t.\tGT\t0/0"}' HRC.PPA_below_0.1.vcf > HRC.PPA_below_0.1.ASEQ.vcf

### ASEQ:
cd ASEQ_above0.5/
for I in /well/got2d/Matthias/atac_seq/Jim_pipeline/all_fully_mapped_unfiltered_up_to_date_04_Sept_2016/sorted*bam; do
	base=`basename $I .bam`
	/well/got2d/agata/bin/ASEQ/linux64/ASEQ vcf=../HRC.PPA_above0.5.ASEQ.vcf bam=$I threads=1 mbq=1 mrq=1 mdc=1 out=$base
done

mkdir ASEQ_0.1_0.5
cd ASEQ_0.1_0.5/
for I in /well/got2d/Matthias/atac_seq/Jim_pipeline/all_fully_mapped_unfiltered_up_to_date_04_Sept_2016/sorted*bam; do
	base=`basename $I .bam`
	/well/got2d/agata/bin/ASEQ/linux64/ASEQ vcf=../HRC.PPA_0.1_0.5.ASEQ.vcf bam=$I threads=1 mbq=1 mrq=1 mdc=1 out=$base
done

mkdir ASEQ_below_0.1
cd ASEQ_below_0.1/
for I in /well/got2d/Matthias/atac_seq/Jim_pipeline/all_fully_mapped_unfiltered_up_to_date_04_Sept_2016/sorted*bam; do
	base=`basename $I .bam`
	/well/got2d/agata/bin/ASEQ/linux64/ASEQ vcf=../HRC.PPA_below_0.1.ASEQ.vcf bam=$I threads=1 mbq=1 mrq=1 mdc=1 out=$base
done


