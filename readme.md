# Mutation detection method

## Data cleaning

* fastp was used to cleaning raw fastq data.

	for filename in *1.fastq.gz; do indexNo=`basename ${filename} 1.fastq.gz`; echo ${indexNo}; fastp -i ${indexNo}1.fastq.gz -I ${indexNo}2.fastq.gz -3 -o ${indexNo}paired_1.fastq.gz -O ${indexNo}paired_2.fastq.gz -h ${indexNo}report.html -j ${indexNo}report.json -q 20 -n 10 -t 1 -T 1 -l 20 -w 16; done；

## Read mapping

* We used bwa-mem2 for read mapping.
* Maping and marking duplicates were conducted in one line.

		REF=Your reference genome
		bwa-mem2 index ${REF};
		for filename in *paired_1.fastq.gz; do indexNo=`basename ${filename} paired_1.fastq.gz`; echo ${indexNo}; bwa-mem2 mem -t 20 -R "@RG\tID:"${indexNo}"\tPL:ILLUMINA\tSM:"${indexNo} ${REF} ${indexNo}paired_1.fastq.gz ${indexNo}paired_2.fastq.gz | samtools fixmate -m - - | samtools sort -| samtools markdup -@8 --reference ${REF} - ${indexNo}.bam ; done;

## SNP calling (samtools)

* We used two SNP caller to call accurate somatic mutation. We extracted shared SNPs between two SNP caller
* Indels were normalized after SNP calling.

		ls *.bam | xargs -P20 -I{} samtools index {};
		ls *.bam > bamlist;
		bcftools mpileup -Ou -q 40 -Q 40 --threads 8 -f ${REF} -b bamlist | bcftools call -vmO z --threads 8 -o shorea1stData_samtools_q40Q40.vcf.gz&
		bcftools mpileup -Ou -q 30 -Q 30 --threads 8 -f ${REF} -b bamlist | bcftools call -vmO z --threads 4 -o shorea1stData_samtools_q30Q30.vcf.gz&
		bcftools mpileup -Ou -q 20 -Q 20 --threads 8 -f ${REF} -b bamlist | bcftools call -vmO z --threads 4 -o shorea1stData_samtools_q20Q20.vcf.gz;
		bcftools norm -d all -f ${REF} -O z -o norm_shorea2ndData_samtools_q40Q40.vcf.gz shorea1stData_samtools_q40Q40.vcf.gz&
		bcftools norm -d all -f ${REF} -O z -o norm_shorea2ndData_samtools_q30Q30.vcf.gz shorea1stData_samtools_q30Q30.vcf.gz&
		bcftools norm -d all -f ${REF} -O z -o norm_shorea2ndData_samtools_q20Q20.vcf.gz shorea1stData_samtools_q20Q20.vcf.gz;

* We divided SNPs and indels in the vcf files.

		vcftools --gzvcf norm_shorea2ndData_samtools_q40Q40.vcf.gz --max-missing 1 --max-alleles 2 --remove-indels --out SNP_norm_shorea2ndData_samtools_q40Q40 --recode&
		vcftools --gzvcf norm_shorea2ndData_samtools_q30Q30.vcf.gz --max-missing 1 --max-alleles 2 --remove-indels --out SNP_norm_shorea2ndData_samtools_q30Q30 --recode&
		vcftools --gzvcf norm_shorea2ndData_samtools_q20Q20.vcf.gz --max-missing 1 --max-alleles 2 --remove-indels --out SNP_norm_shorea2ndData_samtools_q20Q20 --recode;
		vcftools --gzvcf norm_shorea2ndData_samtools_q40Q40.vcf.gz --max-missing 1 --keep-only-indels --out indel_norm_shorea2ndData_samtools_q40Q40 --recode&
		vcftools --gzvcf norm_shorea2ndData_samtools_q30Q30.vcf.gz --max-missing 1 --keep-only-indels --out indel_norm_shorea2ndData_samtools_q30Q30 --recode&
		vcftools --gzvcf norm_shorea2ndData_samtools_q20Q20.vcf.gz --max-missing 1 --keep-only-indels --out indel_norm_shorea2ndData_samtools_q20Q20 --recode;

## SNP calling (gatk)
* We also used gatk for SNP calling.

		gatk CreateSequenceDictionary -R ${REF} -O hypo_assembly.dict;
		seqkit seq -n ${REF} >intervals.list;
		samtools faidx ${REF};
		for fpath in *.bam; do fname=`basename ${fpath} .bam` ; echo ${fname}; done | sort |xargs -P4 -I{} gatk HaplotypeCaller -R ${REF} --emit-ref-confidence GVCF -I {}.bam -O {}.g.vcf;
		gvcf_files="";
		for gvcf_file in *.g.vcf;do gvcf_files=${gvcf_files}"-V ${gvcf_file} ";done;
		gatk GenomicsDBImport -R ${REF} ${gvcf_files} -L intervals.list --genomicsdb-workspace-path gvcfs_db;
		gatk GenotypeGVCFs -R ${REF} -V gendb://gvcfs_db -O merged.vcf;
		
* We divided and filtered SNPs and indels in the vcf files. 

		gatk SelectVariants -R ${REF} -V merged.vcf --select-type-to-include SNP -O merged_snps.vcf;
		gatk VariantFiltration -R ${REF} -V merged_snps.vcf -O 1st_merged_snps_filtered.vcf \
                       -filter "QD < 2.0" --filter-name "QD2"       \
                       -filter "QUAL < 30.0" --filter-name "QUAL30" \
                       -filter "SOR > 4.0" --filter-name "SOR4"     \
                       -filter "FS > 60.0" --filter-name "FS60"     \
                       -filter "MQ < 40.0" --filter-name "MQ40"     \
                       -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                       -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8";
		
		gatk SelectVariants -R ${REF} -V merged.vcf --select-type-to-include INDEL -O merged_indels.vcf;
		gatk VariantFiltration -R ${REF} -V merged_indels.vcf -O 1st_merged_indels_filtered.vcf \
                       -filter "QD < 2.0" --filter-name "QD2"       \
                       -filter "QUAL < 30.0" --filter-name "QUAL30" \
                       -filter "FS > 200.0" --filter-name "FS200"   \
                       -filter "SOR > 10.0" -filter-name "SOR10"    \
                       -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20";
		       
## Extracting shared SNP between two SNP caller and replicates.
* We removed fixed sites from vcf the file using Tassel5.
* After removing fixed sites, we extracted shared SNPs from four vcf files.

		python compare_vcf.py your_gatk_file1.vcf your_gatk_file2.vcf your_samtools_file1.vcf your_samtools_file2.vcf
