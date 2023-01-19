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
