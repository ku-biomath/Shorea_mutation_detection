# Mutation detection method

## Data cleaning

* fastp was used to cleaning raw fastq data.

		for filename in *1.fastq.gz; do indexNo=`basename ${filename} 1.fastq.gz`; echo ${indexNo}; fastp -i ${indexNo}1.fastq.gz -I ${indexNo}2.fastq.gz -3 -o ${indexNo}paired_1.fastq.gz -O ${indexNo}paired_2.fastq.gz -h ${indexNo}report.html -j ${indexNo}report.json -q 20 -n 10 -t 1 -T 1 -l 20 -w 16; done；

## Read mapping

* We used bwa-mem2 for read mapping.
* Maping and marking duplicates were conducted in one line.

		REF=/media/imai2/ssd2/move/1stind/hypo_assembly.fasta
		bwa-mem2 index ${REF};
		for filename in *paired_1.fastq.gz; do indexNo=`basename ${filename} paired_1.fastq.gz`; echo ${indexNo}; bwa-mem2 mem -t 20 -R "@RG\tID:"${indexNo}"\tPL:ILLUMINA\tSM:"${indexNo} ${REF} ${indexNo}paired_1.fastq.gz ${indexNo}paired_2.fastq.gz | samtools fixmate -m - - | samtools sort -| samtools markdup -@8 --reference ${REF} - ${indexNo}.bam ; done;

