#Mutation detection method

##Data cleaning

*fastp was used to cleaning raw fastq data.

		for filename in *1.fastq.gz; do indexNo=`basename ${filename} 1.fastq.gz`; echo ${indexNo}; fastp -i ${indexNo}1.fastq.gz -I ${indexNo}2.fastq.gz -3 -o ${indexNo}paired_1.fastq.gz -O ${indexNo}paired_2.fastq.gz -h ${indexNo}report.html -j ${indexNo}report.json -q 20 -n 10 -t 1 -T 1 -l 20 -w 16; done；
