#!/bin/bash

bwa_bin='/BiO/apps/bwa-0.7.17/bwa'
samtools_bin='/BiO/apps/samtools/samtools'
gatk_bin='/BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar'
gatk3_bin='/BiO/apps/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar'
bcftools_bin='/BiO/apps/bcftools/bcftools'
snpeff_bin='/BiO/apps/snpEff/snpEff.jar'
snpsift_bin='/BiO/apps/snpEff/SnpSift.jar'
refseq='/BiO/data/reference/hg19.fasta'
target_bed='/BiO/data/target/target.bed'
target_interval='BiO/data/target/target.interval_list'
sickle_bin='/BiO/apps/sickle/sickle'
dbsnp_file='/BiO/data/DB/dbSnp151_chr.vcf.gz'
indel_goldstandard='/BiO/data/DB/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz'
clinvar_file='/BiO/data/DB/clinvar_20200728.vcf.gz'
dbnsfp_file='/BiO/data/DB/dbNSFP2.9.3.txt.gz'


if [ -z "$1" ]; then
	echo
	echo "usage: $0 [passward]"
	echo
	exit
else
	if [ "$1" == "wespipe" ];then

		echo "WES pipeline activated."


		for sample in NA12878 NA12891 NA12892;do
			if [ ! -d "${sample}" ];then
				mkdir ${sample}
			fi

#			ln -s /BiO/data/raw_data/${sample}_1.fastq.gz ${sample}/${sample}_1.fastq.gz
#			ln -s /BiO/data/raw_data/${sample}_2.fastq.gz ${sample}/${sample}_2.fastq.gz

#			$sickle_bin pe -t sanger -g -f ${sample}/${sample}_1.fastq.gz -r ${sample}/${sample}_2.fastq.gz -o ${sample}/${sample}-trimmed_1.fastq.gz -p ${sample}/${sample}-trimmed_2.fastq.gz -s ${sample}/${sample}-trimmed_3.fastq.gz

#			$bwa_bin mem -M -R "@RG\tID:BioEdu\tSM:${sample}\tPL:illumina\tLB:WES" -t 6 $refseq ${sample}/${sample}-trimmed_1.fastq.gz ${sample}/${sample}-trimmed_2.fastq.gz | $samtools_bin view -bS -q 20 - | $samtools_bin sort -o ${sample}/${sample}.sorted.bam

			cp /BiO/data/example/${sample}.sorted.bam ${sample}/${sample}.sorted.bam

			java -Xmx4g -jar $gatk_bin MarkDuplicates -I ${sample}/${sample}.sorted.bam -O ${sample}/${sample}.rmdup.bam -M ${sample}/${sample}.rmdup.metrics --REMOVE_DUPLICATES=true

			$samtools_bin index ${sample}/${sample}.rmdup.bam

			java -Xmx4g -jar $gatk_bin BaseRecalibrator -R $refseq -I ${sample}/${sample}.rmdup.bam -L $target_bed --known-sites $dbsnp_file --known-sites $indel_goldstandard -O ${sample}/${sample}_recal_data.table
	
			java -Xmx4g -jar $gatk_bin ApplyBQSR -bqsr ${sample}/${sample}_recal_data.table -I ${sample}/${sample}.rmdup.bam -O ${sample}/${sample}.recal.bam

			$samtools_bin index ${sample}/${sample}.recal.bam

			java -Xmx4g -jar $gatk3_bin -T DepthOfCoverage -R $refseq -I ${sample}/${sample}.recal.bam -o ${sample}/${sample}_target_cov -ct 1 -ct 5 -ct 10 -ct 20 -ct 30 -omitBaseOutput -L $target_bed

			bedtools intersect -wa -a ${sample}/${sample}.recal.bam -b $target_bed > ${sample}/${sample}.target.bam

			$samtools_bin index ${sample}/${sample}.target.bam

			java -Xmx4g -jar $gatk_bin HaplotypeCaller -R $refseq -I ${sample}/${sample}.target.bam -O ${sample}/${sample}.gatk.vcf

			bgzip ${sample}/${sample}.gatk.vcf
	
			tabix -p vcf ${sample}/${sample}.gatk.vcf.gz

			$bcftools_bin mpileup -f $refseq ${sample}/${sample}.target.bam | $bcftools_bin call -mv -Ov -o ${sample}/${sample}.samt.vcf

			bgzip ${sample}/${sample}.samt.vcf
	
			tabix -p vcf ${sample}/${sample}.samt.vcf.gz

			vcf-isec -o -n +2 ${sample}/${sample}.gatk.vcf.gz ${sample}/${sample}.samt.vcf.gz > ${sample}/${sample}.consensus.vcf

			java -Xmx4g -jar $gatk_bin VariantFiltration -R $refseq -O ${sample}/${sample}.consensus.filt.vcf --variant ${sample}/${sample}.consensus.vcf --filter-expression "DP < 10 || FS > 60.0" --filter-name "LOWQUAL"

			cat ${sample}/${sample}.consensus.filt.vcf | awk -F '\t' '($7 != "LOWQUAL") {print}' | bgzip > ${sample}/${sample}.final.vcf.gz

			tabix -p vcf ${sample}/${sample}.final.vcf.gz

		done
	elif [ "$1" == "ann" ];then

		echo "Annotation activated."

		cp /BiO/data/example/WES.vcf.gz .
	
		java -Xmx4g -jar $snpsift_bin filter 'DP>60' WES.vcf.gz | bgzip > WES.vcf.filt.gz

		java -Xmx4g -jar $snpeff_bin -lof -s WES_summary.html hg19 WES.vcf.filt.gz | bgzip > WES.snpeff.vcf.gz

		java -jar $snpsift_bin annotate $dbsnp_file WES.snpeff.vcf.gz | bgzip > WES.snpeff.dbSnp151.vcf.gz

		java -jar $snpsift_bin annotate $clinvar_file WES.snpeff.dbSnp151.vcf.gz | bgzip > WES.snpeff.dbSnp151.clinvar.vcf.gz

		java -jar $snpsift_bin dbnsfp -db $dbnsfp_file WES.snpeff.dbSnp151.clinvar.vcf.gz | bgzip > WES.snpeff.dbSnp151.clinvar.dbnsfp.vcf.gz

	else
		echo 
		echo "( T.T)/ wrong passward!"
		echo
		exit
	fi
fi
