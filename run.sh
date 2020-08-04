sample=$1


#3. copy the sorted bam file
cp /BiO/data/example/${sample}.sorted.bam .


#4.Remove duplicates
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar MarkDuplicates -I  ${sample}.sorted.bam -O ${sample}.rmdup.bam -M ${sample}.rmdup.metrics --REMOVE_DUPLICATES=true


#5. Producing index file of bam removed duplicates
/BiO/apps/samtools/samtools index ${sample}.rmdup.bam


#6.BQSR
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar BaseRecalibrator -R /BiO/data/reference/hg19.fasta -I ${sample}.rmdup.bam -L /BiO/data/target/target.bed --known-sites /BiO/data/DB/dbSnp151_chr.vcf.gz --known-sites /BiO/data/DB/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz -O ${sample}_recal_data.table

#7.Apply BQSR
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar ApplyBQSR -bqsr ${sample}_recal_data.table -I ${sample}.rmdup.bam -O ${sample}.recal.bam

#8.On-target coverage
java -Xmx4g -jar /BiO/apps/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T DepthOfCoverage -R /BiO/data/reference/hg19.fasta -I ${sample}.recal.bam -o ${sample}_target_cov -ct 1 -ct 5 -ct 10 -ct 20 -ct 30 -omitBaseOutput -L /BiO/data/target/target.bed


#9.Remove off-target reads
bedtools intersect -wa -a ${sample}.recal.bam -b /BiO/data/target/target.bed > ${sample}.target.bam

/BiO/apps/samtools/samtools index ${sample}.target.bam

# 10.Variant calling - GATK HaplotypeCaller
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar HaplotypeCaller -R /BiO/data/reference/hg19.fasta -I ${sample}.target.bam -O ${sample}.gatk.vcf

bgzip ${sample}.gatk.vcf

tabix -p vcf ${sample}.gatk.vcf.gz

# 11. Variant calling - samtools (bcftools)
/BiO/apps/bcftools/bcftools mpileup -f /BiO/data/reference/hg19.fasta ${sample}.target.bam  | /BiO/apps/bcftools/bcftools call -mv -Ov -o ${sample}.samt.vcf

bgzip ${sample}.samt.vcf
tabix -p vcf ${sample}.samt.vcf.gz

#12.consensus VCF
vcf-isec -o -n +2 ${sample}.gatk.vcf.gz ${sample}.samt.vcf.gz > ${sample}.consensus.vcf

# Filtration
java -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar VariantFiltration -R /BiO/data/reference/hg19.fasta -O ${sample}.consensus.filt.vcf --variant ${sample}.consensus.vcf --filter-expression 'DP < 10 || FS > 60.0' --filter-name 'LOWQUAL'

# Remove 'LOWQUAL' in filter column
cat ${sample}.consensus.filt.vcf | awk -F '\t' '($7!="LOWQUAL"){print}' | bgzip > ${sample}.final.vcf.gz

tabix -p vcf ${sample}.final.vcf.gz

