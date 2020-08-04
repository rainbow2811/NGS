#!/bin/bash
if [ -z "$1" ]; then
	echo
	echo "usage: $0 [annotated vcf]"
	echo
	exit
fi

inputvcf=$1

snpsift_bin='/BiO/apps/snpEff/SnpSift.jar'

java -jar $snpsift_bin extractFields -s ',' -e '.' $inputvcf \
CHROM POS ID REF ALT \
GEN[0].GT GEN[1].GT GEN[2].GT GEN[0].AD GEN[1].AD GEN[2].AD \
ANN[*].GENE ANN[*].FEATUREID ANN[*].FEATURE ANN[*].BIOTYPE ANN[*].EFFECT ANN[*].IMPACT ANN[*].HGVS_C ANN[*].HGVS_P ANN[*].ERRORS \
CLNSIG dbNSFP_SIFT_pred
