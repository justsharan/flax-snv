#!/bin/bash

echo "Removing INDELS, bases around INDELS, multiallelic sites..."
bcftools filter -Ou -g 5:indel,other $1 | bcftools view -Oz -M 2 -v snps -o without_indels.vcf.gz

echo "Hard filtering (DP, QUAL, MQ, SP)..."
bcftools filter -Oz -e "INFO/DP>2.5*AVG(INFO/DP) | QUAL<30 | MQ<30 | SP>3" without_indels.vcf.gz -o quality_filters.vcf.gz

echo "Creating separate files for HIGH and MODERATE variants..."
java -jar $snpEff/SnpSift.jar filter "(ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE')" quality_filters.vcf.gz | gzip > high_and_moderate.vcf.gz
java -jar $snpEff/SnpSift.jar filter "ANN[*].IMPACT = 'HIGH'" high_and_moderate.vcf.gz | gzip > high.vcf.gz
java -jar $snpEff/SnpSift.jar filter "ANN[*].IMPACT = 'MODERATE'" high_and_moderate.vcf.gz | gzip > moderate.vcf.gz

# Requires vcflib tool
echo "Subsample main variants file for analysis..."
bcftools view $1 | vcfrandomsample -r 0.005 | gzip > variants_subset.vcf.gz
