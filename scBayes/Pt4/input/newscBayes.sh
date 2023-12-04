#!/bin/bash

module load gcc/10.2.0 htslib/1.9 scbayes/1.0.0
module load python/3.9.7

VCFNAME=combined
SAMPLE=$1
VCF=${VCFNAME}.norm.vcf
YML1NAME=pre
YML1=${YML1NAME}SubcloneAssign.yml
YML2NAME=post
YML2=${YML2NAME}SubcloneAssign.yml
YML3NAME=lymphocyte
YML3=${YML3NAME}SubcloneAssign.yml

#scGenotype $1.bam <(zcat $1.barcodes.tsv.gz) $VCF >$1.${VCFNAME}.genotype
#~/script/gt2af.py <$1.${VCFNAME}.genotype >$1.${VCFNAME}.scaf.tsv
#scAssign -q 4 $YML1 $1.${VCFNAME}.genotype >$1.${YML1NAME}.assign
#scAssign -q 4 $YML2 $1.${VCFNAME}.genotype >$1.${YML2NAME}.assign
scAssign -q 4 $YML3 $1.${VCFNAME}.genotype >$1.${YML3NAME}.assign
