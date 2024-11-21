#!/bin/bash

#awk '{print $1, $4, $5}' hg19withID.bed | sort | uniq -c | awk '$1 > 1'

./liftOver /data/shared/biomuta/generated/datasets/2024_10_22/liftover/hg19withID.bed ucscHg19ToHg38.over.chain /data/shared/biomuta/generated/datasets/2024_10_22/liftover/ucsc_hg38withID.bed /data/shared/biomuta/generated/datasets/2024_10_22/liftover/ucsc_unmapped_withID.bed

./liftOver /data/shared/biomuta/generated/datasets/2024_10_22/liftover/ucsc_unmapped_withID.bed ensembl_GRCh37_to_GRCh38.chain /data/shared/biomuta/generated/datasets/2024_10_22/liftover/ensembl_hg38withID.bed /data/shared/biomuta/generated/datasets/2024_10_22/liftover/ensembl_unmapped_withID.bed