import os
import sys


PATH_TO_GenomeAnalysisTK = ""
PATH_TO_impute2 = ""
GENOME="Data/genome/grape12Xv2.fa"


def impute_ancient(BAM, GENE_NAME):
    os.system("java -jar %s -T UnifiedGenotyper -R %s --genotyping_mode GENOTYPE_GIVEN_ALLELES --output_mode EMIT_ALL_SITES -I %s -o %s -nct 20 --alleles Data/Modern.vcf.gz --dbsnp Data/Modern.vcf.gz -L 2.list" % (PATH_TO_GenomeAnalysisTK ,GENOME, BAM, 'output/tmp/'+GENE_NAME+'.vcf'))
    os.system("python convertVCFToImpute2Input.py -i %s -p 10 -o %s" % ('output/tmp/'+GENE_NAME+'.vcf','output/tmp/'+GENE_NAME+'.gens'))
    os.system("%s -g  %s -m Data/mybA1_genetic_map_SNPs.txt -phase -int 14350000 14353000 -Ne 20000 -buffer 2000 -k 400 --k_hap 2000 -h Data/Modern.haps.gz -l Data/Modern.legend.gz -pgs_prob -prob_g -o %s" % (PATH_TO_impute2, 'output/tmp/'+i+'.gens', 'output/'+ GENE_NAME))

