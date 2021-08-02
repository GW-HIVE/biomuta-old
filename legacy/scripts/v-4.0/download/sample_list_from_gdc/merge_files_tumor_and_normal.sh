rm merge_files_tumor_and_normal.log
touch merge_files_tumor_and_normal.log

for study in \
        TCGA-KIRP \
        TCGA-KICH \
        TCGA-KIRC \
        TCGA-BRCA \
        TCGA-PRAD \
        TCGA-LUAD \
        TCGA-READ \
        TCGA-ESCA \
        TCGA-THCA \
        TCGA-UCEC \
        TCGA-STAD \
        TCGA-LIHC \
        TCGA-COAD \
        TCGA-HNSC \
        TCGA-LUSC \
        TCGA-BLCA


do

echo "##################################################################################################" >> merge_files_tumor_and_normal.log
echo "python merge_files_tumor_and_normal.py ../files_tcga/${study}  ../files_tcga/merged_gene_filtered/" >> merge_files_tumor_and_normal.log
python merge_files_tumor_and_normal.py     ../files_tcga/${study}    ../files_tcga/merged_gene_filtered/  >> merge_files_tumor_and_normal.log

done
