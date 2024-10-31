## download the cellranger with previous vision 7.2
wget -O cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1730227655&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=oiI4SUIuS0SVkOUjwjg-m6gQvjs5JQt39rvw2POfVCGYZHQ0cNc2iytrOMN3RdoBoDencFKQ3g7oQxJlNS2B2YWdecYaKcj~3nkVhkY-gVNLTz7oIgX0BcEcIMzPUC2nLEbekTTRyGim1bAsSSZE8N8VXz2JB73ufjv1BxpHetxFvQpxCptBewekn7a8QAoM1ESxoTpINtqRzel8IiRHBUPCd4R2RdO~vP7C6CEsqDnYLauIGeU-Ou7snqD1Ik05uEGSpqzcGCUHYNBnPLJio1YLqgQ5LHKVabJOvIDiSyuFZA4XNhCwf79ZR6hkPQMER~Cb5PgAmtjRYTkwVsBHYw__"
tar -xzvf cellranger-7.2.0.tar.gz
export PATH=/home/ec2-user/count/cellranger-7.2.0:$PATH

## download the ref and fastq file from 10X website
wget https://cf.10xgenomics.com/samples/cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_Controller/500_PBMC_3p_LT_Chromium_Controller_fastqs.tar
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"


## run cellranger count
cellranger count --id=500_PBMC --fastqs=500_PBMC_3p_LT_Chromium_Controller_fastqs --sample=500_PBMC_3p_LT_Chromium_Controller --transcriptome=refdata-gex-GRCh38-2024-A

ls -1 outs/

"analysis
cloupe.cloupe
filtered_feature_bc_matrix
filtered_feature_bc_matrix.h5
metrics_summary.csv
molecule_info.h5
possorted_genome_bam.bam
possorted_genome_bam.bam.bai
raw_feature_bc_matrix
raw_feature_bc_matrix.h5
web_summary.html"
