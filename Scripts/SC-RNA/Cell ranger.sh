## download the cellranger with previous vision 7.2
curl -o cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1730065409&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=RwSrYbRwT3gYlKDzAzUChYoqB6Q~EBcySqOivpZOXdySLZbXhl5iKwV8622J~Zt~6NdyvpitzJo~Hfgu54FhsTyjLe~zZfri212S3L92yknLbmo2uTcNWVjzNf5XFGcAEWvYxsQhiFDrgwse85t1DR-NfxPAveMplgsurcJ3rIApXUfy8INp-X19I0X0wnxXodkTVSbNEtXqn3PTI1wBsU6zT59DSfpvvQOkAOXG4-ZZt2Qj8EWA5jYatCs5-mLwbzJ4vX6z2R2emXh~ujrWMdWhtSmM73kEAPin45IOWlYFRE3WHjSgi0x-XQ3XZIX9qWwwDJcjvMLW7JBkIgAazQ__"
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
