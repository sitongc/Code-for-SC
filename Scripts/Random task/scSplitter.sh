#scSplitter is a preprocessing tool designed to convert scRNA seq results into a format suitable for mutation or SNV calling at the single cell level for the purpose of lineage tracing.

STAR --runMode genomeGenerate --genomeDir /groups/ligrp/sitongchen/scrna/dpi35  --genomeFastaFiles /groups/ligrp/sitongchen/scrna/22/mm10.fa --sjdbGTFfile /groups/ligrp/sitongchen/scrna/22/mm10.refGene.gtf -- sjdbOverhang 199

python3 scSplitter.py --ig True --f /groups/ligrp/sitongchen/scrna/dpi35/gel_barcode1_list.txt --i /groups/ligrp/sitongchen/scrna/dpi35/fastqs --r /groups/ligrp/sitongchen/scrna/dpi35/result --ind /groups/ligrp/sitongchen/scrna/dpi35/sta
STAR --runMode genomeGenerate --genomeDir /groups/ligrp/sitongchen/scrna/dpi35  --genomeFastaFiles /groups/ligrp/sitongchen/scrna/22/mm10.fa --sjdbGTFfile /groups/ligrp/sitongchen/scrna/22/mm10.refGene.gtf -- sjdbOverhang 199
