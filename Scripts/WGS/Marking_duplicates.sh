date
echo "Marking duplicates"
gatk MarkDuplicates \
--INPUT /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted.bam \
--METRICS_FILE /groups/ligrp/sitongchen/scrna/sortbam/${filename}_dupMetrics \
--OUTPUT /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted_dupMarked.bam
if [ /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted_dupMarked.bam -s ]; then
rm -f /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted.bam
fi
