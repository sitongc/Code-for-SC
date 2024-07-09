# Sort sam files
echo "Sort sam file"
gatk SortSam \
--INPUT /groups/ligrp/sitongchen/scrna/sam/${filename}.sam \
--OUTPUT /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted.bam \
--SORT_ORDER coordinate
if [ /groups/ligrp/sitongchen/scrna/sortbam/${filename}_sorted.bam -s ]; then
rm -f /groups/ligrp/sitongchen/scrna/sam/${filename}.sam
fi
