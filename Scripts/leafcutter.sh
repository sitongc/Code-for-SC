i=$(awk '{print " " $1}' data.txt)
for data in $i
do
aws s3 cp s3://zuchnerlab/RNASeq/Schule/leafcutter/${data}.junc .
echo ${data}.junc >> sample.txt

python leafcutter/clustering/leafcutter_cluster.py -j sample.txt  -m 50 -o RNA_${data} -l 500000 
./leafcutter/scripts/leafcutterMD.R --num_threads 8  RNA_${data}_perind_numers.counts.gz 
mv leafcutter_outlier_pVals.txt RNA_${data}_outlier_pVals.txt
mv leafcutter_outlier_effSize.txt RNA_${data}_outlier_effSize.txt
mv leafcutter_outlier_clusterPvals.txt RNA_${data}_outlier_clusterPvals.txt
aws s3 cp RNA_${data}_outlier_pVals.txt s3://zuchnerlab/RNASeq/Schule/leafcutter/
aws s3 cp RNA_${data}_outlier_effSize.txt s3://zuchnerlab/RNASeq/Schule/leafcutter/
aws s3 cp RNA_${data}_outlier_clusterPvals.txt s3://zuchnerlab/RNASeq/Schule/leafcutter/ 
aws s3 cp RNA_${data}_perind_numers.counts.gz s3://zuchnerlab/RNASeq/Schule/leafcutter/
aws s3 cp RNA_${data}_perind.counts.gz s3://zuchnerlab/RNASeq/Schule/leafcutter/
sed -i '$d' sample.txt
done


python leafcutter/clustering/leafcutter_cluster_regtools.py -j test_juncfiles.txt -m 50 -o TUB0805 -l 500000
leafcutter/scripts/leafcutter_ds.R --num_threads 4 TUB0805_perind_numers.counts.gz groups_file.txt -i 3
./leafcutter/leafviz/gtf2leafcutter.pl -o hg37 ref/Homo_sapiens.GRCh37.75.gtf
leafcutter/leafviz/prepare_results.R -m groups_file.txt TUB0805_perind_numers.counts.gz leafcutter_ds_cluster_significance.txt leafcutter_ds_effect_sizes.txt hg37


