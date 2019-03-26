#!/bin/bash
#avant de lancer le script : source activate kipoi-DeepSEA__predict
cd ~/link_epi_to_expr/
#CLEAR PREVIOUS TASKS
for DIRECTORY in temp/ temp/intervals/  temp/predictions/ temp/fa_output/ temp/expression/ temp/temppredictions/;do
   if [ ! -d "$DIRECTORY" ]; then
      mkdir "$DIRECTORY";else
      rm "$DIRECTORY"*
   fi
done
#less nbr_snp_per_gene |sort -nr|less
#TO EXTRACT FA SEQ:
for gene in `less nbr_snp_per_gene |sort -nr|sed "s/^ *//"|cut -f2 -d ' '|head -n 200`;do
   for patient in `more ~/link_epi_to_expr/list_commun_patients`;do
      samtools faidx /data1/antoine/hs37d5.fa.gz 1:`echo $gene|awk '{print $1-500}'`-`echo $gene|awk '{print $1+700}'` |bcftools consensus -s $patient /data1/antoine/genome_seq/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz|sed "s/^>\(.*\)/>chr$patient:\1/";
   done>~/link_epi_to_expr/temp/fa_output/out`echo $gene`.fa &
done
wait

more temp/fa_output/out149230716.fa.fai|cut -f2|sort -n|head -n 1

for file in temp/fa_output/out*.fa;do
   samtools faidx $file;
   if [ 445 != `more "$file".fai | wc -l` ] || [ `more "$file".fai|cut -f2|sort -n|head -n 1` -lt 1000 ]
   then
      rm $file $file.fai
      #mkdir "$DIRECTORY";
   fi
done


#TO CREATE INTERVALS FILE:
for gene in `less nbr_snp_per_gene |sort -nr|sed "s/^ *//"|cut -f2 -d ' '|head -n 200`;do
   more ~/link_epi_to_expr/temp/fa_output/out`echo $gene`.fa |grep "^>"|sed "s/^>//"|sed "s/$/\t1\t1000/"|split - -l 9 ~/link_epi_to_expr/temp/intervals/`echo $gene`
done

#TO PREDICT EPI:
for gene in `less nbr_snp_per_gene |sort -nr|sed "s/^ *//"|cut -f2 -d ' '|head -n 200`;do
   if [ -e temp/fa_output/out`echo $gene`.fa ]
   then
      for x in temp/intervals/`echo $gene`* ;do
         kipoi predict DeepSEA/predict   --dataloader_args='{"intervals_file": '$x', "fasta_file": "temp/fa_output/out'$gene'.fa"}' -o 'temp/temppredictions/predict.'`basename $x`'.tsv' &
      done
      wait
      #MERGE INTERVALS
      cat `ls temp/temppredictions/ -1 | egrep 'predict\.'$gene'..\.tsv'|awk '{print "temp/temppredictions/"$0}'`|grep -v 'm'|cat predictions_header.tsv - > temp/predictions/predictions$gene.tsv
      rm temp/temppredictions/*
      #OUTPUT EXPRESSION
      less /data1/antoine/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |grep -e '^T' -e $'\t1\t'$gene$'\t'>temp/expression/$gene.tsv
   fi
done


rm temp/fa_output/* temp/intervals/*
echo "starting R processing"
Rscript correlation.R
echo "finished!"
rm temp/expression/* temp/predictions/*
