#!/bin/bash
start=`date`
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
for gene in 88199922;do
   for position in `seq -250000 100 249999`;do
      for patient in `more ~/link_epi_to_expr/list_commun_patients`;do
         samtools faidx /data1/antoine/hs37d5.fa.gz 5:`echo $gene $position|awk '{print $1+$2}'`-`echo $gene $position|awk '{print $1+$2+1300}'` |bcftools consensus -s $patient /data1/antoine/genome_seq/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz|sed "s/^>\(.*\)/>chr$patient:$position:\1/";
      done>~/link_epi_to_expr/temp/fa_output/out`echo $gene'_'$position`.fa &
      if [ 15 -lt `jobs -p | grep "^[0-9]"| wc -l` ]
      then
         wait
      fi
   done
   wait
   bcf_end=`date`
   for file in temp/fa_output/out*.fa;do
      samtools faidx $file;
      if [ 445 != `more "$file".fai | wc -l` ] || [ `more "$file".fai|cut -f2|sort -n|head -n 1` -lt 1000 ] || [ `more $file|awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'|sed '1d'|grep -v "^>"|uniq|wc -l` = 1 ]
      then
         rm $file $file.fai
      else
         more $file |grep "^>"|sed "s/^>//"|sed "s/$/\t1\t1000/"|split - -l 9 ~/link_epi_to_expr/temp/intervals/`basename $file`
         for x in temp/intervals/`basename $file`* ;do
            kipoi predict DeepSEA/predict   --dataloader_args='{"intervals_file": '$x', "fasta_file": "'$file'"}' -o 'temp/temppredictions/predict.'`basename $x`'.tsv' &
            if [ 15 -lt `jobs -p | grep "^[0-9]"| wc -l` ]
            then
               wait
            fi
         done
         wait
         rm temp/intervals/*
         #MERGE INTERVALS
         file_base=`basename $file`
         cat `ls temp/temppredictions/ -1 | egrep 'predict\.'$file_base'..\.tsv'|awk '{print "temp/temppredictions/"$0}'`|grep -v 'm'|cat predictions_header.tsv - > temp/predictions/predictions`basename $file`.tsv
         rm temp/temppredictions/*
      fi
   done
   #OUTPUT EXPRESSION
   gene=`basename $file|egrep "[0-9]*" -o|head -n 1`
   less /data1/antoine/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |grep -e '^T' -e $'\t5\t'$gene$'\t'>temp/expression/$gene.tsv
done
wait
kipoi_end=`date`
echo "started at "$start", bcf output ended at "$bcf_end", kipoi output ended at "$kipoi_end


#rm temp/fa_output/* temp/intervals/*
echo "starting R processing"
Rscript correlation_500kb.R
R_end=`date`
echo "finished!"
#rm temp/expression/* temp/predictions/*
echo "started at "$start", bcf output ended at "$bcf_end", kipoi output ended at "$kipoi_end", R output ended at "$R_end

#started at Thu Feb 21 19:52:40 EST 2019, bcf output ended at Thu Feb 21 23:03:27 EST 2019, kipoi output ended at Sat Feb 23 04:48:38 EST 2019, R output ended at Sat Feb 23 05:44:29 EST 2019
