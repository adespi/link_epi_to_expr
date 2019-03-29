#!/bin/bash
#this scrip is the same as output_predictions_500kb_optimized_TF_big_batch.bash except that it skips the fa output part
#avant de lancer le script : source activate kipoi-DeepSEA__predict
overall_start=`date`
#to lauch the script from any directory
cd `dirname $0`
#import configuration settings, not synced with git, independant on every computer
. ./config

#loop to go through all the TFs in the list
for TF in `tac $list_of_genes |sed '$d'|cut -f1`;do
   echo -ne "Analysing gene "`grep "$TF" $list_of_genes |cut -f6`"   "`date`"\t\t\t\t\r"
   read -r chromosome gene <<<$(less ../data/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |cut -f1-5|grep $TF|cut -f3-4)
   gene_name=`grep "$TF" $list_of_genes |cut -f6`
   if [ "$chromosome" = '' ] || [ -e 'correlations_para/correlations_'$chromosome'_'$gene'_'$gene_name'.csv.gz' ] || [ -e 'correlations_para/correlations_'$chromosome'_'$gene'_'$gene_name'.csv' ] || [ -e 'correlations_done/correlations_'$chromosome'_'$gene'_'$gene_name'.csv.gz' ] || [ -e 'correlations_done/correlations_'$chromosome'_'$gene'_'$gene_name'.csv' ] || [ ! -d 'temp/'`echo $chromosome'_'$gene`'/' ]
   then
      continue
   fi
   echo "Computing for gene $gene_name chr $chromosome position $gene, window size = $window_size, window step = $window_step "`date`
   #new shell to filter the stderr
   echo $(for DIRECTORY in temp/ temp/`echo $chromosome'_'$gene`/ temp/`echo $chromosome'_'$gene`/intervals/ temp/`echo $chromosome'_'$gene`/fa_output/ temp/`echo $chromosome'_'$gene`/expression/ ;do
      if [ ! -d "$DIRECTORY" ]; then
         mkdir "$DIRECTORY";
      #else
         #rm "$DIRECTORY"*
      fi
   done

   echo -ne "outputing fa...           \r"
   start=`date`
   : <<'END_COMMENT'
   for patient in `more list_commun_patients`;do
      samtools faidx ../data/hs37d5.fa.gz $chromosome:`echo $gene $window_size |awk '{print $1<$2/2 ? 1 : $1-$2/2}'`-`echo $(($gene + $window_size/2 +1000))` | bcftools consensus -i 'type="snp"' -s $patient ../data/genome_seq/ALL.chr$chromosome.phase3_shapeit2_mvncall_integrated_v*.20130502.genotypes.vcf.gz | sed "s/^>\(.*\)/>$patient:\1/" | bgzip > temp/`echo $chromosome'_'$gene`/fa_output/out`echo $chromosome'_'$gene'_'$patient`.fa.gz 2>/dev/null
   done
   fa_out_end=`date`
   echo -ne "converting fa to kipoi input...           \r"
END_COMMENT
   seq_arg=`echo 1 $window_step $window_size`
   : <<'END_COMMENT'
   #for file in temp/`echo $chromosome'_'$gene`/fa_output/out`echo $chromosome'_'$gene'_'`*.fa;do
   #   samtools faidx $file
   #done
   

   for file in temp/`echo $chromosome'_'$gene`/fa_output/out`echo $chromosome'_'$gene'_'`*.fa.gz;do
      name_seq=`bgzip -cd $file|head -n 1| tr -d ">"`
      for position in `seq $seq_arg`;do
         echo -e $name_seq"\t"$position"\t"$(($position+999))
      done > "temp/`echo $chromosome'_'$gene`/intervals/`basename $file .gz`"
   done
   create_intervals=`date`
   less ../data/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |grep -e '^T' -e $'\t'$chromosome$'\t'$gene$'\t'>temp/`echo $chromosome'_'$gene`/expression/`echo $chromosome'_'$gene`.tsv
END_COMMENT
   if [ ! -d "correlations_para/" ]; then
      mkdir "correlations_para/";
   fi
   echo -ne "predicting epigenome...           \r"
   if [ `command -v nvidia-smi` ]; then
      if [ `nvidia-smi -i 0 --query-gpu=utilization.memory --format=csv,noheader,nounits` -lt `nvidia-smi -i 1 --query-gpu=utilization.memory --format=csv,noheader,nounits` ];
      then
         CUDA_VISIBLE_DEVICES=0 python python_prediction.py `echo $chromosome'_'$gene` $(echo `seq $seq_arg|wc -l`) $gene_name $seq_arg $python_batch_size
      else
         echo `echo $chromosome'_'$gene` $(echo `seq $seq_arg|wc -l`) $gene_name $seq_arg
         CUDA_VISIBLE_DEVICES=1 python python_prediction.py `echo $chromosome'_'$gene` $(echo `seq $seq_arg|wc -l`) $gene_name $seq_arg $python_batch_size
      fi
   else
      python python_prediction.py `echo $chromosome'_'$gene` $(echo `seq $seq_arg|wc -l`) $gene_name $seq_arg $python_batch_size
   fi
   python_end=`date`
   #rm -r temp/`echo $chromosome'_'$gene`/
   echo -e "started at "$start",\nfa_output_end at "$fa_out_end",\ncreate_intervals end at "$create_intervals",\npython ended at "$python_end>correlations_para/time_for_`echo $chromosome'_'$gene'_'$gene_name`.txt) > /dev/null 2>&1 &
   #sleep $wait_time_between_two_batches
   while [ $(($nbr_para_tasks-1)) -lt `jobs -p | grep "^[0-9]"| wc -l` ]
   do
      sleep 60
   done
done
wait
echo -e 'started the procedure at '$overall_start'\nended  the  procedure at '`date`
