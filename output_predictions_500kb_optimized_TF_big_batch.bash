#!/bin/bash
#avant de lancer le script : source activate kipoi-DeepSEA__predict
overall_start=`date`
#to lauch the script from any directory
cd `dirname $0`
#import configuration settings, to modify settings, create a copy of config as config.local in the git root directory and modify it as you want
#to update config during program, create config.update file with variables to change (beware, may cause problem)
if [ -e config.local ]; then
   . ./config.local
else
   . ./config
fi

if [ -e config.update ]; then
   rm config.update
fi

#loop to go through all the TFs in the list
for TF in `tac "gene_info_"$list_of_genes".txt" |sed '$d'|cut -f1`;do
   echo -ne "Analysing gene "`grep "$TF" "gene_info_"$list_of_genes".txt" |cut -f6`"   "`date`"\t\t\t\t\r"
   read -r chromosome gene <<<$(less ../data/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |cut -f1-5|grep $TF|cut -f3-4)
   gene_name=`grep "$TF" "gene_info_"$list_of_genes".txt" |cut -f6`
   if [ "$chromosome" = '' ] #not used anymore because of multistep process || [ -e 'correlations/correlations_'$chromosome'_'$gene'_'$gene_name'.csv.gz' ] || [ -e 'correlations/correlations_'$chromosome'_'$gene'_'$gene_name'.csv' ] || [ -e 'correlations_done/correlations_'$chromosome'_'$gene'_'$gene_name'.csv.gz' ] || [ -e 'correlations_done/correlations_'$chromosome'_'$gene'_'$gene_name'.csv' ]
   then
      continue
   fi
   echo "Computing for gene $gene_name chr $chromosome position $gene, window size = $window_size, window step = $window_step "`date`
   #new shell make in run in parallel and to filter the stderr in output screen
   echo $(
   #output fa for every patient (500kb around the interesting position)
   if [ "$output_fa" = true ] && [ ! -d temp/`echo $chromosome'_'$gene`/fa_output/ ] ; then
      for DIRECTORY in temp/ temp/`echo $chromosome'_'$gene`/ temp/`echo $chromosome'_'$gene`/intervals/ temp/`echo $chromosome'_'$gene`/fa_output/ temp/`echo $chromosome'_'$gene`/expression/ ;do
         if [ ! -d "$DIRECTORY" ]; then
            mkdir "$DIRECTORY";
         #else
            #rm "$DIRECTORY"*
         fi
      done
      echo -ne "outputing fa...           \r"
      start=`date`
      for patient in `more list_commun_patients`;do
         samtools faidx ../data/hs37d5.fa.gz $chromosome:`echo $gene $window_size |awk '{print $1<$2/2 ? 1 : $1-$2/2}'`-`echo $(($gene + $window_size/2 +1000))` | bcftools consensus -i 'type="snp"' -s $patient -H A ../data/genome_seq/ALL.chr$chromosome.phase3_shapeit2_mvncall_integrated_v*.20130502.genotypes.vcf.gz 2>/dev/null | sed "s/^>\(.*\)/>$patient:\1/" | bgzip > temp/`echo $chromosome'_'$gene`/fa_output/out`echo $chromosome'_'$gene'_'$patient`.fa.gz 2>/dev/null
      done
      fa_out_end=`date`
   fi

   if [[ "$compute_predictions_for_one_gene_expression" == true && -d temp/`echo $chromosome'_'$gene`/ && ! -e "correlations/correlations_"`echo $chromosome'_'$gene`"_"$gene_name".csv.gz" ]]; then
      #create intervals files for every patient (needed for kipoi prediction)
      seq_arg=`echo 1 $window_step $window_size`

      for file in temp/`echo $chromosome'_'$gene`/fa_output/out`echo $chromosome'_'$gene'_'`*.fa.gz;do
         name_seq=`bgzip -cd $file|head -n 1| tr -d ">"`
         for position in `seq $seq_arg`;do
            echo -e $name_seq"\t"$position"\t"$(($position+999))
         done > "temp/`echo $chromosome'_'$gene`/intervals/`basename $file .gz`"
      done
      create_intervals=`date`

      #fetch expression for the gene and put it into a file ready for python script
      less ../data/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |grep -e '^T' -e $'\t'$chromosome$'\t'$gene$'\t'>temp/`echo $chromosome'_'$gene`/expression/`echo $chromosome'_'$gene`.tsv


      if [ ! -d "correlations/" ]; then
         mkdir "correlations/";
      fi

      #predict epigenome marks using a python script, if one or more GPUs are avaivable, we use the GPU with the most RAM available
      echo -ne "predicting epigenome...           \r"
      if [ `command -v nvidia-smi` ]; then
         if [ `nvidia-smi -i 0 --query-gpu=utilization.memory --format=csv,noheader,nounits` -lt `nvidia-smi -i 1 --query-gpu=utilization.memory --format=csv,noheader,nounits` ];
         then
            CUDA_VISIBLE_DEVICES=0 python python_prediction.py `echo $chromosome'_'$gene` $(echo `seq $seq_arg|wc -l`) $gene_name $seq_arg $python_batch_size
         else
            CUDA_VISIBLE_DEVICES=1 python python_prediction.py `echo $chromosome'_'$gene` $(echo `seq $seq_arg|wc -l`) $gene_name $seq_arg $python_batch_size
         fi
      else
         python python_prediction.py `echo $chromosome'_'$gene` $(echo `seq $seq_arg|wc -l`) $gene_name $seq_arg $python_batch_size
      fi
      python_end=`date`
   fi


   if [[ "$convert_correlation_table_to_qvalues" == true && -e "correlations/correlations_"`echo $chromosome'_'$gene`"_"$gene_name".csv.gz" && ! -e correlations_small/correlations_$chromosome'_'$gene"_"$gene_name".csv.gz" ]]; then
      #convert correlation table to qvalues for each position
      if [ ! -d "correlations_small/" ]; then
         mkdir "correlations_small/";
      fi
      Rscript correlation_to_qvalue_single.R "correlations/correlations_"`echo $chromosome'_'$gene`"_"$gene_name".csv.gz"
   fi


   if [[ "$compute_predictions_for_all_gene_expression" == true && -e correlations_small/correlations_$chromosome'_'$gene"_"$gene_name".csv.gz" && ! -e "correlations/correlations_"$chromosome"_"$gene"_"$gene_name ]]; then
      if [ `less correlations_small/correlations_$chromosome'_'$gene"_"$gene_name".csv.gz"|wc -l` == 2 ]; then      
         rm temp/`echo $chromosome'_'$gene`/intervals/*
         #create intervals files for every patient only for interesting positions (needed for kipoi prediction)
         seq_arg=`echo 1 $window_step $window_size`


         positions_list=$(for x in `less correlations_small/correlations_$chromosome'_'$gene"_"$gene_name".csv.gz"|head -n 1|tr -d "\""|tr "," " "`;do
            i=`printf %.0f $(echo $x+1+$window_size/2|bc)`
            echo -e "\t"$i"\t"$(($i+999))","
         done)
         positions_list=`echo "${positions_list::-1}"`
         OLDIFS=$IFS
         for file in temp/`echo $chromosome'_'$gene`/fa_output/out`echo $chromosome'_'$gene'_'`*.fa.gz;do
            name_seq=`bgzip -cd $file|head -n 1| tr -d ">"`
            IFS=','
            for x in $positions_list;do
               IFS=$OLDIFS
               echo -e $name_seq$x
            done | tr " " "\t" > "temp/`echo $chromosome'_'$gene`/intervals/`basename $file .gz`"
            IFS=$OLDIFS
         done
         #name_seqs=`for file in 'temp/'$chromosome'_'$gene'/fa_output/out'$chromosome'_'$gene'_'*'.fa.gz'; do bgzip -cd $file |head -n 1| tr -d ">"; done`
         #for x in `less correlations_small/correlations_$chromosome'_'$gene"_"$gene_name".csv.gz"|head -n 1|tr -d "\""|tr "," " "`;do
         #   i=`printf %.0f $(echo $x+1+$window_size/2|bc)`
         #   j=1
         #   for file in temp/`echo $chromosome'_'$gene`/fa_output/out`echo $chromosome'_'$gene'_'`*.fa.gz;do
         #      name_seq=`echo $name_seqs | cut -d " " -f $j`
         #      echo -e $name_seq"\t"$i"\t"$(($i+999))>> "temp/`echo $chromosome'_'$gene`/intervals/`basename $file .gz`"
         #      j=$(($j+1))
         #   done
         #done
         create_intervals2=`date`
         if [ ! -d "correlations/correlations_"$chromosome"_"$gene"_"$gene_name ]; then
            mkdir "correlations/correlations_"$chromosome"_"$gene"_"$gene_name;
         fi
         #predict epigenome marks using a python script, if one or more GPUs are avaivable, we use the GPU with the most RAM available
         echo -ne "predicting epigenome...           \r"
         if [ `command -v nvidia-smi` ]; then
            if [ `nvidia-smi -i 0 --query-gpu=utilization.memory --format=csv,noheader,nounits` -lt `nvidia-smi -i 1 --query-gpu=utilization.memory --format=csv,noheader,nounits` ];
            then
               CUDA_VISIBLE_DEVICES=0 python python_prediction_multiple_genes.py `echo $chromosome'_'$gene` $(more "temp/`echo $chromosome'_'$gene`/intervals/`basename $file .gz`"|wc -l) $gene_name $seq_arg $python_batch_size $list_of_genes
            else
               CUDA_VISIBLE_DEVICES=1 python python_prediction_multiple_genes.py `echo $chromosome'_'$gene` $(more "temp/`echo $chromosome'_'$gene`/intervals/`basename $file .gz`"|wc -l) $gene_name $seq_arg $python_batch_size $list_of_genes
            fi
         else
            python python_prediction_multiple_genes.py `echo $chromosome'_'$gene` $(more "temp/`echo $chromosome'_'$gene`/intervals/`basename $file .gz`"|wc -l) $gene_name $seq_arg $python_batch_size $list_of_genes
         fi
         rm temp/`echo $chromosome'_'$gene`/intervals/*
      fi
   fi


   #if [[ "$convert_correlation_table_to_qvalues_all_gene_expression" == true && -d "correlations/correlations_"`echo $chromosome'_'$gene`"_"$gene_name ]]; then
   #   #convert correlation table to qvalues for each position
   #   if [ ! -d "correlations_small/correlations_"$chromosome"_"$gene"_"$gene_name ]; then
   #      mkdir "correlations_small/correlations_"$chromosome"_"$gene"_"$gene_name;
   #   fi
   #   for file in "correlations/correlations_"`echo $chromosome'_'$gene`"_"$gene_name/*.csv.gz;do
   #      Rscript correlation_to_qvalue_single.R $file
   #   done
   #fi


   #rm -r temp/`echo $chromosome'_'$gene`/
   #echo -e "started at "$start",\nfa_output_end at "$fa_out_end",\ncreate_intervals end at "$create_intervals",\npython ended at "$python_end>correlations/time_for_`echo $chromosome'_'$gene'_'$gene_name`.txt
   ) > /dev/null 2>&1 &
   sleep 5
   sleep $wait_time_between_two_batches

   #pause the script in order to have only $nbr_para_tasks jobs running at the same time
   while [ $(($nbr_para_tasks-1)) -lt `jobs -pr | grep "^[0-9]"| wc -l` ]
   do
      sleep 60
      if [ -e config.update ]; then
         . ./config.update
         rm config.update
      fi
   done
   if [ -e config.update ]; then
      . ./config.update
      rm config.update
   fi
done
wait
echo -e 'started the procedure at '$overall_start'\nended  the  procedure at '`date`
