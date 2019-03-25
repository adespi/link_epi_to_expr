#!/bin/bash
#avant de lancer le script : source activate kipoi-DeepSEA__predict
overall_start=`date`
#cd ~/link_epi_to_expr/
#CREATE DIRECTORIES

window_size=500000
window_step=100

#less nbr_snp_per_gene |sort -nr|less
#TO EXTRACT FA SEQ:
for TF in `tac gene_info_some_genes.txt |sed '$d'|cut -f1`;do
   echo "Computing for gene "`grep "$TF" gene_info_some_genes.txt |cut -f6`"   "`date`
   echo $(read -r chromosome gene <<<$(less /data1/antoine/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |cut -f1-5|grep $TF|cut -f3-4)
   gene_name=`grep "$TF" gene_info_some_genes.txt |cut -f6`
   if [ "$chromosome" = '' ] || [ -e 'correlations_para/correlations_'$chromosome'_'$gene'_'$gene_name'.csv' ]
   then
      continue
   fi
   #for gene in `less nbr_snp_per_gene_500kb |sort -nr|sed "s/^ *//"|cut -f2 -d ' '|head -n 7`;do
   echo "gene $gene_name chr $chromosome position $gene, window size = $window_size, window step = $window_step"
   for DIRECTORY in temp/ temp/`echo $chromosome'_'$gene`/ temp/`echo $chromosome'_'$gene`/intervals/  temp/`echo $chromosome'_'$gene`/predictions/ temp/`echo $chromosome'_'$gene`/fa_output/ temp/`echo $chromosome'_'$gene`/kipoi_input/ temp/`echo $chromosome'_'$gene`/expression/ temp/`echo $chromosome'_'$gene`/temppredictions/ temp/`echo $chromosome'_'$gene`/information_patient_seq/ temp/`echo $chromosome'_'$gene`/predictions_with_names/;do
      if [ ! -d "$DIRECTORY" ]; then
         mkdir "$DIRECTORY";
      #else
         #rm "$DIRECTORY"*
      fi
   done

   echo -ne "outputing fa...           \r"
   start=`date`
   for patient in `more list_commun_patients`;do
      samtools faidx /data1/antoine/hs37d5.fa.gz $chromosome:`echo $gene $window_size |awk '{print $1<$2/2 ? 1 : $1-$2/2}'`-`echo $gene $window_size |awk '{print $1+$2/2+1000}'` | bcftools consensus -i 'type="snp"' -s $patient /data1/antoine/genome_seq/ALL.chr$chromosome.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | sed "s/^>\(.*\)/>chr$patient:\1/" > temp/`echo $chromosome'_'$gene`/fa_output/out`echo $chromosome'_'$gene'_'$patient`.fa 
   done
   fa_out_end=`date`
   echo -ne "converting fa to kipoi input...           \r"
   seq_arg=`echo $window_size $window_step |awk '{print "1 "$2" "$1}'`
   for file in temp/`echo $chromosome'_'$gene`/fa_output/out`echo $chromosome'_'$gene'_'`*.fa;do
      samtools faidx $file
   done
   for position in `seq $seq_arg`;do
      for file in temp/`echo $chromosome'_'$gene`/fa_output/out`echo $chromosome'_'$gene'_'`*.fa;do
         name_seq=`head $file -n 1| tr -d ">"`
         samtools faidx $file $name_seq:$position-`echo $position|awk '{print ($1+999)}'`
      done|awk '/^>/ {printf("\n%s\t",$0);next; } { printf("%s",$0);}  END {printf("\n");}'|sed '1d'|sort -k 2|tee >(cut -f 2|uniq|awk '{printf("%s\n", ">"NR"\n"$0)}'> temp/`echo $chromosome'_'$gene`/kipoi_input/`echo $chromosome'_'$gene`_`echo $position $window_size|awk '{print $1-1-$2/2}'`.fa) >(cut -f 2|uniq -c|sed "s/^ *//"|cut -f1 -d ' '>temp/`echo $chromosome'_'$gene`/information_patient_seq/nbr_patient_per_seq_`echo $chromosome'_'$gene`_`echo $position $window_size |awk '{print $1-1-$2/2}'`) |cut -f 1 >temp/`echo $chromosome'_'$gene`/information_patient_seq/list_patient_ordered_by_seq_`echo $chromosome'_'$gene`_`echo $position $window_size |awk '{print $1-1-$2/2}'`
   done
   convert_fa_to_kipoi_input_end=`date`
   echo -ne "formating input...              \r"
   rm temp/`echo $chromosome'_'$gene`/fa_output/out`echo $chromosome'_'$gene'_*'`.fa
   rm temp/`echo $chromosome'_'$gene`/fa_output/out`echo $chromosome'_'$gene'_*'`.fa.fai
   for position in temp/`echo $chromosome'_'$gene`/kipoi_input/$chromosome'_'$gene_*.fa;do
      samtools faidx $position;
      if [ `more $position | wc -l` == 2 ] || [ `more "$position".fai|cut -f2|sort -n|head -n 1` -lt 1000 ]
      then
         rm $position $position.fai temp/`echo $chromosome'_'$gene`/information_patient_seq/list_patient_ordered_by_seq_`basename $position .fa` temp/`echo $chromosome'_'$gene`/information_patient_seq/nbr_patient_per_seq_`basename $position .fa`
      fi
   done
   for position in `seq $seq_arg`;do
      a=1
      kipoi_nbr=1
      file=temp/`echo $chromosome'_'$gene`/information_patient_seq/nbr_patient_per_seq_`echo $chromosome'_'$gene`_`echo $position $window_size |awk '{print $1-1-$2/2}'`
      if [ -e $file ]
      then
         file2=temp/`echo $chromosome'_'$gene`/information_patient_seq/list_patient_ordered_by_seq_`echo $chromosome'_'$gene`_`echo $position $window_size |awk '{print $1-1-$2/2}'`
         for seq in `more $file`;do
            sed -i "$a,`echo $a $seq|awk '{print $1+$2-1}'` s/$/ $kipoi_nbr/" $file2
            a=$((a+seq))
            ((kipoi_nbr++))
         done
         rm $file
      fi
   done
   format_input_end=`date`
   echo -ne "predicting epigenome...           \r"
   for file in temp/`echo $chromosome'_'$gene`/kipoi_input/`echo $chromosome'_'$gene`_*.fa;do
      more $file |grep "^>"|sed "s/^>//"|sed "s/$/\t1\t1000/"|split - -l 9 temp/`echo $chromosome'_'$gene`/intervals/`basename $file`
      for x in temp/`echo $chromosome'_'$gene`/intervals/`basename $file`* ;do
         (kipoi predict DeepSEA/predict   --dataloader_args='{"intervals_file": '$x', "fasta_file": "'$file'", "num_chr_fasta": "False"}' -o 'temp/'`echo $chromosome'_'$gene`'/temppredictions/predict.'`basename $x`'.tsv' 2>&1 |cat > /dev/null
         rm $x )
      done
   done
   kipoi_end=`date`
   echo -ne "formating output...           \r"
   for file in temp/`echo $chromosome'_'$gene`/kipoi_input/`echo $chromosome'_'$gene`_*.fa;do
      #MERGE INTERVALS
      file_base=`basename $file`
      name_filess=`ls 'temp/'$chromosome'_'$gene'/temppredictions/' -1 | egrep 'predict\.'$file_base'..\.tsv'|awk -vfold="$chromosome"_"$gene" '{print "temp/"fold"/temppredictions/"$0}'`
      IFS=$'\n'
      for line in `cat $name_filess|grep -v 'm'`;do echo $line|cut -f1-4,6-|xargs printf "%.3g\t";printf '\n';done|cat predictions_header.tsv - > temp/`echo $chromosome'_'$gene`/predictions/predictions`basename $file`.tsv
      #cat `ls temp/temppredictions/ -1 | egrep 'predict\.'$file_base'..\.tsv'|awk '{print "temp/temppredictions/"$0}'`|grep -v 'm'|cat predictions_header.tsv - > temp/predictions/predictions`basename $file`.tsv
      gene_position=`basename $file '.fa'`
      IFS=$'\n'
      for patient in $(more temp/`echo $chromosome'_'$gene`/information_patient_seq/list_patient_ordered_by_seq_$gene_position);do
         more temp/`echo $chromosome'_'$gene`/predictions/predictions$gene_position.fa.tsv|grep -P $"^`echo $patient |cut -f2 -d ' '`\t"|sed "s/^`echo $patient |cut -f2 -d ' '`\t/`echo $patient |cut -f1 -d ' '`\t/"
      done|sort |cat predictions_header.tsv - > temp/`echo $chromosome'_'$gene`/predictions_with_names/predictions$gene_position.fa.tsv
      rm temp/`echo $chromosome'_'$gene`/predictions/predictions$gene_position.fa.tsv
      rm temp/`echo $chromosome'_'$gene`/information_patient_seq/list_patient_ordered_by_seq_$gene_position
      rm temp/`echo $chromosome'_'$gene`/temppredictions/predict.$gene_position.fa*.tsv
      unset IFS
   done
   rm temp/`echo $chromosome'_'$gene`/kipoi_input/`echo $chromosome'_'$gene`_*.fa
   rm temp/`echo $chromosome'_'$gene`/kipoi_input/`echo $chromosome'_'$gene`_*.fa.fai
   format_output_end=`date`
   #OUTPUT EXPRESSION
   #gene=`basename $file|egrep "[0-9]*" -o|head -n 1`
   less /data1/antoine/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |grep -e '^T' -e $'\t'$chromosome$'\t'$gene$'\t'>temp/`echo $chromosome'_'$gene`/expression/`echo $chromosome'_'$gene`.tsv
   expression_output_end=`date`
   if [ ! -d "correlations_para/" ]; then
         mkdir "correlations_para/"
      fi
   #CALCULATE CORRELATIONS
   echo -ne "starting R processing...           \r"
   Rscript correlation_500kb_optimized_para.R $chromosome $gene $gene_name
   R_end=`date`
   echo -ne "finished!                \r"
   #rm temp/`echo $chromosome'_'$gene`/expression/`echo $chromosome'_'$gene`.tsv 
   #echo temp/`echo $chromosome'_'$gene`/predictions_with_names/predictions`echo $chromosome'_'$gene`_*.fa.tsv|xargs rm
   #rm -r temp/`echo $chromosome'_'$gene`/
   echo -e "started at "$start",\nfa_output_end at "$fa_out_end",\nconvert_fa_to_kipoi_input_end at "$convert_fa_to_kipoi_input_end",\nformat_input_end at "$format_input_end",\nkipoi output ended at "$kipoi_end",\nformat_output_end at "$format_output_end",\nexpression_output_end at "$expression_output_end",\nR output ended at "$R_end>correlations_para/time_for_`echo $chromosome'_'$gene'_'$gene_name`.txt) > /dev/null 2>&1 &
   sleep 60
   while [ 15 -lt `jobs -p | grep "^[0-9]"| wc -l` ]
   do
      sleep 60
   done
done
wait
echo -e 'started the procedure at '$overall_start'\nended  the  procedure at '`date`
#rm `ls temp/temppredictions/ -1 | egrep 'predict\.'$file_base'..\.tsv'|awk '{print "temp/temppredictions/"$0}'`


########more testest.csv |tail -n1|cut -f1-7|sed "s/\([0-9\.][0-9\.]*\)/`printf \1`/g"
########more testest.csv |tail -n1|cut -f1-4,6-15|xargs printf "%.2g\t"
########for line in `cat <(cat testest.csv)`;do echo $line|cut -f1-4,6-15|xargs printf "%.2g\t";printf '\n';done
########cat `ls temp/temppredictions/ -1 | egrep 'predict\.'$file_base'..\.tsv'|awk '{print "temp/temppredictions/"$0}'`|grep -v 'm'|cat predictions_header.tsv - > temp/predictions/predictions`basename $file`.tsv
