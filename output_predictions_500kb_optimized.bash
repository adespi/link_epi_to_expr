#!/bin/bash
#avant de lancer le script : source activate kipoi-DeepSEA__predict
overall_start=`date`
cd ~/link_epi_to_expr/
#CREATE DIRECTORIES
for DIRECTORY in temp/ temp/intervals/  temp/predictions/ temp/fa_output/ temp/kipoi_input/ temp/expression/ temp/temppredictions/ temp/information_patient_seq/ temp/predictions_with_names/;do
   if [ ! -d "$DIRECTORY" ]; then
      mkdir "$DIRECTORY";
   #else
      #rm "$DIRECTORY"*
   fi
done

window_size=3000
window_step=100

#less nbr_snp_per_gene |sort -nr|less
#TO EXTRACT FA SEQ:
chromosome=1
for gene in `less nbr_snp_per_gene_500kb |sort -nr|sed "s/^ *//"|cut -f2 -d ' '|head -n 7`;do
   echo "gene $gene, window size = $window_size, window step = $window_step"
   echo -ne "outputing fa...           \r"
   start=`date`
   for patient in `more ~/link_epi_to_expr/list_commun_patients`;do
      samtools faidx /data1/antoine/hs37d5.fa.gz $chromosome:`echo $gene $window_size |awk '{print $1<$2/2 ? 1 : $1-$2/2}'`-`echo $gene $window_size |awk '{print $1+$2/2+1000}'` | bcftools consensus -i 'type="snp"' -s $patient /data1/antoine/genome_seq/ALL.chr$chromosome.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | sed "s/^>\(.*\)/>chr$patient:\1/" > ~/link_epi_to_expr/temp/fa_output/out`echo $chromosome'_'$gene'_'$patient`.fa &
      if [ 15 -lt `jobs -p | grep "^[0-9]"| wc -l` ]
      then
         wait
      fi
   done
   wait
   fa_out_end=`date`
   echo -ne "converting fa to kipoi input...           \r"
   seq_arg=`echo $window_size $window_step |awk '{print "1 "$2" "$1}'`
   for file in temp/fa_output/out`echo $chromosome'_'$gene'_'`*.fa;do
      samtools faidx $file
   done
   wait
   for position in `seq $seq_arg`;do
      for file in temp/fa_output/out`echo $chromosome'_'$gene'_'`*.fa;do
         name_seq=`head $file -n 1| tr -d ">"`
         samtools faidx $file $name_seq:$position-`echo $position|awk '{print ($1+999)}'`
      done|awk '/^>/ {printf("\n%s\t",$0);next; } { printf("%s",$0);}  END {printf("\n");}'|sed '1d'|sort -k 2|tee >(cut -f 2|uniq|awk '{printf("%s\n", ">"NR"\n"$0)}'> temp/kipoi_input/`echo $chromosome'_'$gene`_`echo $position $window_size|awk '{print $1-1-$2/2}'`.fa) >(cut -f 2|uniq -c|sed "s/^ *//"|cut -f1 -d ' '>temp/information_patient_seq/nbr_patient_per_seq_`echo $chromosome'_'$gene`_`echo $position $window_size |awk '{print $1-1-$2/2}'`) |cut -f 1 >temp/information_patient_seq/list_patient_ordered_by_seq_`echo $chromosome'_'$gene`_`echo $position $window_size |awk '{print $1-1-$2/2}'` &
      if [ 15 -lt `jobs -p | grep "^[0-9]"| wc -l` ]
      then
         wait
      fi
   done
   wait
   convert_fa_to_kipoi_input_end=`date`
   echo -ne "formating input...              \r"
   rm temp/fa_output/out`echo $chromosome'_'$gene'_*'`.fa
   rm temp/fa_output/out`echo $chromosome'_'$gene'_*'`.fa.fai
   for position in temp/kipoi_input/$chromosome'_'$gene_*.fa;do
      samtools faidx $position;
      if [ `more $position | wc -l` == 2 ] || [ `more "$position".fai|cut -f2|sort -n|head -n 1` -lt 1000 ]
      then
         rm $position $position.fai temp/information_patient_seq/list_patient_ordered_by_seq_`echo $chromosome'_'$gene`_`basename $position .fa` temp/information_patient_seq/nbr_patient_per_seq_`echo $chromosome'_'$gene`_`basename $position .fa`
      fi
   done
   for position in `seq $seq_arg`;do
      a=1
      kipoi_nbr=1
      file=temp/information_patient_seq/nbr_patient_per_seq_`echo $chromosome'_'$gene`_`echo $position $window_size |awk '{print $1-1-$2/2}'`
      if [ -e $file ]
      then
         file2=temp/information_patient_seq/list_patient_ordered_by_seq_`echo $chromosome'_'$gene`_`echo $position $window_size |awk '{print $1-1-$2/2}'`
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
   for file in temp/kipoi_input/`echo $chromosome'_'$gene`_*.fa;do
      more $file |grep "^>"|sed "s/^>//"|sed "s/$/\t1\t1000/"|split - -l 9 ~/link_epi_to_expr/temp/intervals/`basename $file`
      for x in temp/intervals/`basename $file`* ;do
         (kipoi predict DeepSEA/predict   --dataloader_args='{"intervals_file": '$x', "fasta_file": "'$file'", "num_chr_fasta": "False"}' -o 'temp/temppredictions/predict.'`basename $x`'.tsv' 2>&1 |cat > /dev/null
         rm $x )&
         if [ 15 -lt `jobs -p | grep "^[0-9]"| wc -l` ]
         then
            wait
         fi
      done
   done
   wait
   kipoi_end=`date`
   echo -ne "formating output...           \r"
   for file in temp/kipoi_input/`echo $chromosome'_'$gene`_*.fa;do
      #MERGE INTERVALS
      file_base=`basename $file`
      name_filess=`ls temp/temppredictions/ -1 | egrep 'predict\.'$file_base'..\.tsv'|awk '{print "temp/temppredictions/"$0}'`
      IFS=$'\n'
      for line in `cat $name_filess|grep -v 'm'`;do echo $line|cut -f1-4,6-|xargs printf "%.3g\t";printf '\n';done|cat predictions_header.tsv - > temp/predictions/predictions`basename $file`.tsv
      #cat `ls temp/temppredictions/ -1 | egrep 'predict\.'$file_base'..\.tsv'|awk '{print "temp/temppredictions/"$0}'`|grep -v 'm'|cat predictions_header.tsv - > temp/predictions/predictions`basename $file`.tsv
      gene_position=`basename $file '.fa'`
      IFS=$'\n'
      for patient in $(more temp/information_patient_seq/list_patient_ordered_by_seq_$gene_position);do
         more temp/predictions/predictions$gene_position.fa.tsv|grep -P $"^`echo $patient |cut -f2 -d ' '`\t"|sed "s/^`echo $patient |cut -f2 -d ' '`\t/`echo $patient |cut -f1 -d ' '`\t/"
      done|sort |cat predictions_header.tsv - > temp/predictions_with_names/predictions$gene_position.fa.tsv
      rm temp/predictions/predictions$gene_position.fa.tsv
      rm temp/information_patient_seq/list_patient_ordered_by_seq_$gene_position
      rm temp/temppredictions/predict.$gene_position.fa*.tsv
      unset IFS
   done
   wait
   rm temp/kipoi_input/`echo $chromosome'_'$gene`_*.fa
   rm temp/kipoi_input/`echo $chromosome'_'$gene`_*.fa.fai
   format_output_end=`date`
   #OUTPUT EXPRESSION
   #gene=`basename $file|egrep "[0-9]*" -o|head -n 1`
   less /data1/antoine/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |grep -e '^T' -e $'\t'$chromosome'\t'$gene$'\t'>temp/expression/`echo $chromosome'_'$gene`.tsv
   expression_output_end=`date`
   if [ ! -d "correlations" ]; then
         mkdir "correlations"
      fi
   #CALCULATE CORRELATIONS
   echo -ne "starting R processing...           \r"
   Rscript correlation_500kb_optimized.R $chromosome $gene
   R_end=`date`
   echo -ne "finished!                \r"
   rm temp/expression/`echo $chromosome'_'$gene`.tsv 
   echo temp/predictions_with_names/predictions`echo $chromosome'_'$gene`_*.fa.tsv|xargs rm
   echo -e "started at "$start",\nfa_output_end at "$fa_out_end",\nconvert_fa_to_kipoi_input_end at "$convert_fa_to_kipoi_input_end",\nformat_input_end at "$format_input_end",\nkipoi output ended at "$kipoi_end",\nformat_output_end at "$format_output_end",\nexpression_output_end at "$expression_output_end",\nR output ended at "$R_end>correlations/time_for_`echo $chromosome'_'$gene`.txt
done
echo 'started the procedure at '$overall_start'ended  the  procedure at '`date`
#rm `ls temp/temppredictions/ -1 | egrep 'predict\.'$file_base'..\.tsv'|awk '{print "temp/temppredictions/"$0}'`


########more testest.csv |tail -n1|cut -f1-7|sed "s/\([0-9\.][0-9\.]*\)/`printf \1`/g"
########more testest.csv |tail -n1|cut -f1-4,6-15|xargs printf "%.2g\t"
########for line in `cat <(cat testest.csv)`;do echo $line|cut -f1-4,6-15|xargs printf "%.2g\t";printf '\n';done
########cat `ls temp/temppredictions/ -1 | egrep 'predict\.'$file_base'..\.tsv'|awk '{print "temp/temppredictions/"$0}'`|grep -v 'm'|cat predictions_header.tsv - > temp/predictions/predictions`basename $file`.tsv
