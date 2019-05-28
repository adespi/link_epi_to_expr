#Before Starting, always :
#source activate kipoi-DeepSEA__predict

#WORK ON SEQ:
#bcftools view -G /data1/antoine/genome_seq/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz|less
#less /data1/antoine/hs37d5.fa.gz |grep 'TGTGTTTCCCGATC' -m 3
#samtools faidx /data1/antoine/hs37d5.fa.gz 1:10176-10236 | bcftools consensus /data1/antoine/genome_seq/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -o out.fa
#bcftools consensus -s 'HG00096' -f /data1/antoine/hs37d5.fa.gz /data1/antoine/genome_seq/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > out.fa
#samtools faidx /data1/antoine/hs37d5.fa.gz 1:10176-10236 |bcftools consensus -s 'HG00096' /data1/antoine/genome_seq/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > out.fa

#bcftools view -G /data1/antoine/genome_seq/#ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '^#'|cut -f1-7| less

#WORK ON EXPR
#less /data1/antoine/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |cut -f1-5|sort -k3,3n -k4,4n|grep $'\t1\t'|less

#LIST COMMON PATIENTS:
#for expr in `more temp/expr.txt`;do grep $expr temp/genes.txt -o ;done|less>list_commun_patients

#LIST OF GENE POSITIONS:
#less /data1/antoine/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |cut -f1-5|sort -k3,3n -k4,4n|grep $'\t1\t'|cut -f4> list_gene_positions_chr1

#LIST OF CODING GENES (TF)
#less Network_compendium/Other_networks/Protein-protein_interaction/biogrid-3.2.116.txt.gz |cut -f1|sort|uniq>name_of_all_genes
#less Network_compendium/*/*/*.gz |cut -f1|sort|uniq>name_of_all_genes_bigger






#RM INTERVALS & PREDICTIONS
rm intervals/* predictions/*

#less nbr_snp_per_gene |sort -nr|less
#TO EXTRACT FA SEQ:
for patient in `more ~/link_epi_to_expr/list_commun_patients`;do samtools faidx /data1/antoine/hs37d5.fa.gz 1:568526-569626 |bcftools consensus -s $patient /data1/antoine/genome_seq/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz|sed "s/^>\(.*\)/>chr$patient:\1/";done>~/link_epi_to_expr/out.fa
#COUNT NBR OF BP
#for patient in `more ~/link_epi_to_expr/list_commun_patients`;do samtools faidx /data1/antoine/hs37d5.fa.gz 1:568526-569626 |bcftools consensus -s $patient /data1/antoine/genome_seq/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz|sed "s/^>\(.*\)/>chr$patient:\1/"|grep -v '^>'|tr -d "\n"|wc;done | sort -n|uniq -c

#COUNT NBR OF PREDICTIONS TO DO
#for file in temp/fa_output/out*.fa;do
#more $file|awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'|sed '1d'|#grep -v "^>"|sort|uniq|wc -l
#done|sort|uniq -c



#TO CREATE INTERVALS FILE:
#more ~/link_epi_to_expr/out.fa |grep "^>"|sed "s/^>//"|sed "s/$/\t1\t1000/">~/link_epi_to_expr/intervals.inter
#TO SPLIT INTERVALS:
#split intervals.inter -l 5 intervals/
#ALL AT ONCE
more ~/link_epi_to_expr/out.fa |grep "^>"|sed "s/^>//"|sed "s/$/\t1\t1000/"|split - -l 9 ~/link_epi_to_expr/intervals/XX


#TO PREDICT EPI:
#kipoi predict DeepSEA/predict   --dataloader_args='{"intervals_file": "intervals.inter", "fasta_file": "out.fa"}' -o '/tmp/DeepSEA|predict.example_pred.tsv'
#WARNINNG: at prediction n5 there are some problems, stops at prediction n8
#PREDICT SEVERAL INTERVALS
for x in intervals/* ;do kipoi predict DeepSEA/predict   --dataloader_args='{"intervals_file": '$x', "fasta_file": "out.fa"}' -o 'predictions/predict.'`basename $x`'.tsv'; done
#MERGE INTERVALS
cat predictions/*|grep -v 'm'> temp/predictions.tsv

for line in `more temp/predictions.tsv`;do 
awk -v var=`cut -c4-10 $line` '{for(n=1;n<=NF;n++)if($n == var)print n}' predictions_header.tsv
done

awk '{print $0;OUI}' fichier

#COUNT ACCURACY if n times the same prediction to make
#uniq -f3 predictions.tsv -c
cat predictions_header.tsv temp/predictions.tsv > predictions.tsv

#OUTPUT EXPRESSION
less /data1/antoine/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |grep -e '^T' -e $'\t1\t569076\t'>expression/569076.tsv
