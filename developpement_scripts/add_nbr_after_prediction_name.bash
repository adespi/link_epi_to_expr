for x in correlations_*.gz; do
   cp $x "$x~" &&
   echo -ne $x &&
   gzip -cd "$x~" |sed  -e 's/^\([^,][^,]*\)/\1(###)/'|awk '{gsub("###",NR-1,$0);print}'|gzip > $x ;
done

for x in correlations*/*.gz; do
   cp $x "$x~" &&
   echo -ne $x &&
   gzip -cd "$x~" |sed  -e 's/^\([^,][^,]*\)/\1(###)/'|awk '{gsub("###",NR-1,$0);print}'|gzip > $x ;
done
