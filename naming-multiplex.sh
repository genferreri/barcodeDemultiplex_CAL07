#!/bin/bash
NAME3=$(echo *.fastq.gz | awk -v FS='_' '{print $1}')

# Add name and index
for f in *.txt
 do NAME=${f%%_*} && T=${f%%.*} && INDEX=$(echo "$T".txt | grep -oP '(?<=_)\d+(?=\.)')
 echo $NAME;
 echo $INDEX;
 echo $T;
 awk -F, -v s="$NAME" '{$2="\t"s;print}' "$NAME"_"$INDEX".txt > "$NAME"_"$INDEX"-temp.txt
 tr -d '-' < "$NAME"_"$INDEX"-temp.txt > "$NAME"_"$INDEX"-tempII.txt
 rm "$NAME"_"$INDEX"-temp.txt
 mv "$NAME"_"$INDEX"-tempII.txt "$NAME"_"$INDEX".temp.txt
 awk -F, -v s="$INDEX" '{$3="\t"s;print}' "$NAME"_"$INDEX".temp.txt > "$NAME"_"$INDEX".tempII.txt
 rm "$NAME"_"$INDEX".temp.txt
 mv "$NAME"_"$INDEX".tempII.txt "$NAME"_"$INDEX"_final.txt
done

# This produces multi.txt 
rm "$NAME3"_??.txt

# Produces "multi", "multiPlot" and "Shannon"
Rscript ~/scripts/multiplex/Shannon.andPlot.r

##### Here we have produced filtered_final.txt
# Generate amino acid frequencies and logos plot
Rscript ~/scripts/multiplex/Translate.CAL07.barcode.R

rm Rplots.pdf

mv multi.plot.pdf "$NAME3".plot.pdf

# Add name and index to Shannon and name to multi
for i in *.txt
 do NAME=${f%%_*} && T=${f%%.*} && INDEX=$(echo "$T".txt | grep -oP '(?<=_)\d+(?=\.)')
 awk -F, -v s="$NAME" '{$2="\t"s;print}' Shannon > "$NAME"_Shannon-temp.txt
 tr -d '-' < "$NAME"_Shannon-temp.txt > "$NAME"_Shannon-tempII.txt
 rm "$NAME"_Shannon-temp.txt
 mv "$NAME"_Shannon-tempII.txt "$NAME"_Shannon-temp.txt
 # for multi
 awk -F, -v s="$NAME" '{$5="\t"s;print}' multi > "$NAME"_multi_temp.txt
 tr -d '-' < "$NAME"_multi_temp.txt > "$NAME"_"$INDEX"_multi_final.txt
 rm "$NAME"_multi_temp.txt
 # Add index to Shannon
 awk -F, -v s="$INDEX" '{$3="\t"s;print}' "$NAME"_Shannon-temp.txt > "$NAME"_Shannon_final.txt
 rm "$NAME"_Shannon-temp.txt

done

rm Shannon multi

 # Run Richness
 sh Richness.sh
 
 mv Richness.txt "$NAME3"_richness.txt

# AWK below subtracts 1 from the total cont which considers the header.
cat filtered_final.txt | wc -l | awk '{print $1-1}' > "$NAME3"_Count.txt

total=$(cat "$NAME3"_Count.txt)
awk -F, -v s="$total" '{$4="\t"s;print}' "$NAME"_richness.txt > "$NAME"_richness_final.txt
rm "$NAME"_richness.txt
rm "$NAME3"_Count.txt

# Fix name for filtered_final.txt
mv filtered_final.txt "$NAME3"_filtered_final.txt 

# mv files to folder
mkdir "$NAME3"_"$INDEX"/
mv -v *_final.txt "$NAME3"_"$INDEX"/
mv -v *_"$INDEX".fastq "$NAME3"_"$INDEX"/ 
mv -v "$NAME3".plot.pdf "$NAME3"_"$INDEX"/ 
mv -v "$NAME3"_richness.txt "$NAME3"_"$INDEX"/ 

mv -v Total.count.txt "$NAME3"_total.count.txt
mv -v "$NAME3"_total.count.txt "$NAME3"_"$INDEX"/

mv -v multi.aa.txt "$NAME3"_multi.aa.txt
mv -v "$NAME3"_multi.aa.txt "$NAME3"_"$INDEX"/

mv logos.plot.pdf "$NAME3"_logos.plot.pdf
mv "$NAME3"_logos.plot.pdf "$NAME3"_"$INDEX"/

mv filtered_final.txt "$NAME3"_filtered_final.txt
mv "$NAME3"_filtered_final.txt "$NAME3"_"$INDEX"/

