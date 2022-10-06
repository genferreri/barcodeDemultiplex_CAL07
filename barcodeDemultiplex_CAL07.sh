#!/bin/bash

NAME1=$(find . -maxdepth 1 -name "*_S*_L001_R1_001.fastq.gz" | xargs -I {} basename {} .fastq.gz)
NAME2=$(find . -maxdepth 1 -name "*_S*_L001_R2_001.fastq.gz" | xargs -I {} basename {} .fastq.gz)
NAME3=$(echo *.fastq.gz | awk -v FS='_' '{print $1}')

# Evaluate quality of fastq files
fastqc "$NAME1".fastq.gz "$NAME2".fastq.gz
mkdir 00_Original 
mv *.zip 00_Original 
mv *.html 00_Original 

#filter reads Q30
bbduk.sh in1="$NAME1".fastq.gz in2="$NAME2".fastq.gz out1="$NAME1"_clean.fastq.gz out2="$NAME2"_clean.fastq.gz qtrim=rl trimq=30

# Evaluate clean fastq files
fastqc "$NAME1"_clean.fastq.gz "$NAME2"_clean.fastq.gz
mkdir 01_Clean_fastq
mv *.zip 01_Clean_fastq
mv *.html 01_Clean_fastq

# Conserve original file
mv "$NAME1".fastq.gz 00_Original
mv "$NAME2".fastq.gz 00_Original

# Rename "cleaned" file
mv "$NAME1"_clean.fastq.gz "$NAME1".fastq.gz
mv "$NAME2"_clean.fastq.gz "$NAME2".fastq.gz

# Make reverse complement for R2
seqtk seq -r "$NAME2".fastq.gz > "$NAME2".rev.fastq.gz

# Parse barcode
# 1
zgrep ^CTACACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_01.fastq
zgrep ^CTACACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_01.fastq
cat "$NAME1"_01.fastq "$NAME2"_01.fastq > "$NAME3"_01.fastq

grep AGTTC "$NAME3"_01.fastq > "$NAME3"_01.clean.fastq
mv "$NAME3"_01.clean.fastq "$NAME3"_01.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_01.fastq > "$NAME3"_01-clean.txt
rm "$NAME3"_01.fastq
mv "$NAME3"_01-clean.txt "$NAME3"_01.txt

sh naming-multiplex.sh

# Parse barcode
# 2
zgrep ^CAATACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_02.fastq
zgrep ^CAATACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_02.fastq
cat "$NAME1"_02.fastq "$NAME2"_02.fastq > "$NAME3"_02.fastq

grep AGTTC "$NAME3"_02.fastq > "$NAME3"_02.clean.fastq
mv "$NAME3"_02.clean.fastq "$NAME3"_02.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_02.fastq > "$NAME3"_02-clean.txt
rm "$NAME3"_02.fastq
mv "$NAME3"_02-clean.txt "$NAME3"_02.txt

sh naming-multiplex.sh

# Parse barcode
# 3
zgrep ^TTATACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_03.fastq
zgrep ^TTATACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_03.fastq
cat "$NAME1"_03.fastq "$NAME2"_03.fastq > "$NAME3"_03.fastq

grep AGTTC "$NAME3"_03.fastq > "$NAME3"_03.clean.fastq
mv "$NAME3"_03.clean.fastq "$NAME3"_03.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_03.fastq > "$NAME3"_03-clean.txt
rm "$NAME3"_03.fastq
mv "$NAME3"_03-clean.txt "$NAME3"_03.txt

sh naming-multiplex.sh

# Parse barcode
# 4
zgrep ^TGCCACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_04.fastq
zgrep ^TGCCACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_04.fastq
cat "$NAME1"_04.fastq "$NAME2"_04.fastq > "$NAME3"_04.fastq

grep AGTTC "$NAME3"_04.fastq > "$NAME3"_04.clean.fastq
mv "$NAME3"_04.clean.fastq "$NAME3"_04.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_04.fastq > "$NAME3"_04-clean.txt
rm "$NAME3"_04.fastq
mv "$NAME3"_04-clean.txt "$NAME3"_04.txt

sh naming-multiplex.sh

# Parse barcode
# 5
zgrep ^TTGCACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_05.fastq
zgrep ^TTGCACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_05.fastq
cat "$NAME1"_05.fastq "$NAME2"_05.fastq > "$NAME3"_05.fastq

grep AGTTC "$NAME3"_05.fastq > "$NAME3"_05.clean.fastq
mv "$NAME3"_05.clean.fastq "$NAME3"_05.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_05.fastq > "$NAME3"_05-clean.txt
rm "$NAME3"_05.fastq
mv "$NAME3"_05-clean.txt "$NAME3"_05.txt

sh naming-multiplex.sh

# Parse barcode
# 6
zgrep ^GACAACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_06.fastq
zgrep ^GACAACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_06.fastq
cat "$NAME1"_06.fastq "$NAME2"_06.fastq > "$NAME3"_06.fastq

grep AGTTC "$NAME3"_06.fastq > "$NAME3"_06.clean.fastq
mv "$NAME3"_06.clean.fastq "$NAME3"_06.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_06.fastq > "$NAME3"_06-clean.txt
rm "$NAME3"_06.fastq
mv "$NAME3"_06-clean.txt "$NAME3"_06.txt

sh naming-multiplex.sh

# Parse barcode
# 7
zgrep ^GACTACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_07.fastq
zgrep ^GACTACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_07.fastq
cat "$NAME1"_07.fastq "$NAME2"_07.fastq > "$NAME3"_07.fastq

grep AGTTC "$NAME3"_07.fastq > "$NAME3"_07.clean.fastq
mv "$NAME3"_07.clean.fastq "$NAME3"_07.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_07.fastq > "$NAME3"_07-clean.txt
rm "$NAME3"_07.fastq
mv "$NAME3"_07-clean.txt "$NAME3"_07.txt

sh naming-multiplex.sh

# Parse barcode
# 8
zgrep ^TCCTACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_08.fastq
zgrep ^TCCTACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_08.fastq
cat "$NAME1"_08.fastq "$NAME2"_08.fastq > "$NAME3"_08.fastq

grep AGTTC "$NAME3"_08.fastq > "$NAME3"_08.clean.fastq
mv "$NAME3"_08.clean.fastq "$NAME3"_08.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_08.fastq > "$NAME3"_08-clean.txt
rm "$NAME3"_08.fastq
mv "$NAME3"_08-clean.txt "$NAME3"_08.txt

sh naming-multiplex.sh

# Parse barcode
# 9
zgrep ^ACTCACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_09.fastq
zgrep ^ACTCACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_09.fastq
cat "$NAME1"_09.fastq "$NAME2"_09.fastq > "$NAME3"_09.fastq

grep AGTTC "$NAME3"_09.fastq > "$NAME3"_09.clean.fastq
mv "$NAME3"_09.clean.fastq "$NAME3"_09.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_09.fastq > "$NAME3"_09-clean.txt
rm "$NAME3"_09.fastq
mv "$NAME3"_09-clean.txt "$NAME3"_09.txt

sh naming-multiplex.sh

# Parse barcode
# 10
zgrep ^ACACACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_10.fastq
zgrep ^ACACACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_10.fastq
cat "$NAME1"_10.fastq "$NAME2"_10.fastq > "$NAME3"_10.fastq

grep AGTTC "$NAME3"_10.fastq > "$NAME3"_10.clean.fastq
mv "$NAME3"_10.clean.fastq "$NAME3"_10.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_10.fastq > "$NAME3"_10-clean.txt
rm "$NAME3"_10.fastq
mv "$NAME3"_10-clean.txt "$NAME3"_10.txt

sh naming-multiplex.sh

# Parse barcode
# 11
zgrep ^TCCCACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_11.fastq
zgrep ^TCCCACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_11.fastq
cat "$NAME1"_11.fastq "$NAME2"_11.fastq > "$NAME3"_11.fastq

grep AGTTC "$NAME3"_11.fastq > "$NAME3"_11.clean.fastq
mv "$NAME3"_11.clean.fastq "$NAME3"_11.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_11.fastq > "$NAME3"_11-clean.txt
rm "$NAME3"_11.fastq
mv "$NAME3"_11-clean.txt "$NAME3"_11.txt

sh naming-multiplex.sh

# Parse barcode
# 12
zgrep ^TCACACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_12.fastq
zgrep ^TCACACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_12.fastq
cat "$NAME1"_12.fastq "$NAME2"_12.fastq > "$NAME3"_12.fastq

grep AGTTC "$NAME3"_12.fastq > "$NAME3"_12.clean.fastq
mv "$NAME3"_12.clean.fastq "$NAME3"_12.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_12.fastq > "$NAME3"_12-clean.txt
rm "$NAME3"_12.fastq
mv "$NAME3"_12-clean.txt "$NAME3"_12.txt

sh naming-multiplex.sh

# Parse barcode
# 13
zgrep ^ATTGACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_13.fastq
zgrep ^ATTGACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_13.fastq
cat "$NAME1"_13.fastq "$NAME2"_13.fastq > "$NAME3"_13.fastq

grep AGTTC "$NAME3"_13.fastq > "$NAME3"_13.clean.fastq
mv "$NAME3"_13.clean.fastq "$NAME3"_13.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_13.fastq > "$NAME3"_13-clean.txt
rm "$NAME3"_13.fastq
mv "$NAME3"_13-clean.txt "$NAME3"_13.txt

sh naming-multiplex.sh

# Parse barcode
# 14
zgrep ^AGCTACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_14.fastq
zgrep ^AGCTACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_14.fastq
cat "$NAME1"_14.fastq "$NAME2"_14.fastq > "$NAME3"_14.fastq

grep AGTTC "$NAME3"_14.fastq > "$NAME3"_14.clean.fastq
mv "$NAME3"_14.clean.fastq "$NAME3"_14.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_14.fastq > "$NAME3"_14-clean.txt
rm "$NAME3"_14.fastq
mv "$NAME3"_14-clean.txt "$NAME3"_14.txt

sh naming-multiplex.sh

# Parse barcode
# 15
zgrep ^TCTCACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_15.fastq
zgrep ^TCTCACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_15.fastq
cat "$NAME1"_15.fastq "$NAME2"_15.fastq > "$NAME3"_15.fastq

grep AGTTC "$NAME3"_15.fastq > "$NAME3"_15.clean.fastq
mv "$NAME3"_15.clean.fastq "$NAME3"_15.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_15.fastq > "$NAME3"_15-clean.txt
rm "$NAME3"_15.fastq
mv "$NAME3"_15-clean.txt "$NAME3"_15.txt

sh naming-multiplex.sh

# Parse barcode
# 16
zgrep ^GATGACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_16.fastq
zgrep ^GATGACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_16.fastq
cat "$NAME1"_16.fastq "$NAME2"_16.fastq > "$NAME3"_16.fastq

grep AGTTC "$NAME3"_16.fastq > "$NAME3"_16.clean.fastq
mv "$NAME3"_16.clean.fastq "$NAME3"_16.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_16.fastq > "$NAME3"_16-clean.txt
rm "$NAME3"_16.fastq
mv "$NAME3"_16-clean.txt "$NAME3"_16.txt

sh naming-multiplex.sh

# Parse barcode
# 17
zgrep ^GGCGACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_17.fastq
zgrep ^GGCGACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_17.fastq
cat "$NAME1"_17.fastq "$NAME2"_17.fastq > "$NAME3"_17.fastq

grep AGTTC "$NAME3"_17.fastq > "$NAME3"_17.clean.fastq
mv "$NAME3"_17.clean.fastq "$NAME3"_17.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_17.fastq > "$NAME3"_17-clean.txt
rm "$NAME3"_17.fastq
mv "$NAME3"_17-clean.txt "$NAME3"_17.txt

sh naming-multiplex.sh

# Parse barcode
# 18

zgrep ^ATTAACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_18.fastq
zgrep ^ATTAACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_18.fastq
cat "$NAME1"_18.fastq "$NAME2"_18.fastq > "$NAME3"_18.fastq

grep AGTTC "$NAME3"_18.fastq > "$NAME3"_18.clean.fastq
mv "$NAME3"_18.clean.fastq "$NAME3"_18.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_18.fastq > "$NAME3"_18-clean.txt
rm "$NAME3"_18.fastq
mv "$NAME3"_18-clean.txt "$NAME3"_18.txt

sh naming-multiplex.sh

# Parse barcode
# 19

zgrep ^GCTTACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_19.fastq
zgrep ^GCTTACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_19.fastq
cat "$NAME1"_19.fastq "$NAME2"_19.fastq > "$NAME3"_19.fastq

grep AGTTC "$NAME3"_19.fastq > "$NAME3"_19.clean.fastq
mv "$NAME3"_19.clean.fastq "$NAME3"_19.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_19.fastq > "$NAME3"_19-clean.txt
rm "$NAME3"_19.fastq
mv "$NAME3"_19-clean.txt "$NAME3"_19.txt

sh naming-multiplex.sh

# Parse barcode
# 20

zgrep ^GACCACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_20.fastq
zgrep ^GACCACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_20.fastq
cat "$NAME1"_20.fastq "$NAME2"_20.fastq > "$NAME3"_20.fastq

grep AGTTC "$NAME3"_20.fastq > "$NAME3"_20.clean.fastq
mv "$NAME3"_20.clean.fastq "$NAME3"_20.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_20.fastq > "$NAME3"_20-clean.txt
rm "$NAME3"_20.fastq
mv "$NAME3"_20-clean.txt "$NAME3"_20.txt

sh naming-multiplex.sh

# Parse barcode
# 21

zgrep ^CGCTACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_21.fastq
zgrep ^CGCTACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_21.fastq
cat "$NAME1"_21.fastq "$NAME2"_21.fastq > "$NAME3"_21.fastq

grep AGTTC "$NAME3"_21.fastq > "$NAME3"_21.clean.fastq
mv "$NAME3"_21.clean.fastq "$NAME3"_21.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_21.fastq > "$NAME3"_21-clean.txt
rm "$NAME3"_21.fastq
mv "$NAME3"_21-clean.txt "$NAME3"_21.txt

sh naming-multiplex.sh

# Parse barcode
# 22

zgrep ^CTGTACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_22.fastq
zgrep ^CTGTACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_22.fastq
cat "$NAME1"_22.fastq "$NAME2"_22.fastq > "$NAME3"_22.fastq

grep AGTTC "$NAME3"_22.fastq > "$NAME3"_22.clean.fastq
mv "$NAME3"_22.clean.fastq "$NAME3"_22.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_22.fastq > "$NAME3"_22-clean.txt
rm "$NAME3"_22.fastq
mv "$NAME3"_22-clean.txt "$NAME3"_22.txt

sh naming-multiplex.sh

# Parse barcode
# 23

zgrep ^GAGCACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_23.fastq
zgrep ^GAGCACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_23.fastq
cat "$NAME1"_23.fastq "$NAME2"_23.fastq > "$NAME3"_23.fastq

grep AGTTC "$NAME3"_23.fastq > "$NAME3"_23.clean.fastq
mv "$NAME3"_23.clean.fastq "$NAME3"_23.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_23.fastq > "$NAME3"_23-clean.txt
rm "$NAME3"_23.fastq
mv "$NAME3"_23-clean.txt "$NAME3"_23.txt

sh naming-multiplex.sh

# Parse barcode
# 24

zgrep ^CAGTACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_24.fastq
zgrep ^CAGTACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_24.fastq
cat "$NAME1"_24.fastq "$NAME2"_24.fastq > "$NAME3"_24.fastq

grep AGTTC "$NAME3"_24.fastq > "$NAME3"_24.clean.fastq
mv "$NAME3"_24.clean.fastq "$NAME3"_24.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_24.fastq > "$NAME3"_24-clean.txt
rm "$NAME3"_24.fastq
mv "$NAME3"_24-clean.txt "$NAME3"_24.txt

sh naming-multiplex.sh

# Parse barcode
# 25

zgrep ^ACGTACCTTCTTCTTGACTCAA "$NAME1".fastq.gz | cat > "$NAME1"_25.fastq
zgrep ^ACGTACCTTCTTCTTGACTCAA "$NAME2".rev.fastq.gz | cat > "$NAME2"_25.fastq
cat "$NAME1"_25.fastq "$NAME2"_25.fastq > "$NAME3"_25.fastq

grep AGTTC "$NAME3"_25.fastq > "$NAME3"_25.clean.fastq
mv "$NAME3"_25.clean.fastq "$NAME3"_25.fastq

# Filter reads based on size
awk '{ if (length($0) >= 109) print }' "$NAME3"_25.fastq > "$NAME3"_25-clean.txt
rm "$NAME3"_25.fastq
mv "$NAME3"_25-clean.txt "$NAME3"_25.txt

sh naming-multiplex.sh

## Move cleaned 

mv *.fastq.* 01_Clean_fastq
