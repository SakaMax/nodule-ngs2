#!/bin/bash
# workflow for nodule amplicon NGS data V2
# Usage:
# all_in_one.sh R1.fastq.gz R2.fastq.gz 
# Note:
# You can look at the latest result of fastp in http://<address-to-machine>:80/latest.html


# create folder
# ./folder_temp/folder_time
folder_time=`date +%Y_%m%d_%H%M`
folder_temp="tmp_file"
folder_name="${folder_temp}/${folder_time}"

# step 1
# rename reads by tag and cut tag seqences
echo "step 1"
mkdir -p "${folder_name}/tag_removed" && python3 -m cutadapt --no-indels --discard-untrimmed \
    -g file:./source/R1tag.fasta -G file:./source/R2tag.fasta -y ' {name}' \
    -o "${folder_name}/tag_removed/R1.fastq" -p "${folder_name}/tag_removed/R2.fastq" $1 $2
    
# step 2
# remove primer
echo "step 2"
mkdir -p "${folder_name}/primer_removed" && python3 -m cutadapt --discard-untrimmed \
    -g file:./source/R1primer.fasta -G file:./source/R2primer.fasta \
    -o "${folder_name}/primer_removed/R1.fastq" -p "${folder_name}/primer_removed/R2.fastq" \
    "${folder_name}/tag_removed/R1.fastq" "${folder_name}/tag_removed/R2.fastq"

# step 3
# filtering by fastp
echo "step 3"
mkdir -p "${folder_name}/fastp" && fastp -i "${folder_name}/primer_removed/R1.fastq" -I "${folder_name}/primer_removed/R2.fastq" \
    -3 -o "${folder_name}/fastp/R1.fastq" -O "${folder_name}/fastp/R2.fastq" \
    -h "${folder_name}/fastp/report.html" -q 30 -n 5 -A && cp "${folder_name}/fastp/report.html" /var/www/html/latest.html
