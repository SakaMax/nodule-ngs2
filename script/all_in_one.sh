#!/bin/bash
# workflow for nodule amplicon NGS data V2
# Usage:
# all_in_one.sh R1.fastq.gz R2.fastq.gz 
# Note:
# You can look at the latest result of fastp in http://<address-to-machine>:80/latest.html

# Check where this script is.
# see https://codechacha.com/ja/how-to-get-path-of-bash-script/
# and https://qiita.com/koara-local/items/04d3efd1031ea62d8db5
relative_dir=`dirname "$0"`
canocnical_dir=`readlink -f $relative_dir`
repo_dir=${canocnical_dir%/*}

# create folder
# ./folder_temp/folder_time
folder_time=`date +%Y_%m%d_%H%M`
folder_data="data"
out_folder_name="${repo_dir}/${folder_data}/${folder_time}"

# step 1
# rename reads by tag and cut tag seqences
echo "step 1"
mkdir -p "${out_folder_name}/tag_removed" && python3 -m cutadapt --no-indels --discard-untrimmed \
    -g "file:${repo_dir}/source/R1tag.fasta" -G file:"${repo_dir}/source/R2tag.fasta" -y ' {name}' \
    -o "${out_folder_name}/tag_removed/R1.fastq" -p "${out_folder_name}/tag_removed/R2.fastq" $1 $2
    
# step 2
# remove primer
echo "step 2"
mkdir -p "${out_folder_name}/primer_removed" && python3 -m cutadapt --discard-untrimmed \
    -g "file:${repo_dir}/source/R1primer.fasta" -G "file:${repo_dir}/source/R2primer.fasta" \
    -o "${out_folder_name}/primer_removed/R1.fastq" -p "${out_folder_name}/primer_removed/R2.fastq" \
    "${out_folder_name}/tag_removed/R1.fastq" "${out_folder_name}/tag_removed/R2.fastq"

# step 3
# filtering by fastp
echo "step 3"
mkdir -p "${out_folder_name}/fastp" && fastp -i "${out_folder_name}/primer_removed/R1.fastq" -I "${out_folder_name}/primer_removed/R2.fastq" \
    -3 -o "${out_folder_name}/fastp/R1.fastq" -O "${out_folder_name}/fastp/R2.fastq" \
    -h "${out_folder_name}/fastp/report.html" -q 30 -n 5 -A && cp "${out_folder_name}/fastp/report.html" /var/www/html/latest.html
