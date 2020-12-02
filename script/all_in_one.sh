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

# create folders
# ./folder_temp/folder_time
folder_time=`date +%Y_%m%d_%H%M`
folder_data="data"
out_folder_name="${repo_dir}/${folder_data}/${folder_time}"
mkdir -p $out_folder_name

# create the log file
all_in_one_log="${out_folder_name}/all_in_one.log"
touch $all_in_one_log

echo "All-in-One script for Nodule NGS!" |& tee -a $all_in_one_log
echo "Run started in ${folder_time}" |& tee -a $all_in_one_log

# step 1
# rename reads by tag and cut tag seqences
echo  "----step 1: cutadapt (pass 1/2)" |& tee -a $all_in_one_log
mkdir -p "${out_folder_name}/tag_removed" 
python3 -m cutadapt --no-indels --discard-untrimmed \
    -g "file:${repo_dir}/source/R1tag.fasta" -G file:"${repo_dir}/source/R2tag.fasta" -y ' {name}' \
    -o "${out_folder_name}/tag_removed/R1.fastq" -p "${out_folder_name}/tag_removed/R2.fastq" $1 $2 \
    |& tee -a $all_in_one_log
    
# step 2
# remove primer
echo "----step 2: cutadapt (pass 2/2)" |& tee -a $all_in_one_log
mkdir -p "${out_folder_name}/primer_removed"
python3 -m cutadapt --discard-untrimmed \
    -g "file:${repo_dir}/source/R1primer.fasta" -G "file:${repo_dir}/source/R2primer.fasta" \
    -o "${out_folder_name}/primer_removed/R1.fastq" -p "${out_folder_name}/primer_removed/R2.fastq" \
    "${out_folder_name}/tag_removed/R1.fastq" "${out_folder_name}/tag_removed/R2.fastq" \
    |& tee -a $all_in_one_log

# step 3
# filtering by fastp, then copy report to html folder
echo "----step 3: fastp" |& tee -a $all_in_one_log
mkdir -p "${out_folder_name}/fastp"
fastp -i "${out_folder_name}/primer_removed/R1.fastq" -I "${out_folder_name}/primer_removed/R2.fastq" \
    -3 -o "${out_folder_name}/fastp/R1.fastq" -O "${out_folder_name}/fastp/R2.fastq" \
    -h "${out_folder_name}/fastp/report.html" -q 30 -n 5 -A \
    |& tee -a $all_in_one_log
cp "${out_folder_name}/fastp/report.html" /var/www/html/latest.html
echo "Now you can find the report at <docker machine's IP>:8080/latest.html" |& tee -a $all_in_one_log

# step 4
# demultiplex(python)
echo "----step 4: demultiplex" |& tee -a $all_in_one_log
echo "This step may take a while." |& tee -a $all_in_one_log
python ./script/demultiplex.py \
    "${out_folder_name}/fastp/R1.fastq" "${out_folder_name}/fastp/R2.fastq" "${out_folder_name}" \
    |& tee -a $all_in_one_log