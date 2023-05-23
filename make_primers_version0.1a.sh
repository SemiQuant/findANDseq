#!/bin/bash
usage () { #echo -e to read \n
  echo "
findANDseq by SemiQuant

Usage Options
  -t|--threads = number of threads to use (default = 4)
  -n|--name = name for primer library
  -m|--miss_lib = a multifasta containing all the rRNAs that may be present in the sample (required)
  -b|--background = a multifasta containing all the background transcripts that may be present in teh sample (required)
  -f|--foreground = a multifasta containing all the transcripts of intrest present in the sample, should have rRNAs removed althoug this is also done using a basic grep (required)
  -mt|--min_tm = minimum primer Tm (default = 45)
  -ot|--opt_tm =  optimum primer Tm (default = 48)
  -mxt|--max_tm =  maximum primer Tm (default = 55)
  -ml|--min_len = minimum primer length (default = 8)
  -ol|--opt_len = optimum primer length (default = 15)
  -mxl|--max_len = maximum primer length (default = 30)
  -o|--out_dir = path to output directory
  -s|--Script_dir = path to script directory
  -r|--ret = primers to return for each gene, increase to find overlapping/degen primers (default = 10)
  -tm|--trim = trim backgroud fasta to be noly 1000bp max (default is on, set this flag to turn off)

# Required Software
  bowtie
  samtools
  primer3
  python3
    pandas
    Bio.Seq
  R
    strinR
    openprimeR
    biostrings
"
}

function show_progress {
    current="$1"
    total="$2"

    # calculate the progress in percentage 
    percent=$(bc <<< "scale=1; 100 * $current / $total" )
    # The number of done and todo characters
    done=$(bc <<< "scale=0; 40 * $percent / 100" )
    todo=$(bc <<< "scale=0; 40 - $done" )

    # build the done and todo sub-bars
    done_sub_bar=$(printf "%${done}s" | tr " " "#")
    todo_sub_bar=$(printf "%${todo}s" | tr " " "-")

    # output the bar
    echo -ne "\rProgress : [${done_sub_bar}${todo_sub_bar}] ${percent}%"

    if [ $total -eq $current ]; then
        echo -e "\n"
    fi
}


declare_globals () {
    # if same thing twice will take second one
    while [[ "$#" -gt 0 ]]
    do
        case $1 in
        -s|--Script_dir)
        Script_dir="$2"
        ;;
        -t|--threads)
        threads="$2"
        ;;
        -n|--name)
        name="$2"
        ;;
        -m|--miss_lib)
        miss_lib="$2"
        ;;
        -b|--background)
        background="$2"
        ;;
        -f|--foreground)
        foreground="$2"
        ;;
        -mt|--min_tm)
        min_tm="$2"
        ;;
        -ot|--opt_tm)
        opt_tm="$2"
        ;;
        -mxt|--max_tm)
        max_tm="$2"
        ;;
        -ml|--min_len)
        min_len="$2"
        ;;
        -ol|--opt_len)
        opt_len="$2"
        ;;
        -mxl|--max_len)
        max_len="$2"
        ;;
        -mxr|--max_rRNA)
        max_rRNA="$2"
        ;;
        -o|--out_dir)
        out_dir="$2"
        ;;
        -r|--ret)
        ret="$2"
        ;;
        -mm|--miss_match)
        miss_match="$2"
        ;;
        -tm|--trim)
        trim="N"
        ;;
        -k|--k_in)
        k_in="$2" #min = 1
        ;;
    esac
        shift
    done
}



get_primers () {
  echo "id=${id}
PRIMER_NUM_RETURN=${ret}
SEQUENCE_ID=${id}
PRIMER_TASK=generic
SEQUENCE_INCLUDED_REGION=${region}
SEQUENCE_TEMPLATE=${template}
PRIMER_MISPRIMING_LIBRARY=${miss_lib}
PRIMER_MAX_LIBRARY_MISPRIMING=${1}
PRIMER_MIN_SIZE=${min_len}
PRIMER_OPT_SIZE=${opt_len}
PRIMER_MAX_SIZE=${max_len}
PRIMER_MIN_TM=${min_tm}
PRIMER_OPT_TM=${opt_tm}
PRIMER_MAX_TM=${max_tm}
PRIMER_GC_CLAMP=0
PRIMER_MAX_GC=80
PRIMER_OPT_GC_PERCENT=50
PRIMER_MIN_GC=30
PRIMER_PICK_LEFT_PRIMER=0
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_MAX_HAIRPIN_TH=47.0
PRIMER_INTERNAL_MAX_HAIRPIN_TH=47.0
PRIMER_MAX_END_STABILITY=100.0
PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=-1
PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=-1
PRIMER_EXPLAIN_FLAG=0
PRIMER_LIBERAL_BASE=0
PRIMER_FIRST_BASE_INDEX=0
PRIMER_MAX_TEMPLATE_MISPRIMING=-1.00
PRIMER_MAX_TEMPLATE_MISPRIMING_TH=-1.00
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=-1.00
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=-1.00
PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=0
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/Users/semiquant/bioinfomatics/primer3-2.4.0/src/primer3_config/
=" > "${id}.primer.tmp"
}


if [ $# == 0 ]
then
    usage
    exit 1
fi


command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed. Aborting."; exit 1; }
command -v primer3_core >/dev/null 2>&1 || { echo >&2 "I require primer3_core but it's not installed. Aborting."; exit 1; }
command -v bowtie >/dev/null 2>&1 || { echo >&2 "I require bowtie but it's not installed. Aborting."; exit 1; }


if [[ ! -z "$miss_lib"  || ! -z "$background" || ! -z "$foreground" ]]
then
    usage
    exit 1
fi


if [[ ! -z "$out_dir" ]]
then
    cd "$out_dir"
fi


# set defults
declare_globals "$@"
min_tm=${min_tm:-48}
opt_tm=${opt_tm:-50}
max_tm=${max_tm:-80}
min_len=${min_len:-8}
opt_len=${opt_len:-13}
max_len=${max_len:-30}
max_rRNA=${max_rRNA:-10}
name="${name:-STOPseq_primers}"
threads=${threads:-4}
ret=${ret:-1000000000}
miss_match=${miss_match:-1}
for_genes="${name}_gene_region_table.tsv"
k_in=${k_in:-50}

Script_dir_tmp="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
Script_dir="${Script_dir:-$Script_dir_tmp}"


# trim background fasta
if [ -z ${trim+x} ]
then
    awk '/^>/ {print;next} {print substr($0,1,1000)}' "$background" > "background_truncated.fasta"
    background="background_truncated.fasta"
fi

# make bowtie index
if [ $(ls "${background/.fasta/}"*.ebwt | wc -l) -lt 1 ]
then 
    bowtie-build --threads $threads "$background" "${background/.fasta/}"
else
    echo "found index for $background"
fi

cat "$foreground" | tr -dc '>[:alnum:]\n\r' > "forground_noSpecial.fasta"
foreground="forground_noSpecial.fasta"

if [ $(ls "${foreground/.fasta/}"*.ebwt | wc -l) -lt 1 ]
then 
    bowtie-build --threads $threads "$foreground" "${foreground/.fasta/}"
else
    echo "found index for $foreground"
fi

if [ $(ls "${miss_lib/.fasta/}"*.ebwt | wc -l) -lt 1 ]
then 
    bowtie-build --threads $threads "$miss_lib" "${miss_lib/.fasta/}"
else
    echo "found index for $miss_lib"
fi


Rscript "${Script_dir}/make_pmr_table.R" "$foreground" "$for_genes"

# make primers
current_task=0
tasks_in_total=$(wc -l < "$for_genes")
echo "Picking initial primers" | tee "${name}.log"
{
read # skip header
while read -r id template str len #region
do
    current_task=$((current_task+1))
    region="${str},${len}"
    primer_found=0
    get_primers $max_rRNA 
    primer3_core --format_output < "${id}.primer.tmp" | grep "RIGHT_PRIMER" | awk 'BEGINFILE{printf " 0 "}{print}' > "${id}.primers.tsv"
    rm "${id}.primer.tmp"
    primer_found=$(wc -l < "${id}.primers.tsv")

    if [ $primer_found -gt 0 ]
    then
        awk -v id=$id '{printf ">%s_%s_%s_%s_%s_%s_%s_%s_%s\n%s\n",id,$3,$4,$5,$6,$7,$8,$9,$10,$11}' "${id}.primers.tsv" >> "${name}.allPs.fasta"
        rm "${id}.primers.tsv"
    else
        if [ $((len+str-50)) -gt $len ]
        then
            echo -e "\n trying with relaxed binding region for ${id}"
            region="50,$((len+str-50))"
            get_primers $max_rRNA
            primer3_core --format_output < "${id}.primer.tmp" | grep "RIGHT_PRIMER" | awk 'BEGINFILE{printf " 0 "}{print}' > "${id}.primers.tsv"
            primer_found=$(wc -l < "${id}.primers.tsv")
        fi
        if [ $primer_found -gt 0 ]
        then
            awk -v id=$id '{printf ">%s_%s_%s_%s_%s_%s_%s_%s_%s\n%s\n",id,$3,$4,$5,$6,$7,$8,$9,$10,$11}' "${id}.primers.tsv" >> "${name}.allPs.fasta"
            rm "${id}.primers.tsv"
        else
            echo -e "\n trying with relaxed rRNA binding for ${id}"
            get_primers 20
            primer3_core --format_output < "${id}.primer.tmp" | grep "RIGHT_PRIMER" | awk 'BEGINFILE{printf " 0 "}{print}' > "${id}.primers.tsv"
            rm "${id}.primer.tmp"
            primer_found=$(wc -l < "${id}.primers.tsv")
            if [ $primer_found -gt 0 ]
            then
                awk -v id=$id '{printf ">%s_%s_%s_%s_%s_%s_%s_%s_%s\n%s\n",id,$3,$4,$5,$6,$7,$8,$9,$10,$11}' "${id}.primers.tsv" >> "${name}.allPs.fasta"
                rm "${id}.primers.tsv"
            else
                echo -e "\n no primer found for ${id}"
            fi
        fi
    fi
    show_progress $current_task $tasks_in_total
done 
} < "$for_genes" | tee "${name}.log" 

echo -e "\n Selecting best primerset" | tee "${name}.log"
# Now check alignment to background genome
bowtie -S --all --nofw -k $k_in -v $miss_match --threads $threads -x "${background/.fasta/}" -f "${name}.allPs.fasta" --un "${name}.unalined.fasta" > "${name}.aligned.sam"
# moving this into the loop as otherwise the file can get massive...

# the primers are already the reverse complment, so changed --nofw to --norc
# if taking long, remove -a and change to -k 10, which sets the limit for secondary alignments at 10


# if [[ $(cat "${name}.unalined.fasta" | wc -l) == 0 ]]
rm "${name}.aligned.info.tsv" 2> /dev/null
if [[ -e "${name}.unaligned.fasta" ]]
then
    awk 'BEGIN{RS=">"}{print 0" "$1"\t"$2;}' "${name}.unalined.fasta" | tail -n+2 > "${name}.aligned.info.tsv"
    # samtools view -SF 4 "aligned.sam" | awk '{ print $1"\t"$3"\t"$10 }' | uniq | awk '{$0=gensub(/\s*\S+/,"",2)}1' | uniq -c  > "aligned.info.tsv"
fi
samtools view -SF 4 "${name}.aligned.sam" | awk '$4 > 35 && $5 == 255 { print $1"\t"$3"\t"$10 }' | sort | uniq | awk '{$0=gensub(/\s*\S+/,"",2)}1' | awk '{a[$0]++}END{for(x in a)print a[x], x}' >> "${name}.aligned.info.tsv"

rm "${name}.aligned.sam"

bowtie --all --nofw -a -v 0 --threads $threads -x "${foreground/.fasta/}" -f "${name}.allPs.fasta" > "${name}.aligned.info.ref.tsv" | tee "${name}.log"


# change this R section to be in the python script
if [[ $(cat "${name}.aligned.info.ref.tsv" | wc -l) -gt 0 ]]
then
    Rscript "${Script_dir}/remHost.R" "${name}.aligned.info.tsv" "${name}.aligned.info.ref.tsv" "$for_genes" "${name}_primers.tmp.tsv" 2> "${name}.log"
else
    cat "${name}.aligned.info.tsv" > "${name}_primers.tmp.tsv"
fi


python "${Script_dir}/pickBest.py" "${name}_primers.tmp.tsv" "${name}_primers.tsv" "$opt_tm" "${name}_noTails_primers.fasta" | tee "${name}.log"

bowtie --all --nofw -a -v 0 --threads $threads -x "${foreground/.fasta/}"  \
    -f "${name}_noTails_primers.fasta" > "${name}.alignments.tsv"

echo "genes covered in primerset ${name}.alignments.tsv" | tee "${name}.log"
awk '{print $3}' "${name}.alignments.tsv" | sort | uniq | wc -l | tee "${name}.log"

echo "total bases with tailed primerset" | tee "${name}.log"
wc -c "${name}_primers.tsv" | tee "${name}.log"

echo "lengths of primers in set" | tee "${name}.log"
grep -v ">" "${name}_noTails_primers.fasta" | awk '{print length}' | sort -n | uniq -c | tee "${name}.log"

awk '{print $3"\t"$4}' "${name}.alignments.tsv" > "${name}.alignments.positions.tsv"
echo "info on lengths of products" | tee "${name}.log"
range=$(awk 'NR==1{min=max=$2;next}{if($2>max)max=$2;if($2<min)min=$2}END{print min, max}' "${name}.alignments.positions.tsv")
echo "Range: $range" | tee "${name}.log"
mean=$(awk '{ sum += $2 } END { if (NR > 0) print sum / NR }' "${name}.alignments.positions.tsv")
echo "Mean: $mean" | tee "${name}.log"
median=$(awk '{ a[i++]=$2; } END { n=i; if (n%2==0) print (a[n/2-1]+a[n/2])/2; else print a[int(n/2)]; }' "${name}.alignments.positions.tsv")
echo "Median: $median" | tee "${name}.log"
IQR=$(awk 'BEGIN{FS="\t"}{a[NR]=$2}END{n=asort(a);q1=(a[int(n/4)]+a[int(n/4)+1])/2;q3=(a[int(n*3/4)]+a[int(n*3/4)+1])/2;print q1", " q3}' "${name}.alignments.positions.tsv")
echo "IQR: $IQR" | tee "${name}.log"


echo "off target alignments to rRNAs with upto $miss_match mismatches" | tee "${name}.log"
bowtie --all --nofw -a -v $miss_match --threads $threads -x "${miss_lib/.fasta/}"  \
    -f "${name}_noTails_primers.fasta" | awk '{print $3"\t"$4}' > "${name}.miss_lib.alignments.positions.tsv"
echo "info on lengths of products" | tee "${name}.log"
range=$(awk 'NR==1{min=max=$2;next}{if($2>max)max=$2;if($2<min)min=$2}END{print min, max}' "${name}.miss_lib.alignments.positions.tsv")
echo "Range: $range" | tee "${name}.log"
mean=$(awk '{ sum += $2 } END { if (NR > 0) print sum / NR }' "${name}.miss_lib.alignments.positions.tsv")
echo "Mean: $mean" | tee "${name}.log"
median=$(awk '{ a[i++]=$2; } END { n=i; if (n%2==0) print (a[n/2-1]+a[n/2])/2; else print a[int(n/2)]; }' "${name}.miss_lib.alignments.positions.tsv")
echo "Median: $median" | tee "${name}.log"
IQR=$(awk 'BEGIN{FS="\t"}{a[NR]=$2}END{n=asort(a);q1=(a[int(n/4)]+a[int(n/4)+1])/2;q3=(a[int(n*3/4)]+a[int(n*3/4)+1])/2;print q1", " q3}' "${name}.miss_lib.alignments.positions.tsv")
echo "IQR: $IQR" | tee "${name}.log"



echo "off target alignments to background genome with upto $miss_match mismatches" | tee "${name}.log"
bowtie --all --nofw -a -v $miss_match --threads $threads -x "${background/.fasta/}"  \
    -f "${name}_noTails_primers.fasta" | awk '{print $3"\t"$4}' > "${name}.background.alignments.positions.tsv"
echo "numer of unique genes" | tee "${name}.log"
awk '{ count[$1]++ } END { print length(count) }' "${name}.background.alignments.positions.tsv"
echo "info on lengths of products" | tee "${name}.log"
range=$(awk 'NR==1{min=max=$2;next}{if($2>max)max=$2;if($2<min)min=$2}END{print min, max}' "${name}.background.alignments.positions.tsv")
echo "Range: $range" | tee "${name}.log"
mean=$(awk '{ sum += $2 } END { if (NR > 0) print sum / NR }' "${name}.background.alignments.positions.tsv")
echo "Mean: $mean"
median=$(awk '{ a[i++]=$2; } END { n=i; if (n%2==0) print (a[n/2-1]+a[n/2])/2; else print a[int(n/2)]; }' "${name}.background.alignments.positions.tsv")
echo "Median: $median" | tee "${name}.log"
IQR=$(awk 'BEGIN{FS="\t"}{a[NR]=$2}END{n=asort(a);q1=(a[int(n/4)]+a[int(n/4)+1])/2;q3=(a[int(n*3/4)]+a[int(n*3/4)+1])/2;print q1", " q3}' "${name}.background.alignments.positions.tsv")
echo "IQR: $IQR" | tee "${name}.log"

rm "${name}_primers.tmp.tsv" "forground_noSpecial.1.ebwt" "forground_noSpecial.2.ebwt" "forground_noSpecial.3.ebwt" "forground_noSpecial.4.ebwt" "forground_noSpecial.fasta" "forground_noSpecial.rev.1.ebwt" "forground_noSpecial.rev.2.ebwt"

if [ -z ${trim+x} ]
then
  rm "background_truncated.fasta"
fi
