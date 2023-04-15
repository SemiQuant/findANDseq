# slective transcriptomics of pathogens
# STOPseq

# doesnt really check if it found all primers
# Would be nice to weight the backgroud genome, kind of like I am with the rRNA


usage () { #echo -e to read \n
  echo "
Usage Options
  -t|--threads = number of threads to use (default = 4)
  -n|--name = name for primer library
  -m|--miss_lib = a multifasta containing all the rRNAs that may be present in the sample (required)
  -b|--background = a multifasta containing all the background transcripts that may be present in teh sample (required)
  -f|--foreground = a multifasta containing all the transcripts of intrest present in the sample, should have rRNAs removed althoug this is also done using a basic grep (required)
  -mt|--min_tm = minimum primer Tm (default = 50)
  -ot|--opt_tm =  optimum primer Tm (default = 50)
  -mxt|--max_tm =  maximum primer Tm (default = 80)
  -ml|--min_len = minimum primer length (default = 8)
  -ol|--opt_len = optimum primer length (default = 13)
  -mxl|--max_len = maximum primer length (default = 30)
  -o|--out_dir = path to output directory
  -s|--Script_dir = path to script directory

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
    esac
        shift
    done
}



get_primers () {
  echo "id=${id}
PRIMER_NUM_RETURN=10000
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
min_tm=${min_tm:-50}
opt_tm=${opt_tm:-50}
max_tm=${max_tm:-80}
min_len=${min_len:-8}
opt_len=${opt_len:-13}
max_len=${max_len:-30}
max_rRNA=${max_rRNA:-0}
name="${name:-STOPseq_primers}"
threads=${threads:-4}
for_genes="${name}_gene_region_table.tsv"

Script_dir_tmp="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
Script_dir="${Script_dir:-$Script_dir_tmp}"



# # should have all the rRNAs, like also the human ones
# miss_lib="/Users/semiquant/Bioinformatics/Projects/selectivecDNA/rRNAs/rRNAs.fasta"
# background="/Users/semiquant/Bioinformatics/Projects/selectivecDNA/gencode.v38.transcripts_rRNAatend.fasta"
# foreground="/Users/semiquant/Bioinformatics/Projects/selectivecDNA/Mycobacterium_tuberculosis_H37Rv_genes_v4.0.fasta"
# for_genes="/Users/semiquant/Bioinformatics/Projects/selectivecDNA/gene_region_table.tsv"

# make bowtie index
if [ $(ls "${foreground/.fasta/}"*.ebwt | wc -l) -lt 1 ]
then 
    bowtie-build --threads $threads "$foreground" "${foreground/.fasta/}"
else
    echo "found index for $foreground"
fi

if [ $(ls "${background/.fasta/}"*.ebwt | wc -l) -lt 1 ]
then 
    bowtie-build --threads $threads "$background" "${background/.fasta/}"
else
    echo "found index for $background"
fi

cat "$foreground" | tr -dc '>[:alnum:]\n\r' > "forground_noSpecial.fasta"
foreground="forground_noSpecial.fasta"

Rscript "${Script_dir}/make_pmr_table.R" "$foreground" "$for_genes"


# make primers
while read -r id template str len #region
do
    region="${str},${len}"
    primer_found=0
    if [ $primer_found -lt 1 ]
    then
        get_primers $max_rRNA 
        primer3_core --format_output < "${id}.primer.tmp" | grep "RIGHT_PRIMER" | awk 'BEGINFILE{printf " 0 "}{print}' > "${id}.primers.tsv"

        rm "${id}.primer.tmp"
        primer_found=$(wc -l < "${id}.primers.tsv")
        if [ $primer_found -gt 0 ]
        then
            awk -v id=$id '{printf ">%s_%s_%s_%s_%s_%s_%s_%s_%s\n%s\n",id,$3,$4,$5,$6,$7,$8,$9,$10,$11}' "${id}.primers.tsv" >> "${name}.allPs.fasta"
            rm "${id}.primers.tsv"
        else
            echo "trying with relaxed binding region for ${id}"
            region="50,$((len+str-50))"
            get_primers $max_rRNA
            if [ $primer_found -gt 0 ]
            then
                    echo "trying with relaxed rRNA binding for ${id}"
                    get_primers 10
                    primer3_core --format_output < "${id}.primer.tmp" | grep "RIGHT_PRIMER" | awk 'BEGINFILE{printf " 0 "}{print}' > "${id}.primers.tsv"
                    rm "${id}.primer.tmp"
                    primer_found=$(wc -l < "${id}.primers.tsv")
                    if [ $primer_found -gt 0 ]
                    then
                        awk -v id=$id '{printf ">%s_%s_%s_%s_%s_%s_%s_%s_%s\n%s\n",id,$3,$4,$5,$6,$7,$8,$9,$10,$11}' "${id}.primers.tsv" >> "${name}.allPs.fasta"
                        rm "${id}.primers.tsv"
                    else
                        echo "no primer found for ${id}"
                    fi
            fi
        fi
    fi
done < "$for_genes"


# Now check alignment to human genome
bowtie -S --all --nofw -a -v 2 --threads $threads -x "${background/.fasta/}" -f "${name}.allPs.fasta" --un "${name}.unalined.fasta" > "${name}.aligned.sam"
# the primers are already the reverse complment, so changed --nofw to --norc
# if taking long, remove -a and change to -k 10, which sets the limit for secondary alignments at 10


if [[ $(cat "${name}.unalined.fasta" | wc -l) == 0 ]]
then
    # samtools view -SF 4 "aligned.sam" | awk '{ print $1"\t"$3"\t"$10 }' | uniq | awk '{$0=gensub(/\s*\S+/,"",2)}1' | uniq -c  > "aligned.info.tsv"
    samtools view -SF 4 "${name}.aligned.sam" | awk '$4 > 35 && $5 == 255 { print $1"\t"$3"\t"$10 }' | sort | uniq | awk '{$0=gensub(/\s*\S+/,"",2)}1' | awk '{a[$0]++}END{for(x in a)print a[x], x}' > "${name}.aligned.info.tsv"
else
    awk 'BEGIN{RS=">"}{print 0" "$1"\t"$2;}' "${name}.unalined.fasta" | tail -n+2 > "${name}.aligned.info.tsv"
fi

rm "${name}.aligned.sam"

bowtie --all --nofw -a -v 0 --threads $threads -x "${foreground/.fasta/}" -f "${name}.allPs.fasta" > "${name}.aligned.info.ref.tsv"


# change this R section to be in the python script
Rscript "${Script_dir}/remHost.R" "${name}.aligned.info.tsv" "${name}.aligned.info.ref.tsv" "$for_genes" "${name}_primers.tmp.tsv"

python "${Script_dir}/pickBest.py" "${name}_primers.tmp.tsv" "${name}_primers.tsv" 

rm "${name}_primers.tmp.tsv" "forground_noSpecial.fasta"

# bowtie --all --nofw -a -v 2 --threads 6 -x "/Users/semiquant/Bioinformatics/Projects/selectivecDNA/rRNAs/rRNAs" \
#     -f "/Users/semiquant/Downloads/del/primer_notails_v7.fasta" --un "test.unalined.fasta" > "test.rRNA.aligned.txt"



# bowtie --all --nofw -a -v 0 --threads 6 -x "/Users/semiquant/Bioinformatics/Projects/selectivecDNA/Mycobacterium_tuberculosis_H37Rv_genes_v4.0" \
#     -f "/Users/semiquant/Downloads/del/primer_notails_v7.fasta" --un "test.unalined.fasta" > "test.forgrounnd.aligned.txt"


# GRAVEYARD
# or both together
# samtools view -SF 4 "${final_fasta/.fasta/aligned.sam}" | awk '{ print $1"\t"$3"\t"$10 }' | sort | uniq | awk '{$0=gensub(/\s*\S+/,"",2)}1' |  awk '{a[$0]++}END{for(x in a)print a[x], x}' > "${final_fasta/.fasta/.aligned.info.tsv}"
# awk 'BEGIN{RS=">"}{print 0" "$1"\t"$2;}' "${final_fasta/.fasta/unalined.fasta}" | tail -n+2 >> "${final_fasta/.fasta/.aligned.info.tsv}"

# neet to store which genes..so can use that info later
#samtools view -SF 4 "${name}.aligned_ref.sam" | awk '$4 > 50 && $5 == 255 { print $1"\t"$3"\t"$10 }' | sort | uniq > "${name}.aligned.info.ref.tsv"
#rm "${name}.aligned_ref.sam"

# bowtie2 -a --nofw --ma 7 --score-min L,0,7 \
#     --local -p 6 -x "/Users/semiquant/Bioinformatics/Projects/selectivecDNA/Mycobacterium_tuberculosis_H37Rv_genes_v4.0" \
#     -f "primers_out_14_Mar_2023_noTails.fasta" > "final.aligned_ref.sam"

# and now to the reference genome, there will be at least one per seq
# perfect matches only
# --ma 7 --score-min L,0,7
# bowtie2 -a --nofw --ma 7 --score-min L,0,7 --local -p $threads -x "/Users/semiquant/Bioinformatics/Projects/selectivecDNA/Mycobacterium_tuberculosis_H37Rv_genes_v4.0" -f "${name}.allPs.fasta" > "${name}.aligned_ref.sam"
# bowtie2 -p $threads -a --nofw --end-to-end --score-min L,0,0 -x "/Users/semiquant/Bioinformatics/Projects/selectivecDNA/Mycobacterium_tuberculosis_H37Rv_genes_v4.0" \
    # -f "${name}.allPs.fasta" > "${name}.aligned_ref.sam"

# samtools view -bS "final.aligned_ref.sam" | samtools sort -o "${name}.aligned_ref.bam"
# samtools index "final.aligned_ref.bam"
# samtools view /Users/semiquant/Bioinformatics/Projects/selectivecDNA/result/using_primer3/final.aligned_ref.bam | awk '{print $1"\t"$3}' > alignments.tsv


# in R
# require(tidyverse)
# dat <- read_tsv("/Users/semiquant/Bioinformatics/Projects/selectivecDNA/result/using_primer3/alignments.tsv", col_names = F)
# dat$X2 <- gsub("\\|.*", "", dat$X2)
# length(unique(dat$X2))
# genes <- read_tsv("/Users/semiquant/Bioinformatics/Projects/selectivecDNA/gene_region_table.tsv")
# a=genes$ID[!genes$ID %in% unique(dat$X2)]





# make a not tail version if you want to check here
# bowtie --all --nofw -a -v 1 --threads $threads -x "/Users/semiquant/Bioinformatics/Projects/selectivecDNA/Mycobacterium_tuberculosis_H37Rv_genes_v4.0" \
#     -f "primers_out_14_Mar_2023_noTails.fasta" > "${name}.alignments.tsv"

# echo "genes covered in primerset ${name}.alignments.tsv"
# awk '{print $3}' "${name}.alignments.tsv" | sort | uniq | wc -l












