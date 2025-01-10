#! /bin/bash


sum_otu_table_reads() {
    table_file=$1
    with_taxonomy=$2
    if [ -z "$table_file" ]; then
        echo "Usage: sum_otu_table_reads <table_file>"
        exit 1
    fi
    # if with_taxonomy is true, then sum the reads in the last column
    if [ "$with_taxonomy" = true ]; then
        awk -F'\t' 'BEGIN {sum=0} NR>1 {for(i=2;i<NF;i++) sum+=$i} END {print FILENAME " Reads Count: " sum}' $table_file
    else
        awk -F'\t' 'BEGIN {sum=0} NR>1 {for(i=2;i<NF+1;i++) sum+=$i} END {print FILENAME " Reads Count: " sum}' $table_file
    fi
}

sum_fasta_reads() {
    fasta_file=$1
    if [ -z "$fasta_file" ]; then
        echo "Usage: sum_fasta_reads <fasta_file>"
        exit 1
    fi
    echo -e $fasta_file "Reads Count: " $(egrep -E "size=[0-9]+[\s;]" $fasta_file | awk -F'size=' '{print $2}' | awk -F'[ ;]' 'BEGIN {sum=0} {sum+=$1} END {print sum}')

}


main() {
    # Main script goes here
    help_message="Usage: sum_reads.sh [-t <table_file> [-wt|--with_taxonomy] | -f|--fasta_file <fasta_file>]"
    if [ "$#" -lt 1 ]; then
        echo "$help_message"
        exit 1
    fi
    with_taxonomy=false
    # Parse command line arguments
    while [ "$#" -gt 0 ]; do
        case "$1" in
            -t)
                table_file="$2"
                shift 2
                ;;
            -wt|--with_taxonomy)
                with_taxonomy=true
                shift
                ;;
            -f|--fasta_file)
                fasta_file="$2"
                shift 2
                ;;
            *)
                echo "$help_message"
                exit 1
                ;;
        esac
    done

    if [ -n "$table_file" ]; then
        sum_otu_table_reads $table_file $with_taxonomy
    fi
    if [ -n "$fasta_file" ]; then
        sum_fasta_reads $fasta_file
    fi
    if [ -z "$table_file" ] && [ -z "$fasta_file" ]; then
        echo "Usage: sum_reads.sh [-t <table_file> [-wt|--with_taxonomy] | -f <fasta_file>]"
        exit 1
    fi

}

main "$@"