#!/bin/bash
#author: Yuan Gao
#contact: gy.james@163.com
source ~/.bashrc

usage() { echo "usage: $0 <-I> <-O> <-D> <-L> <-S> [-T] " 1>&2; exit 1; }

[ $# -eq 0 ] && usage
while getopts ":hI:O:T:L:D:S:" arg; do
  	case $arg in
	    I) # Specify input
			I=${OPTARG}
			echo "input of fastq is ${I}"
			;;
	    O) # Specify output
			O=${OPTARG}
			echo "output directory is ${O}"
			;;
		L) # Specify length of UMI
			#L=${OPTARG}
			#((L > 30)) || usage
			#;;
			L=${OPTARG}
			if [ $L -lt 0 ]; then
				echo "length of UMI cannot be less than 0."
				usage 
			else
				echo "length of UMI is $L."
			fi
			;;
		D) # Specify database of splint sequence
			D=${OPTARG}
			echo "database of splint is ${D}"
			;;
		S) # Specify source code path
			S=${OPTARG}
			echo "source code path is ${S}"
			;;
	    T) # Specify thread
			T=${OPTARG}
			#((T > 0)) || usage
			#;;
			if [ $T -gt 0 ]; then
			echo "thread is $T."
			else
				echo "thread needs to be larger than 0."
				usage 
			fi
			;;
	    H | *) # Display help.
			usage
			;;
	esac
done

if [ -z $T ]; then
	T=1
	echo "thread is 1"
fi

if [ -z ${I+x} ]; then
	echo "input needs to designated"
	usage
fi

if [ ! -d "$I" ]; then
	echo "input cannot be found"
	usage
fi

if [ -z ${O+x} ]; then
	echo "output needs to designated"
	usage
fi

if [ ! -d "$O" ]; then
	mkdir $O
	if [ $? -eq 0 ]; then
	    echo "output ${O} is made"
	else
	    echo "output ${O} cannot be made"
	    exit 1;
	fi
fi

HMMER_CHECK=$(nhmmer -h | awk '$2=="HMMER"&&$3=="3.3.1"{print "nhmmer OK"}') #{print $3}' | awk -F '.' '$1==3&&$2==3&&$3==1
if [ "$HMMER_CHECK" = "nhmmer OK" ];
then
	echo "nhmmer OK"
else
	echo "nhmmer not OK"
	exit 2
fi

ABPOA_CHECK=$(abpoa -v)
if [ "$ABPOA_CHECK" = "1.0.6" ];
then
	echo "abpoa OK"
else
	echo "abpoa not OK"
	exit 2
fi

LINK_LN=$(awk '{print length($0)}' $D | tail -1)
echo "length of splint is ${LINK_LN}"

NAME='out'

date +%T

if sed -n '1~4s/^@/>/p;2~4p' $I/*.fastq > $O/$NAME.fa; then
    echo "fastq converted to fasta"
else
    echo "fastq cannot be converted to fasta"
    exit 1;
fi
date +%T

if nhmmer -T 8 --max --cpu $T --tblout $O/${NAME}_T8_max.hmmer1 -o $O/${NAME}_T8_max.hmmer --qformat fasta --qsingle_seqs $D $O/$NAME.fa; then
    echo "nhmmer for raw reads done"
else
    echo "nhmmer for raw reads error"
    exit 1;
fi
date +%T

if sort $O/${NAME}_T8_max.hmmer1 > $O/${NAME}_T8_max.hmmer1_sort; then
    echo "sorting nhmmer done"
else
    echo "sorting nhmme error"
    exit 1;
fi

if perl $S/repeat_extract.pl $O/$NAME.fa $O/${NAME}_T8_max.hmmer1_sort $LINK_LN > $O/extract_unit_splints_RAW_versatile_$NAME.fa0; then
    echo "repeat_extract done"
else
    echo "repeat_extract error"
    exit 1;
fi
date +%T
rm $O/extract_unit_splints_consensus_versatile_$NAME.fa

if perl $S/parallel_abpoa.pl $O/extract_unit_splints_RAW_versatile_$NAME.fa0 $O/extract_unit_splints_consensus_versatile_$NAME.fa $T; then
    echo "parallel_abpoa done"
else
    echo "parallel_abpoa error"
    exit 1;
fi
date +%T

if nhmmer -T 8 --max --cpu $T --tblout $O/extract_unit_splints_consensus_versatile_$NAME.hmmer1 -o $O/extract_unit_splints_consensus_versatile_$NAME.hmmer --qformat fasta --qsingle_seqs $D $O/extract_unit_splints_consensus_versatile_$NAME.fa; then
    echo "nhmmer for consensus done"
else
    echo "nhmmer for consensus error"
    exit 1;
fi
date +%T

if [ $L -gt 0 ]; then
	if perl $S/UMI_extract.pl $O/extract_unit_splints_consensus_versatile_$NAME.fa $O/extract_unit_splints_consensus_versatile_$NAME.hmmer1 $LINK_LN $L 0 GGGAA ACTAT > $O/extract_UMI_cDNA_versatile_hmmer_$NAME.txt; then
	    echo "UMI_extract done"
	else
	    echo "UMI_extract error"
	    exit 1;
	fi
	date +%T
fi
#cd /mnt/isilon/xing_lab/aspera/Feng/test_LRCA_pipeline2/
#qsub -l h_vmem=8G -l m_mem_free=8G -pe smp 20 src2/spirit.sh -I TEST -O OUT -D /mnt/isilon/xing_lab/aspera/Feng/analysis/completeness_evaluation/splint_50nt.fa -L 60 -S src2 -T 10