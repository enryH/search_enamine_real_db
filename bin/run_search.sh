#!/bin/bash

usage()
{   
    echo "Search Enamine database using stored fingerprints in ./blobs/"
    echo "./run_search.sh -e rdkit-env --cpus 4 | [-h]]"
    echo "--cpus: 1, 2, 3, 4, 6 or 12"
}

while [ "$1" != "" ]; do
    case $1 in
        -c | --cpus )           shift
                                cpus=$1
                                ;;
        -q | --query)           shift
                                query=$1
                                ;;
        -f | --force)           shift
                                overwrite=$1
                                ;;
        -t | --threshold)       shift
                                threshold=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

echo "Started script in folder: $PWD"
cpus=${cpus:-2}
threshold=${threshold:-0.7}
overwrite=${overwrite:-false}

echo "Use $cpus cpus"
echo "Threshold: $threshold"
echo "Overwrite previous results: $overwrite"

if ((12 % $cpus)); then
    echo "Pleas select for --cups either 1, 2, 3, 4, 6 or 12"
    exit 1
fi

n_iter=$((12/$cpus))
echo "Perform $n_iter loops."

cmd=''
for NO in {1..12} #  0{1..9} {10..12};
do
    if (($NO < 10)); then
	    cmd="$cmd python search_database_singleprocess.py --input_folder ./blobs/part_0$NO"
    else
	    cmd="$cmd python search_database_singleprocess.py --input_folder ./blobs/part_$NO"
    fi
    cmd="$cmd --pattern pkl --reference_mol '$query' --tanimoto_threshold '$threshold' --force $overwrite & "
    if ! (($NO % $cpus)); then
	echo "$cmd wait"
        eval "$cmd wait"
        cmd=''
    fi
done
eval python combine_results.py --resultpath 'results/$query' --threshold $threshold
