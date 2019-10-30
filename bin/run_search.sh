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
        -e | --env )            shift
                                env=$1
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
echo "Use $env and $cpus cpus"
env=${env:-rdkit}
cpus=${cpus:-2}
#conda activate $env
echo "Use $env and $cpus cpus"

if ((12 % $cpus)); then
    echo "Pleas select for --cups either 1, 2, 3, 4, 6 or 12"
    exit 1
fi

n_iter=$((12/$cpus))
echo "Perform $n_iter loops."

cmd=''
for NO in 0{1..9} {10..12}
do
    cmd="$cmd python search_database_singleprocess.py --input_folder ./blobs/part_$NO --pattern pkl &"
    if ! (($NO % $cpus)); then
	    echo $cmd wait
        eval $cmd wait
        cmd=''
    fi
done
sleep 30