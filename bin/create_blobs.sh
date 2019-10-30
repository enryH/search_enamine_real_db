#!/bin/bash

usage()
{   
    echo "Creates binary files for Enamine project."
    echo "create_blobs -e rdkit-env --cpus 4 | [-h]]"
    echo "-n: 2, 3, 4, 6 or 12"
    echo "Example:"
    # echo "qsub_executor.sh --script alignment.sh --numbers \"13 14 15 16\""
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

echo "Use $env and $cpus cpus"
env=${env:-rdkit}
cpus=${cpus:-2}
#conda activate $env
echo "Use $env and $cpus cpus"

if ((12 % $cpus)); then
    echo "Pleas select for --cups either 2, 3, 4, 6 or 12"
    exit 1
fi

n_iter=$((12/$cpus))
echo "Perform $n_iter loops."

#python create_blobs.py -n 1 & python create_blobs.py -n 2 & python create_blobs.py -n 3 & python create_blobs.py -n 4 & python create_blobs.py -n 5 & python create_blobs.py -n 6 & python create_blobs.py -n 7 & python create_blobs.py -n 8 & python create_blobs.py -n 9 & python create_blobs.py -n 10 & python create_blobs.py -n 11 & python create_blobs.py -n 12 & wait

cmd=''
for NO in {1..12..1}
do
    cmd="$cmd python create_blobs.py -n $NO &"
    if ! (($NO % $cpus)); then
	echo "$cmd wait"
        #bash #!/bin/bash; $cmd wait
	eval $cmd wait
        cmd=''
    fi
done
