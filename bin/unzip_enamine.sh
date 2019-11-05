

#!/bin/bash

usage()
{   
    echo "Unzip Enamine database to smiles"
    echo "./unzip_enamine.sh --cpus 4 | [-h]]"
    echo "--cpus: 1, 2, 3, 4, 6 or 12"
}

while [ "$1" != "" ]; do
    case $1 in
        -c | --cpus )           shift
                                cpus=$1
                                ;;
        -q | --inputfolder)     shift
                                folder=$1
                                ;;
        # -f | --force)           shift
        #                         overwrite=$1
        #                         ;;
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
folder=${folder:-data}
# overwrite=${overwrite:-false}

echo "Use $cpus cpus"
echo "Use input-folder: $inputfolder"
cd $folder

# echo "Overwrite previous results: $overwrite"

if ((12 % $cpus)); then
    echo "Pleas select for --cups either 1, 2, 3, 4, 6 or 12"
    exit 1
fi

# unzip -n '*.zip'
find . -name '*.zip' -print0 | xargs -0 -I {} -P $cpus  unzip -n {}
# https://askubuntu.com/questions/431478/decompressing-multiple-files-at-once
