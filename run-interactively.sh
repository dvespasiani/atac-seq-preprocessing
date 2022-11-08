#!/usr/bin/env bash

## some default variables
cf='/config/project_config.yaml'
basedir=$(pwd)

usage(){
    echo ""
    echo "Usage: $0 -b <base dir> -p <project name> -cf <config yaml>"
    echo -e "\t-b specifies the path of your base directory. This SHOULD NOT BE your $HOME" 
    echo -e "\t-p specifies the name of your project directory, this is relative to your base dir"
    echo -e "\t-cf specifies the name of your config yaml containing the salloc directives which, by default, is $cf"
    exit 1 # Exit script after printing help
}


while getopts ":b:p:cf:" opt; do
    case $opt in
        b) 
            basedir="$OPTARG"
            ;;
        p) 
            projectname="$OPTARG"
            ;;
        cf) 
            config_yaml="$OPTARG"
            ;;
        *)
            echo "invalid command $OPTARG"
            ;;
        ?) 
            usage 
            ;; # Print helpFunction
    esac
done
shift $((OPTIND -1))

projectdir="$basedir/$projectname"

echo "changing working directory to $projectdir"
cd $projectdir

echo
echo "Parsing variables from $cf file"
echo "-------------------------------------"
eval "$(python3 ${projectdir}/utils/parse-yaml.py -y $cf)"
echo
cmd="salloc --ntasks=$ntasks --threads=$threads --cpus-per-task=$cpus_per_task \
--mem=$mem --time=$time --partition=$partition"
echo
echo "Launching interactive session"
echo
echo '!!!!After this step remember to load all contained in the modules.txt file as eval "$(cat "$modules")" !!!!'
echo "-------------------------------------"
echo $cmd
echo
$cmd


