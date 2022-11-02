#!/bin/bash
usage(){
   echo ""
   echo "Usage: $0 -p <project directory>"
   echo -e "\t-p specifies the name of your project directory"
   exit 1 # Exit script after printing help
}

while getopts ":p:" opt; do
    case $opt in
        p) 
            projectDir="$OPTARG"
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

baseDir='/stornext/General/data/user_managed/grpu_jchoi_0/projects/davide'

echo
echo "Setting env vars from the config file located in the ${baseDir}/${projectDir} project dir"
echo "-------------------------------------"
eval "$(python3 ${baseDir}/${projectDir}/utils/parse_configyaml.py -pd $projectDir)"
# echo
# env | egrep "project=|account=|nstasks=|threads=|cpus_per_task=|time=|partition=|project_dir=|modules_configfile="
# echo
# echo
echo 
echo "Redirecting to project directory: $projectDir"
echo 
cd ${baseDir}/${projectDir}
echo
cmd="salloc --ntasks=$ntasks --threads=$threads --cpus-per-task=$cpus_per_task \
--mem=$mem --time=$time --partition=$partition"
echo
echo "Launching interactive session"
echo
echo '!!!!After this step remember to load all contained in the modules.txt file as eval "$(cat "$modules")" !!!!'
echo $cmd
echo
$cmd