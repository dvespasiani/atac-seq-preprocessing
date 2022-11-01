import argparse
import os
import yaml

homeDir = os.path.expanduser("/stornext/General/data/user_managed/grpu_jchoi_0/projects/davide/") # your home directory

def config_arguments():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-bd', '--baseDir', type=str, default=homeDir, help='path to base working directory, default is the home directory')
    parser.add_argument('-pd', '--projectDir', type=str, default='atac-pipeline/', help='path to project dir (relative to base dir)')
    
    args = parser.parse_args()  
    args.interSessYAML = os.path.join(args.baseDir,args.projectDir, "config/cluster_config.yaml")
    args.projectYAML = os.path.join(args.baseDir,args.projectDir, "config/project_config.yaml") 
    return args


config = config_arguments()

def parseYAML(y):
    with open(y, 'r') as file:
        yamlFile = yaml.safe_load(file)  
    for k, v in yamlFile.items():    
        print('export {}="{}"'.format(k,v)) 

parseYAML(config.interSessYAML)
parseYAML(config.projectYAML)


