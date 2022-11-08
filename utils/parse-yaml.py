import argparse
import os
import yaml as ym

def config_arguments():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--project_dir', type=str, default=os.getcwd(), help='Root of your project directory')
    parser.add_argument('-y', '--yaml_file', type=str, default=None, help='the yaml file containing the directives you want to parse')
 
    args = parser.parse_args()  
    return(args)

def parseYAML(y):
    with open(y, 'r') as file:
        out = ym.safe_load(file)
    return(out)

config = config_arguments()
yaml_file = os.path.join(config.project_dir + config.yaml_file)

parsed_yaml = parseYAML(yaml_file)

if "interactive_sess" in parsed_yaml:
    interactive_sess = parsed_yaml['interactive_sess']
    interactive_sess['modules'] = config.project_dir + '/' + interactive_sess['modules']
    for k, v in interactive_sess.items():    
        print('export {}="{}"'.format(k,v)) 


# def parseYAML(y):
#     with open(y, 'r') as file:
#         yamlFile = yaml.safe_load(file)  
    # for k, v in yamlFile.items():    
        # print('export {}="{}"'.format(k,v)) 

# parseYAML(config.interSessYAML)
# parseYAML(config.projectYAML)




# parser.add_argument('-b', '--baseDir', type=str, default=os.getcwd(), help='directory where this script is stored, I recommend it to be your projects root dir')
# parser.add_argument('-p', '--projectName', type=str, default=None, help='Name of your project directory, this is relative to your base working dir')
# parser.add_argument('-m', '--modules', type=str, default='config/modules.txt', help='Path to your modules.txt file relative to the project root dir')

# args = parser.parse_args()  

# project_dir = args.baseDir + '/' + args.projectName
# modules = project_dir + args.modules

# project_config_yaml = os.path.join(project_dir + "config/project_config.yaml")



# homeDir = os.path.expanduser("/stornext/General/data/user_managed/grpu_jchoi_0/projects/davide/") # your home directory

# def config_arguments():
    
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-bd', '--baseDir', type=str, default=homeDir, help='path to base working directory, default is the home directory')
#     parser.add_argument('-pd', '--projectDir', type=str, default='atac-pipeline/', help='path to project dir (relative to base dir)')
    
#     args = parser.parse_args()  
#     args.interSessYAML = os.path.join(args.baseDir,args.projectDir, "config/cluster_config.yaml")
#     args.projectYAML = os.path.join(args.baseDir,args.projectDir, "config/project_config.yaml") 
#     return args


# config = config_arguments()

# def parseYAML(y):
#     with open(y, 'r') as file:
#         yamlFile = yaml.safe_load(file)  
#     for k, v in yamlFile.items():    
#         print('export {}="{}"'.format(k,v)) 

# parseYAML(config.interSessYAML)
# parseYAML(config.projectYAML)






