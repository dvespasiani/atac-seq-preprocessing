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
