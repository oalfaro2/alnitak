import json
import os
import sys
#Writes given key to config.json to avoid having to nano the file
#For use in command line only
#python /path/alnitak/api_key.py insertkeyhere
key = sys.argv[1]

cwd, path = os.path.split(sys.argv[0])
filename = os.path.join(cwd, 'config.json')
with open(filename, 'r') as f:
    config = json.load(f)
    config['API_Key'] = key

os.remove(filename)
with open(filename, 'w') as f:
    json.dump(config, f, indent=4)
