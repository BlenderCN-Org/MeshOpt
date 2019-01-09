import json



def update_sources():

    for source in data['source']:
        new_heat = new_heatflow*source['rate']/old_heatflow
        source['rate'] = new_heat
        if source['name'][-2] == ' ':
            source['name'] = source['name'][:-2] + '0' + source['name'][-1]


PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Model1/Waiwera/CRU/'
json_name = 'm1_hex.json'

# MPATH = '/home/lmar626/Documents/Meshes/Basicmodels/Model1/Waiwera/CRU/'
# mesh_name = 'gm1_hex.exo'

with open(PATH+json_name) as f:
    data = json.load(f)

data['mesh']['filename'] = '/home/lmar626/Documents/Meshes/Basicmodels/Model1/Waiwera/CRU/gm1_hex.exo'

old_heatflow = 1
new_heatflow = 5e6
total_surf = 4640000.0

update_sources()

with open(PATH+'new'+json_name, 'w') as output:
    json.dump(data, output, sort_keys=True, indent=4, separators=(', ', ': '))