import meshio
import Optimization_w_CS as ocs
import numpy as np

PATH = '/home/mltn/Documents/Meshes/Basicmodels/Simplest/bigger/'
# bound = 'bdata_m1.txt'
names = ['tet.vtk', 'otet_m1.vtk']

mesh = meshio.read(PATH+names[0])
points = mesh.points
# print len(points)
el = mesh.cells['tetra']
# print len(el)
vn, eln = ocs.connection_mapping(el)
vars = np.arange(0,len(points)+1,1)
cons = ocs.Con(vn, eln, el, vars)
indexes = cons.get_indexes()
node = ocs.build_node_objects(points, vars, {})
devs = ocs.quality_control(node,cons)
# print len(devs), indexes
max_dev = np.zeros(len(el))
min_dev = np.ones(len(el))
avg_dev_num = np.zeros(len(el))
avg_dev_count = np.zeros(len(el))
i_ind = 0
for i_dev, dev in enumerate(devs):
    while len(indexes[i_ind]) < 2:
        i_ind += 1
    for index in indexes[i_ind]:
        # print index
        if max_dev[index] < dev:
            max_dev[index] = dev
        if min_dev[index] > dev:
            min_dev[index] = dev
        avg_dev_num[index] += dev
        avg_dev_count[index] += 1
    i_ind += 1

zeros = 0
for count in avg_dev_count:
    if count == 0:
        zeros += 1
print zeros
avg_dev = np.zeros(len(el))
for ind, num in enumerate(avg_dev_num):
    avg_dev[ind] = num/avg_dev_count[ind]

mesh.cell_data['tetra']['max_dev'] = max_dev
mesh.cell_data['tetra']['min_dev'] = min_dev
mesh.cell_data['tetra']['avg_dev'] = avg_dev
meshio.write(PATH+'test.vtk', mesh)