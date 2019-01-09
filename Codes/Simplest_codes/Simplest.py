import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import meshio
from scipy.spatial import Delaunay
from collections import OrderedDict
from itertools import permutations
import json


def connection_mapping(tet_list):
    "Dict that make correspond every faces to the two missing vertices of the elements sharing it"

    dict_faces_opti = OrderedDict()
    dict_faces = OrderedDict()
    for i in range(len(tet_list)):
        permutations_list = permutations(tet_list[i], 3)
        faces_list = []
        for perm in permutations_list:
            face = frozenset(perm)
            if face not in faces_list:
                faces_list.append(face)
                if face in dict_faces_opti:
                    dict_faces_opti[face] += [vertice for vertice in tet_list[i] if vertice not in face]
                    dict_faces[face] += [i]
                else:
                    dict_faces_opti[face] = [vertice for vertice in tet_list[i] if vertice not in face]
                    dict_faces[face] = [i]

    return dict_faces_opti, dict_faces


def fill_boundaries():
    j_bound = []
    if type == 'hex':
        for index, el in enumerate(elements):
            test = False
            for p in el:
                if points[p][0] == 0.:
                    test = True
            if test:
                j_bound.append({'cells': [index],'normal': [-1.,0.,0.]})
    if type == 'tet':
        for con in el_neighbours:
            if len(el_neighbours[con]) == 1:
                face = [vert for vert in elements[el_neighbours[con][0]] if vert != neighbours[con][0]]
                current_x = [points[face[0]][0], points[face[1]][0], points[face[2]][0]]
                if all(abs(x) < 0.001 for x in current_x):
                    j_bound.append({'cells': [el_neighbours[con][0]], 'normal': [-1., 0., 0.]})
    return j_bound


def create_sources():
    gen = []
    if type == 'hex':
        for index, el in enumerate(elements):
            test = False
            for p in el:
                if points[p][0] == 1.:
                    test = True
            if test:
                cid = str(index)
                if len(cid) == 1:
                    cur_name = '   p' + cid
                elif len(cid) == 2:
                    cur_name = '  p' + cid
                else:
                    cur_name = ' p' + cid
                gen.append({'cell': index, 'name': cur_name, 'rate': q/16})
    if type == 'tet':
        for con in el_neighbours:
            if len(el_neighbours[con]) == 1:
                face = [vert for vert in elements[el_neighbours[con][0]] if vert != neighbours[con][0]]
                current_surf = compute_surface_triangle(face)
                current_x = [points[face[0]][0], points[face[1]][0], points[face[2]][0]]
                if all(abs(x-1) < 0.001 for x in current_x):
                    cid = str(el_neighbours[con][0])
                    if len(cid) == 1:
                        cur_name = '   p' + cid
                    elif len(cid) == 2:
                        cur_name = '  p' + cid
                    else:
                        cur_name = ' p' + cid
                    gen.append({'cell': el_neighbours[con][0], 'name': cur_name, 'rate': current_surf*q})
    return gen


def compute_surface_triangle(triangle):

    a = np.linalg.norm([points[triangle[0]][i] - points[triangle[1]][i] for i in xrange(3)])
    b = np.linalg.norm([points[triangle[0]][i] - points[triangle[2]][i] for i in xrange(3)])
    c = np.linalg.norm([points[triangle[1]][i] - points[triangle[2]][i] for i in xrange(3)])
    s = (a + b + c) / 2
    return (s*(s-a)*(s-b)*(s-c)) ** 0.5


PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Simplest/bigger/'
model = 'tet_opti.msh'
out = model[:-4]
title = 'bigger_mesh_opti'

P = 5e5
T = 20.0
q = -0.5e-5
K = 1e-15
poro = 0.1


mesh = meshio.read(PATH+model)
points = mesh.points
if 'hexahedron' in mesh.cells:
    elements = mesh.cells['hexahedron']
    type = 'hex'
elif 'tetra' in mesh.cells:
    elements = mesh.cells['tetra']
    neighbours, el_neighbours = connection_mapping(elements)
    type = 'tet'


data = {}
data['title'] = title
data['boundaries'] = [{'region': 1, 'primary': [P], 'faces': fill_boundaries()}]
data['eos'] = {'name': 'w', 'temperature': T}
data['gravity'] = 0.
data['initial'] = {'primary': [P]}
data['thermodynamics'] = 'ifc67'
data['mesh'] = {'filename': PATH + model, 'permeability_angle': 0.0}
data['output'] = {'filename': PATH + out + '.h5', 'final': True, 'frequency': 0}
dom = {'cells': [ID for ID in xrange(len(elements))], 'name': '  dom', 'permeability': [K,K,K], 'porosity': poro}
data['rock'] = {'types': [dom]}
data['source'] = create_sources()
data['time'] = {"start": 0, "step": {"adapt": {"amplification": 2, "maximum": 8,  "method": "iteration", "minimum": 5, "on": True,"reduction": 0.2},"maximum": {"number": 500,"size": None},"method": "beuler","size": 100000,"solver": {"linear": {"preconditioner": {"type": "bjacobi"},"type": "bcgs"},"nonlinear": {"maximum": {"iterations": 8},"minimum": {"iterations": 1},"tolerance": {"function": {"absolute": None,"relative": None}}}}},"stop": 1e+15}



with open(PATH + out + '.json', 'w') as output:
    json.dump(data, output, sort_keys=True, indent=4, separators=(', ', ': '))

# PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Simplest/'
# model = 'tet'
# exts = ['.msh', '.vtk']
# base = [0.,0.5,1.]

# mesh = meshio.read(PATH+model+'.exo')
# for ext in exts:
#     meshio.write(PATH+model+ext,mesh)


#Hex Cube mesh
# indexes = np.linspace(0,26,27)
# X,Y,Z = np.meshgrid(base,base,base, indexing='ij')
# indexes = indexes.reshape(X.shape)
# XX, YY, ZZ = X.flatten(),Y.flatten(),Z.flatten()
# points = np.array([[XX[p],YY[p],ZZ[p]] for p in xrange(27)])
# hex = []
# for i in xrange(2):
#     for j in xrange(2):
#         for k in xrange(2):
#             hex.append([int(indexes[i,j,k]),int(indexes[i+1,j,k]),int(indexes[i+1,j+1,k]),int(indexes[i,j+1,k]),int(indexes[i,j,k+1]),int(indexes[i+1,j,k+1]),int(indexes[i+1,j+1,k+1]),int(indexes[i,j+1,k+1])])
# hex = np.array(hex)
#
# for ext in exts:
#     meshio.write_points_cells(PATH+model+ext, points, {'hexahedron': hex})

#Tetnotopti
# X,Y,Z = np.meshgrid(base,base,base, indexing='ij')
# XX, YY, ZZ = X.flatten(),Y.flatten(),Z.flatten()
# points = [[XX[p],YY[p],ZZ[p]] for p in xrange(27)]
#
# for i in base:
#     rc = [0.1*np.random.randint(1,5) for it in xrange(8)]
#     print rc
#     points += [[i,rc[0], rc[1]],[i,0.5+rc[2],rc[3]],[i,0.5+rc[4],0.5+rc[5]],[i,rc[6],0.5+rc[7]]]
#     rc = [0.1*np.random.randint(1,5) for it in xrange(8)]
#     print rc
#     points += [[rc[0],i,rc[1]],[0.5+rc[2],i,rc[3]],[0.5+rc[4],i,0.5+rc[5]],[rc[6],i,0.5+rc[7]]]
#     rc = [0.1 * np.random.randint(1, 5) for it in xrange(8)]
#     print rc
#     points += [[rc[0], rc[1], i],[0.5 + rc[2], rc[3], i], [0.5 + rc[4], 0.5 + rc[5], i], [rc[6], 0.5 + rc[7], i]]
#
# rc = [0.1 * np.random.randint(1, 5) for it in xrange(12)]
# print rc
# points += [[rc[0],rc[1],rc[2]],[0.5+rc[3],rc[4],rc[5]],[0.5+rc[6],0.5+rc[7],rc[8]],[rc[9],0.5+rc[10],rc[11]]]
# rc = [0.1 * np.random.randint(1, 5) for it in xrange(12)]
# print rc
# points += [[rc[0],rc[1],0.5+rc[2]],[0.5+rc[3],rc[4],0.5+rc[5]],[0.5+rc[6],0.5+rc[7],0.5+rc[8]],[rc[9],0.5+rc[10],0.5+rc[11]]]
# points = np.array(points)
# tets = Delaunay(points).simplices
# print max(points[:,0]), max(points[:,1]), max(points[:,2])
# for ext in exts:
#     meshio.write_points_cells(PATH+model+ext, points, {'tetra': tets})
