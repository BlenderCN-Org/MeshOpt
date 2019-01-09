import h5py
import meshio
import matplotlib.pyplot as plt
import json
import t2listing
from mulgrids import *
from mpl_toolkits.mplot3d import Axes3D
import numpy
from collections import OrderedDict
from itertools import permutations


def plot_Tcell_vs_time(x_pos, y_pos, z_pos):

    target = None
    i = 0
    while target is None:
        if abs(points[i][0]-x_pos)<10000. and abs(points[i][1]-y_pos)<10000. and abs(points[i][2]-z_pos)<10000.:
            target = i
        i += 1
    i = 0
    cell_index = None
    while cell_index is None:
        if target in elements[i]:
            cell_index = i
        i += 1

    T_cell = [T[time][cell_index] for time in xrange(len(times))]
    plt.figure()
    plt.plot(times, T_cell)
    plt.show()


def build_isoT(T_target, time):

    data_x = []
    data_y = []
    data_z = []
    test =0
    for index in xrange(len(elements)):
        if abs(T[time][index] - T_target) < 0.0001:
            test += 1
            current_points = [points[vert] for vert in elements[index]]
            data_x.append(sum(current_points[j][0] for j in xrange(4))/4)
            data_y.append(sum(current_points[j][1] for j in xrange(4))/4)
            data_z.append(sum(current_points[j][2] for j in xrange(4))/4)

    return [data_x,data_y,data_z]


def compute_surfaces():

    surfs = numpy.zeros(len(tets))

    for con in el_neighbours:
        if len(el_neighbours[con]) == 1:
            face = [vert for vert in tets[el_neighbours[con][0]] if vert != neighbours[con][0]]
            current_z = [list_points[face[0]][2],list_points[face[1]][2],list_points[face[2]][2]]
            if all(abs(z+4000)<0.0001 for z in current_z):
                centroid = [sum(list_points[face[i]][j] for i in xrange(3))/3 for j in xrange(3)]
                if numpy.linalg.norm([centroid[i]-origin[i] for i in xrange(2)]) <= upflow_radius:
                    face_area = compute_surface_triangle(face)
                    surfs[el_neighbours[con][0]] = face_area

    return surfs


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


def compute_surface_triangle(triangle):

    a = numpy.linalg.norm([list_points[triangle[0]][i] - list_points[triangle[1]][i] for i in xrange(3)])
    b = numpy.linalg.norm([list_points[triangle[0]][i] - list_points[triangle[2]][i] for i in xrange(3)])
    c = numpy.linalg.norm([list_points[triangle[1]][i] - list_points[triangle[2]][i] for i in xrange(3)])
    s = (a + b + c) / 2
    return (s*(s-a)*(s-b)*(s-c)) ** 0.5

PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Model1/Tet/cyl/'
geofile = 'Basic_Model_Geo.dat'
filename = 'otet_m1.h5'
output = filename[:-3]
meshname = 'cyl_m1_opti.msh'
H5 = True
LISTING = False

if H5:
    file = h5py.File(PATH + filename, 'r')
    mesh = meshio.read(PATH + meshname)
    cells = mesh.cells
    points = mesh.points
    if 'hexahedron' in mesh.cells:
        elements = mesh.cells['hexahedron']
        type = 'hex'
    elif 'tetra' in mesh.cells:
        elements = mesh.cells['tetra']
        neighbours, el_neighbours = connection_mapping(elements)
        type = 'tet'


    cell_fields = file['cell_fields']
    cell_indexes = file['cell_index']
    cell_int_index = file['cell_interior_index']
    fields = file['fields']
    bulk_times = file['time']
    times = []
    for t in bulk_times:
        times.append(t[0])

    print len(times)
    T = cell_fields['fluid_temperature']
    Sg = cell_fields['fluid_vapour_saturation']
    P = cell_fields['fluid_pressure']
    if type == 'tet':
        mesh.cells['tetra'] = elements
        print len(elements), len(T[1])
        mesh.cell_data['tetra'] = {}
        mesh.cell_data['tetra']['T_0'] = T[0]
        mesh.cell_data['tetra']['T_1'] = T[1]
        mesh.cell_data['tetra']['P_0'] = P[0]
        mesh.cell_data['tetra']['P_1'] = P[1]
        mesh.cell_data['tetra']['Sg_0'] = Sg[0]
        mesh.cell_data['tetra']['Sg_1'] = Sg[1]
    if type == 'hex':
        mesh.cell_data['hexahedron'] = {}
        mesh.cell_data['hexahedron']['T_0'] = T[0]
        mesh.cell_data['hexahedron']['T_1'] = T[1]
        mesh.cell_data['hexahedron']['P_0'] = P[0]
        mesh.cell_data['hexahedron']['P_1'] = P[1]
        mesh.cell_data['hexahedron']['Sg_0'] = Sg[0]
        mesh.cell_data['hexahedron']['Sg_1'] = Sg[1]
    ####

    meshio.write(PATH + output + '.vtk', mesh)

elif LISTING:
    # file = t2listing.t2listing(PATH+filename, skip_tables=['generation', 'connection'])
    # geo = mulgrid(PATH + geofile)
    # file.write_vtk(geo, PATH + 'model1_hex')
    # mesh = meshio.read(PATH + meshname)
    # cells = mesh.cells
    # points = mesh.points
    # elements = cells['hexahedron']

    times = file.time
    print times

# PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Model1/Tet/'
# geofile = 'Basic_Model_Geo.dat'
# filename = 'm1_tet.json'
# output = 'heatflows.vtk'
# meshname = 'model1.msh'
# H5 = True
# LISTING = False
#
# origin = [5000.,5000.]
# upflow_radius = 1200.
#
#
# if H5:
#     with open(PATH + filename) as f:
#         data = json.load(f)
#
#     sources = data['source']
#     mesh = meshio.read(PATH + meshname)
#     list_points = mesh.points
#     tets = mesh.cells['tetra']
#     neighbours, el_neighbours = connection_mapping(tets)
#     surfaces = compute_surfaces()
#     hf = numpy.zeros(len(mesh.cells['tetra']))
#     for source in sources:
#         hf[source['cell']] = source['rate']
#     mesh.cell_data = {'tetra':{'heatflow': hf, 'surface': surfaces}}
#     meshio.write(PATH + output, mesh)
