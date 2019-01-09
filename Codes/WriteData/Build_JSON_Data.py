import json
import meshio
import numpy
from collections import OrderedDict
from itertools import permutations


def get_cross_product(vector1, vector2):
    "Returns the cross product of two vector to access the normal vector"

    n_vect = [vector1[1] * vector2[2] - vector1[2] * vector2[1],
              -1 * (vector1[0] * vector2[2] - vector1[2] * vector2[0]),
              vector1[0] * vector2[1] - vector1[1] * vector2[0]]
    return n_vect


def get_plane_equation(point_list):
    "Returns the coefs (a,b,c,d) of the equation of a given plane"

    vect1 = [point_list[1][i] - point_list[0][i] for i in xrange(3)]
    vect2 = [point_list[1][i] - point_list[2][i] for i in xrange(3)]
    norm1 = numpy.linalg.norm(vect1)
    norm2 = numpy.linalg.norm(vect2)
    vect1 = [coord/norm1 for coord in vect1]
    vect2 = [coord/norm2 for coord in vect2]
    cross_product = get_cross_product(vect1, vect2)
    coefs = [cross_product[0], cross_product[1], cross_product[2], - (
            cross_product[0] * point_list[1][0] + cross_product[1] * point_list[1][1] + cross_product[2] *
            point_list[1][2])]
    # print coefs, sympy.Plane(sympy.Point3D(point_list[1]),sympy.Point3D(point_list[0]),sympy.Point3D(point_list[2])).equation()
    return coefs


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
    j_bounds = []
    if type == 'tet':
        for con in el_neighbours:
            if len(el_neighbours[con]) == 1:
                face = [vert for vert in elements[el_neighbours[con][0]] if vert != neighbours[con][0]]
                current_z = [list_points[face[0]][2], list_points[face[1]][2], list_points[face[2]][2]]
                if all(abs(z) < 0.001 for z in current_z):
                    j_bounds.append({'cells':[el_neighbours[con][0]],'normal': [0.,0.,1.]})
    else:
        for index, el in enumerate(elements):
            test = False
            for p in el:
                if list_points[p][2] == 0.:
                    test = True
            if test:
                j_bounds.append({'cells': [index],'normal': [1.,0.,0.]})
    return j_bounds


def fill_rocktypes():
    # for i, dat in enumerate(cell_data['tetra']['gmsh:physical']):
    #     if dat == 3:
    #         res['cells'].append(i)
    #     elif dat == 2:
    #         cap['cells'].append(i)
    #     else:
    #         surf['cells'].append(i)
    if type == 'tet':
        for i, el in enumerate(elements):
            centroid = [sum([list_points[v][vv] for v in el])/4 for vv in xrange(3)]
            found = False
            if faults:
                #m3
                # d1 = eq1[0] * centroid[0] + eq1[1] * centroid[1] + eq1[2] * centroid[2] + eq1[3]
                # d2 = eq2[0] * centroid[0] + eq2[1] * centroid[1] + eq2[2] * centroid[2] + eq2[3]
                # if d1 <= 0 and d2 >= 0:
                #     fault['cells'].append(i)
                #     found = True
                #m2
                if 4900. <= centroid[0] <= 5100. or 4900. <= centroid[1] <= 5100.:
                    fault['cells'].append(i)
                    found = True
            if centroid[2] < -1800. and not found:
                res['cells'].append(i)
            elif -1800. <= centroid[2] < -600. and not found:
                cap['cells'].append(i)
            elif -600. <= centroid[2] and not found:
                surf['cells'].append(i)

    else:
        for i, el in enumerate(elements):
            dat = numpy.array([list_points[ind] for ind in el])
            centroid = [(max(dat[:,0])+min(dat[:,0]))/2, (max(dat[:,1])+min(dat[:,1]))/2, (max(dat[:,2])+min(dat[:,2]))/2]
            found = False
            if faults:
                # m3
                # d1 = eq1[0] * centroid[0] + eq1[1] * centroid[1] + eq1[2] * centroid[2] + eq1[3]
                # d2 = eq2[0] * centroid[0] + eq2[1] * centroid[1] + eq2[2] * centroid[2] + eq2[3]
                # if d1 <= 0 and d2 >= 0:
                #     fault['cells'].append(i)
                #     found = True
                # m2
                if 4900. <= centroid[0] <= 5100. or 4900. <= centroid[1] <= 5100.:
                    fault['cells'].append(i)
                    found = True
            if centroid[2] < -1800. and not found:
                res['cells'].append(i)
            elif -1800. <= centroid[2] < -600.and not found:
                cap['cells'].append(i)
            elif -600. <= centroid[2] and not found:
                surf['cells'].append(i)
    if faults:
        return [res,cap, surf, fault]
    else:
        return [res, cap, surf]


def create_sources():

    gen = []
    faces = []
    cell_ID = []
    if type == 'tet':
        for con in el_neighbours:
            if len(el_neighbours[con]) == 1:
                face = [vert for vert in elements[el_neighbours[con][0]] if vert != neighbours[con][0]]
                current_z = [list_points[face[0]][2],list_points[face[1]][2],list_points[face[2]][2]]
                if all(abs(z+4000)<10. for z in current_z):
                    centroid = [sum(list_points[face[i]][j] for i in xrange(3))/3 for j in xrange(3)]
                    if numpy.linalg.norm([centroid[i]-origin[i] for i in xrange(2)]) <= upflow_radius:
                        faces.append(face)
                        cell_ID.append(el_neighbours[con][0])
        total_area = sum([compute_surface_triangle(face) for face in faces])
        for id,face in enumerate(faces):
            cur_name = str(len(gen))
            if len(cur_name) == 1:
                cur_name = 's  0' + cur_name
            elif len(cur_name) == 2:
                cur_name = 's  ' + cur_name
            elif len(cur_name) == 3:
                cur_name = 's ' + cur_name
            else:
                cur_name = 's' + cur_name
            face_area = compute_surface_triangle(face)
            heat_rate = total_heat * face_area / total_area
            gen.append({'cell': cell_ID[id],'name': cur_name, 'rate': heat_rate, 'component': 2})
    else:
        for i, el in enumerate(elements):
            test = False
            dat = numpy.array([list_points[ind] for ind in el])
            for p in dat:
                if abs(p[2]+4000.) < 10.:
                    test = True
            if test:
                centroid = [(max(dat[:,0])+min(dat[:,0]))/2, (max(dat[:,1])+min(dat[:,1]))/2]
                if numpy.linalg.norm([centroid[cc] - origin[cc] for cc in xrange(2)]) <= upflow_radius:
                    faces.append((max(dat[:,0])-min(dat[:,0]))*(max(dat[:,1])-min(dat[:,1])))
                    cell_ID.append(i)
        total_area = sum(faces)
        for id, face in enumerate(faces):
            cur_name = str(len(gen))
            if len(cur_name) == 1:
                cur_name = 's  0' + cur_name
            elif len(cur_name) == 2:
                cur_name = 's  ' + cur_name
            elif len(cur_name) == 3:
                cur_name = 's ' + cur_name
            else:
                cur_name = 's' + cur_name
            heat_rate = total_heat * face / total_area
            gen.append({'cell': cell_ID[id], 'name': cur_name, 'rate': heat_rate, 'component': 2})
    return gen


def compute_surface_triangle(triangle):

    a = numpy.linalg.norm([list_points[triangle[0]][i] - list_points[triangle[1]][i] for i in xrange(3)])
    b = numpy.linalg.norm([list_points[triangle[0]][i] - list_points[triangle[2]][i] for i in xrange(3)])
    c = numpy.linalg.norm([list_points[triangle[1]][i] - list_points[triangle[2]][i] for i in xrange(3)])
    s = (a + b + c) / 2
    return (s*(s-a)*(s-b)*(s-c)) ** 0.5


PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Model2/Tet/cyl/'
type = 'tet'
faults = True
# f1 = [[2000.,0.,0.],[6000.,0.,-4000.],[6000.,10000.,-4000.]]
# f2 = [[2400.,0.,0.],[6400.,0.,-4000.],[6400.,10000.,-4000.]]
# eq1 = get_plane_equation(f1)
# eq2 = get_plane_equation(f2)


filename = 'cyl_m2.msh'
output = 'tet_m2'
title = 'tet_m2_naturalstate'
T_res = 20.
P_res = 1.013e5
upflow_radius = 1200.
origin = [5000.,5000.]
total_heat = 5e6


res = {'cells': [], 'name': '  res', 'permeability': [1e-15, 1e-15, 1e-15], 'porosity': 0.2,"specific_heat": 900.0,"density": 2600.0,"wet_conductivity": 1.5,"dry_conductivity": 1.5}
cap = {'cells': [], 'name': '  cap', 'permeability': [0.1e-15, 0.1e-15, 0.1e-15], 'porosity': 0.1,"specific_heat": 900.0,"density": 2600.0,"wet_conductivity": 1.5,"dry_conductivity": 1.5}
# cap = {'cells': [], 'name': '  cap', 'permeability': [1e-15, 1e-15, 1e-15], 'porosity': 0.2,"specific_heat": 900.0,"density": 2600.0,"wet_conductivity": 1.5,"dry_conductivity": 1.5}
surf = {'cells': [], 'name': ' surf', 'permeability': [1e-15, 1e-15, 1e-15], 'porosity': 0.2,"specific_heat": 900.0,"density": 2600.0,"wet_conductivity": 1.5,"dry_conductivity": 1.5}
if faults:
    fault = {'cells':[], 'name': 'fault','permeability': [1e-14, 1e-14, 1e-14], 'porosity': 0.1,"specific_heat": 900.0,"density": 2600.0,"wet_conductivity": 1.5,"dry_conductivity": 1.5}


mesh = meshio.read(PATH + filename)
cells = mesh.cells
list_points = mesh.points
cell_data = mesh.cell_data
# print cell_data['tetra']['gmsh:physical'][2]
if type == 'tet':
    elements = cells['tetra']
    neighbours, el_neighbours = connection_mapping(elements)
else:
    elements = cells['hexahedron']
# list_bounds = {'plane': [[0., 0., 1., 0.]]}
data = {}


data['boundaries'] = [{'region': 1, 'primary':[P_res,T_res], 'faces': fill_boundaries()}]
data['eos'] = {'name': 'we'}
data['gravity'] = 9.81
data['initial'] = {'primary':[P_res, T_res], 'region': 1}
data['thermodynamics'] = 'ifc67'
# data['initial'] = {'filename': PATH + 'model1.h5', 'index': 7}
data['mesh'] = {'filename': PATH + filename, 'permeability_angle': 0.0}
data['output'] = {'filename': PATH + output + '.h5', 'final': True, 'frequency': 0}
data['rock'] = {'capillary_pressure': {'pressure': 0.0, 'saturation_limits': [0., 1.], 'type': 'linear'}, 'types': fill_rocktypes()}
data['source'] = create_sources()
data['time'] = {"start": 0, "step": {"adapt": {"amplification": 2, "maximum": 8,  "method": "iteration", "minimum": 5, "on": True,"reduction": 0.2},"maximum": {"number": 500,"size": None},"method": "beuler","size": 100000,"solver": {"linear": {"preconditioner": {"type": "bjacobi"},"type": "bcgs"},"nonlinear": {"maximum": {"iterations": 8},"minimum": {"iterations": 1},"tolerance": {"function": {"absolute": None,"relative": None}}}}},"stop": 1e+15}
data['title'] = title

with open(PATH + output + '.json', 'w') as output:
    json.dump(data, output, sort_keys=True, indent=4, separators=(', ', ': '))
