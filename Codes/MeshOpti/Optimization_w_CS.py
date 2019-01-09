import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import meshio
from sys import stdout
from itertools import permutations
import numpy
from scipy.optimize import minimize
import scipy.interpolate
import scipy.signal
from datetime import datetime
from collections import OrderedDict
#from pathos.multiprocessing import ProcessingPool
from multiprocessing import Process #, Queue
from fmq import Queue
from math import cos, sin, acos
from time import clock
from itertools import combinations
from sklearn.cluster import KMeans
from scipy.spatial import distance_matrix
# import sympy

# Input variables

# Inputs
toggle_complex = False
if toggle_complex:
    TOPO_PATH = '/home/lmar626/Documents/Meshes/TestProblems/Scaling/'
    complex_bounds = [TOPO_PATH + 'tdata.txt',
                      [[0., 0., 0.], [10000., 0., 0.], [10000., 10000., 0.], [0., 10000., 0.]]]
    # complex_bounds = [TOPO_PATH + 'TVZ_tdata.txt',[[0., 0., 0.], [160965., 0., 0.], [160965., 86605., 0.], [0., 86605., 0.]]]
else:
    complex_bounds = None

internal_bound = True
BOUND_DAT_PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Model3/Tet/'
plane_bounds = 'bdata_m3.txt'
MESH_PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Model3/Tet/cyl/'
filenames = ['cyl_m3.msh'] #
# filename = '1733.msh'
# Maximum number of iteration during optimization
multiprocess = True
nb_process = 3
# chunk_size = 10000
iter = 1000
epsilon = 10.
#results
SAVE_PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Model3/Tet/cyl/'
fname = 'Opti_m3_result_'
ext = '.txt'


# Geometry


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


def get_vect_direction(point_list):

    # tasks = [[0,1],[1,2],[2,3],[3,0]]
    # vectors = []
    # for task in tasks:
    #     origin = point_list[task[0]]
    #     bulk_vector = [point_list[task[1]][i] - origin[i] for i in xrange(3)]
    #     for ind,coef in enumerate(bulk_vector):
    #         if coef > 0:
    #             bulk_vector[ind] = 1.0
    #         elif coef < 0:
    #             bulk_vector[ind] = -1.0
    #     vectors.append([origin,bulk_vector])
    u = [point_list[1][uu] - point_list[0][uu] for uu in xrange(3)]
    v = [point_list[1][vv] - point_list[2][vv] for vv in xrange(3)]
    norm_u = numpy.linalg.norm(u)
    norm_v = numpy.linalg.norm(v)
    u = [cu/norm_u for cu in u]
    v = [cv/norm_v for cv in v]
    u1 = numpy.dot(u,point_list[1])
    u2 = numpy.dot(u,point_list[0])
    v1 = numpy.dot(v,point_list[1])
    v2 = numpy.dot(v,point_list[2])
    if u1 < u2:
        u_check = [u, u1, u2]
    else:
        u_check = [u, u2, u1]
    if v1 < v2:
        v_check = [v, v1, v2]
    else:
        v_check = [v, v2, v1]

    return [u_check, v_check]


# Class Con

class Con:

    def __init__(self, dict_con, el_con, el_list,vars_list):
        self.faces_ = []
        self.vertices_ = []
        self.indexes_ = []
        self.build_cons(dict_con, el_con, el_list, vars_list)

    def build_cons(self, dict_con, el_con, el_list, vars_list):
        vars_list = set(vars_list)
        size = len(dict_con)
        counter = 0
        self.faces_ = [[]] * size
        self.vertices_ = [[]] * size
        self.indexes_ = [[]] * size
        for con in dict_con:
            vertices = set(dict_con[con])
            if not vars_list.isdisjoint(vertices) or not vars_list.isdisjoint(con):
                self.vertices_[counter] = dict_con[con]
                self.indexes_[counter] = el_con[con]
                face = [vert for vert in el_list[el_con[con][0]] if vert != dict_con[con][0]]
                face.sort()
                self.faces_[counter] = face
                counter += 1

    def get_faces(self): return self.faces_

    def get_vertices(self): return self.vertices_

    def get_indexes(self): return self.indexes_


# Class Nodes


class Node:

    def __init__(self, point_index, point_coord, node_bound_list, node_bound_types, vars_list, dict_bounds):
        self.global_index_ = point_index
        self.x_ = point_coord[0]
        self.y_ = point_coord[1]
        self.z_ = point_coord[2]
        self.alloc_vars_ = []
        self.neighbours_ = set([])
        self.bound_ = node_bound_list
        self.bound_types_ = node_bound_types

        if point_index in vars_list:
            if 0 < len(self.bound_) < 3:
                if 'complex' in self.bound_types_:
                    self.determine_complex_freedom(dict_bounds)
                else:
                    self.determine_freedom(dict_bounds)
            elif len(self.bound_) == 0:
                self.deg_freedom_ = 3
                self.freedom_ = [1, 1, 1]
            else:
                self.deg_freedom_ = 0
                self.freedom_ = None
        else:
            self.deg_freedom_ = 0
            self.freedom_ = None

        if self.freedom_ is not None:
            if self.deg_freedom_ == 3:
                self.coefs_ = None
            else:
                if 'complex' in self.bound_types_ and self.deg_freedom_ == 2:
                    self.coefs_ = None
                elif 'complex' not in self.bound_types_:
                    self.determine_coefs(dict_bounds)

    def determine_complex_freedom(self,dict_bounds):
        #TODO change alg so that it works with self.bound_ only holding the index of the corresponding bound (Normally done)
        self.freedom_ = [1,1,1]
        nb_bound = len(self.bound_)
        if nb_bound == 1:
            self.freedom_[dict_bounds[self.bound_types_[0]][self.bound_[0]][2]] = 0
            self.deg_freedom_ = 2
        elif nb_bound == 2:
            if 'plane' in self.bound_types_:
                loc_axis = [0, 1, 2]
                i = 0
                while dict_bounds[self.bound_types_[0]][self.bound_[0]][0][i] == 0:
                    i += 1

                loc_axis = loc_axis[0:i] + loc_axis[i + 1:3]
                coefs = [dict_bounds[self.bound_types_[0]][self.bound_[0]][0][j] * (-dict_bounds[self.bound_types_[1]][self.bound_[1]][1][i] / dict_bounds[self.bound_types_[0]][self.bound_[0]][0][i]) - dict_bounds[self.bound_types_[1]][self.bound_[1]][1][j] for j
                         in range(4) if j != i]
                j = 0
                while coefs[j] == 0:
                    j += 1
                self.freedom_[i] = 0
                self.freedom_[loc_axis[j]] = 0
                self.deg_freedom_ = 1
            else:
                loc_axis = [0, 1, 2]
                i = 0
                while dict_bounds[self.bound_types_[0]][self.bound_[0]][1][i] == 0:
                    i += 1
                loc_axis = loc_axis[0:i] + loc_axis[i + 1:3]
                coefs = [dict_bounds[self.bound_types_[0]][self.bound_[0]][1][j] * (-dict_bounds[self.bound_types_[1]][self.bound_[1]][1][i] / dict_bounds[self.bound_types_[0]][self.bound_[0]][1][i]) - dict_bounds[self.bound_types_[1]][self.bound_[1]][1][j] for j
                         in range(4) if j != i]
                j = 0
                while coefs[j] == 0:
                    j += 1
                self.freedom_[i] = 0
                self.freedom_[loc_axis[j]] = 0
                self.deg_freedom_ = 1

    def determine_freedom(self,dict_bounds):
        # TODO change alg so that it works with self.bound_ only holding the index of the corresponding bound (Normally done)
        self.freedom_ = [1, 1, 1]
        nb_bound = len(self.bound_)
        if nb_bound == 1:
            i = 0
            while dict_bounds[self.bound_types_[0]][self.bound_[0]][0][i] == 0:
                i += 1
            self.freedom_[i] = 0
            self.deg_freedom_ = 2
        elif nb_bound == 2:
            loc_axis = [0, 1, 2]
            i = 0
            while dict_bounds[self.bound_types_[0]][self.bound_[0]][0][i] == 0:
                i += 1
            loc_axis = loc_axis[0:i] + loc_axis[i + 1:3]
            coefs = [dict_bounds[self.bound_types_[0]][self.bound_[0]][0][j] * (-dict_bounds[self.bound_types_[1]][self.bound_[1]][0][i] / dict_bounds[self.bound_types_[0]][self.bound_[0]][0][i]) - dict_bounds[self.bound_types_[1]][self.bound_[1]][0][j] for j
                     in range(4) if j != i]
            j = 0
            while coefs[j] == 0:
                j += 1
            # print self.bound_, loc_axis, j, coefs, i
            self.freedom_[i] = 0
            self.freedom_[loc_axis[j]] = 0
            self.deg_freedom_ = 1

    def determine_coefs(self,dict_bounds):
        # TODO change alg so that it works with self.bound_ only holding the index of the corresponding bound (Normally done)
        index_1 = 0
        index_2 = 0
        for i in range(3):
            if self.freedom_[i] == 1:
                index_1 = i
            else:
                index_2 = i

        if self.deg_freedom_ == 1:
            k = 0
            while dict_bounds[self.bound_types_[0]][self.bound_[0]][0][k] == 0:
                k += 1
            m = 3
            for j in range(3):
                if j != index_1 and j != k:
                    m = j
            l_coefs2 = [
                (dict_bounds[self.bound_types_[1]][self.bound_[1]][0][k] * dict_bounds[self.bound_types_[0]][self.bound_[0]][0][rang] / dict_bounds[self.bound_types_[0]][self.bound_[0]][0][k] - dict_bounds[self.bound_types_[1]][self.bound_[1]][0][rang]) / (
                        -dict_bounds[self.bound_types_[1]][self.bound_[1]][0][k] * dict_bounds[self.bound_types_[0]][self.bound_[0]][0][m] / dict_bounds[self.bound_types_[0]][self.bound_[0]][0][k] + dict_bounds[self.bound_types_[1]][self.bound_[1]][0][m]) for rang in
                [index_1, 3]]
            l_coefs1 = [
                -dict_bounds[self.bound_types_[0]][self.bound_[0]][0][m] * l_coefs2[0] / dict_bounds[self.bound_types_[0]][self.bound_[0]][0][k] - dict_bounds[self.bound_types_[0]][self.bound_[0]][0][index_1] / dict_bounds[self.bound_types_[0]][self.bound_[0]][0][k],
                -dict_bounds[self.bound_types_[0]][self.bound_[0]][0][m] * l_coefs2[1] / dict_bounds[self.bound_types_[0]][self.bound_[0]][0][k] - dict_bounds[self.bound_types_[0]][self.bound_[0]][0][3] / dict_bounds[self.bound_types_[0]][self.bound_[0]][0][k]]
            if k < m:
                self.coefs_ = [l_coefs1, l_coefs2]
            else:
                self.coefs_ = [l_coefs2, l_coefs1]

        else:
            loc_vars = [0, 0]
            loc_index = 0
            for j in range(3):
                if j != index_2:
                    loc_vars[loc_index] = j
                    loc_index += 1
            loc_vars.sort()
            self.coefs_ = [-dict_bounds[self.bound_types_[0]][self.bound_[0]][0][loc_vars[0]] / dict_bounds[self.bound_types_[0]][self.bound_[0]][0][index_2],
                           -dict_bounds[self.bound_types_[0]][self.bound_[0]][0][loc_vars[1]] / dict_bounds[self.bound_types_[0]][self.bound_[0]][0][index_2],
                           -dict_bounds[self.bound_types_[0]][self.bound_[0]][0][3] / dict_bounds[self.bound_types_[0]][self.bound_[0]][0][index_2]]

    def update_coefs(self, new_coefs):
        self.coefs_ = new_coefs

    def allocation_var_indexes(self, current_lenght):

        if self.deg_freedom_ == 3:
            self.alloc_vars_ = [current_lenght, current_lenght + 1, current_lenght + 2]
        elif self.deg_freedom_ == 2:
            self.alloc_vars_ = [current_lenght, current_lenght + 1]
        else:
            self.alloc_vars_ = [current_lenght]

    def update_coords(self, current_vars, dict_bounds):
        new_coords = [0, 0, 0]
        if self.deg_freedom_ == 3:
            new_coords[0] = current_vars[self.alloc_vars_[0]]
            new_coords[1] = current_vars[self.alloc_vars_[1]]
            new_coords[2] = current_vars[self.alloc_vars_[2]]
        elif self.deg_freedom_ == 2:
            loc_index = 0
            if 'complex' in self.bound_types_:
                for i in xrange(3):
                    if self.freedom_[i] == 1:
                        new_coords[i] = current_vars[self.alloc_vars_[loc_index]]
                        loc_index += 1
                    else:
                        new_coords[i] = scipy.interpolate.bisplev(current_vars[self.alloc_vars_[0]], current_vars[self.alloc_vars_[1]], dict_bounds[self.bound_types_[0]][self.bound_[0]][0])
            else:
                for i in range(3):
                    if self.freedom_[i] == 1:
                        new_coords[i] = current_vars[self.alloc_vars_[loc_index]]
                        loc_index += 1
                    else:
                        new_coords[i] = self.coefs_[0] * current_vars[self.alloc_vars_[0]] + self.coefs_[1] * current_vars[
                            self.alloc_vars_[1]] + self.coefs_[2]
        else:
            loc_index = 0
            if 'complex' in self.bound_types_:
                if 'plane' in self.bound_types_:
                    # TODO change the workflow so that the new point is properly projected onto the intersection curve (Normally done)
                    for i in xrange(3):
                        w_new = scipy.interpolate.splev(current_vars[self.alloc_vars_[0]], self.coefs_[2])
                        if self.freedom_[i] == 1:
                            new_coords[i] = current_vars[self.alloc_vars_[0]]
                        else:
                            new_coords[i] = self.coefs_[loc_index][0]*current_vars[self.alloc_vars_[0]] + self.coefs_[loc_index][1]*w_new + self.coefs_[loc_index][2]
                            loc_index += 1

                else:
                    #TODO change the workflow so that the new point is properly projected onto the intersection curve
                    for i in range(3):
                        if self.freedom_[i] == 1:
                            new_coords[i] = current_vars[self.alloc_vars_[0]]
                        else:
                            new_coords[i] = self.coefs_[loc_index][0] * current_vars[self.alloc_vars_[0]] + \
                                            self.coefs_[loc_index][1]
                            loc_index += 1
                    for ind in xrange(1,3,1):
                        loc_index = [i for i in xrange(3) if i != self.bound_[ind][2]]
                        new_coords[self.bound_[ind][2]] = scipy.interpolate.bisplev(new_coords[loc_index[0]],
                                                                              new_coords[loc_index[1]],
                                                                              self.bound_[ind][0])
            else:
                for i in range(3):
                    if self.freedom_[i] == 1:
                        new_coords[i] = current_vars[self.alloc_vars_[0]]
                    else:
                        new_coords[i] = self.coefs_[loc_index][0] * current_vars[self.alloc_vars_[0]] + \
                                        self.coefs_[loc_index][1]
                        loc_index += 1

        self.x_ = new_coords[0]
        self.y_ = new_coords[1]
        self.z_ = new_coords[2]

    def add_neighbours(self, list_neighbours):
        for pot_neighbour in list_neighbours:
            if pot_neighbour not in self.neighbours_:
                self.neighbours_ = self.neighbours_.union(set([pot_neighbour]))

    def force_coord(self, new_coord, axe):
        if axe == 'x':
            self.x_ = new_coord
        elif axe == 'y':
            self.y_ == new_coord
        else:
            self.z_ = new_coord

    def get_x(self):
        return self.x_

    def get_y(self):
        return self.y_

    def get_z(self):
        return self.z_

    def get_index(self):
        return self.global_index_

    def get_bounds(self):
        return self.bound_

    def get_bounds_type(self):
        return self.bound_types_

    def get_deg_freedom(self):
        return self.deg_freedom_

    def get_freedom(self):
        return self.freedom_

    def get_coefs(self):
        return self.coefs_

    def get_alloc_vars(self):
        return self.alloc_vars_

    def get_neighbours(self):
        return self.neighbours_


def build_node_objects(point_list, vars_list, dict_bounds):
    "Build a list of node objects from the list of points"
    nb_points = len(point_list)
    node_list = [0] * nb_points
    for i in xrange(nb_points):
        current_bounds = []
        current_bounds_type = []
        if 'plane' in dict_bounds:
            for j,bound in enumerate(dict_bounds['plane']):
                current_eq = abs(
                    point_list[i][0] * bound[0][0] + point_list[i][1] * bound[0][1] + point_list[i][2] * bound[0][2] + bound[0][3])
                if current_eq < epsilon:
                    if internal_bound:
                        if bound[1][0][1] <= numpy.dot(bound[1][0][0],point_list[i]) <= bound[1][0][2] and bound[1][1][1] <= numpy.dot(bound[1][1][0],point_list[i]) <= bound[1][1][2]:
                            current_bounds.append(j)
                            current_bounds_type.append('plane')
                    else:
                        current_bounds.append(j)
                        current_bounds_type.append('plane')
        if 'complex' in dict_bounds:
            for j,topo in enumerate(dict_bounds['complex']):
                loc_index = [ind for ind in xrange(3) if ind != topo[2]]
                current_eq = abs(point_list[i][topo[2]]-scipy.interpolate.bisplev(point_list[i][loc_index[0]],point_list[i][loc_index[1]],topo[0]))
                if current_eq < 0.001:
                    current_bounds.append(j)
                    current_bounds_type.append('complex')
        node_list[i] = Node(i, point_list[i], current_bounds, current_bounds_type, vars_list, dict_bounds)
    return node_list


def compute_complex_intersect(node_list, dict_bounds):
    #TODO tackle case where bound_types is ['complex','complex']
    intersects = {}
    out = {}
    for i,node in enumerate(node_list):

        current_deg = node.get_deg_freedom()
        if current_deg == 1:
            current_bounds_type = node.get_bounds_type()
            if current_bounds_type[0] is 'plane' and current_bounds_type[1] is 'complex':
                current_bounds = node.get_bounds()
                if tuple(current_bounds) in intersects:
                    intersects[tuple(current_bounds)] += [i]
                else:
                    intersects[tuple(current_bounds)] = [i]
        # elif current_deg == 0:
        #     current_bounds_type = node.get_bounds_type()
        #     # print i, current_bounds_type
        #     if current_bounds_type[-1] is 'complex':
        #         for t in xrange(len(current_bounds_type)-1):
        #             if current_bounds_type[t] is 'plane':
        #                 current_bounds = node.get_bounds()
        #                 cur_tuple = (current_bounds[t],current_bounds[-1])
        #                 if cur_tuple in intersects:
        #                     intersects[cur_tuple] += [i]
        #                 else:
        #                     intersects[cur_tuple] = [i]

    for key in intersects:

        current_freedom = None
        tick = 0
        while current_freedom is None:
            current_freedom = node_list[intersects[key][tick]].get_freedom()
            tick +=1
        free = 0
        while current_freedom[free] == 0:
            free += 1
        l = dict_bounds['complex'][key[1]][2]
        k = 0
        while dict_bounds['plane'][key[0]][0][k] == 0 and dict_bounds['plane'][key[0]][0][k] != l:
            k += 1

        test = False
        for i_test in xrange(3):
            if dict_bounds['plane'][key[0]][0][i_test] != 0 and i_test != k:
                test = True
        if test:
            w0 = dict_bounds['plane'][key[0]][0][k] # + dict_bounds['plane'][key[0]][3]
            tet = acos(1 / w0)
        else:
            tet = numpy.pi / 2

        if free == 0:
            if k == 1:
                data_u = numpy.array([node_list[index].get_x() for index in intersects[key]])
                data_w = numpy.array([node_list[index].get_y() * cos(tet) + node_list[index].get_z() * sin(tet) for index in intersects[key]])
            else:
                data_u = numpy.array([node_list[index].get_x() for index in intersects[key]])
                data_w = numpy.array([node_list[index].get_z() * cos(tet) + node_list[index].get_y() * sin(tet) for index in intersects[key]])
        elif free == 1:
            if k == 0:
                data_u = numpy.array([node_list[index].get_y() for index in intersects[key]])
                data_w = numpy.array([node_list[index].get_x() * cos(tet) + node_list[index].get_z() * sin(tet) for index in intersects[key]])
            else:
                data_u = numpy.array([node_list[index].get_y() for index in intersects[key]])
                data_w = numpy.array([node_list[index].get_z() * cos(tet) + node_list[index].get_x() * sin(tet) for index in intersects[key]])
        else:
            if k == 0:
                data_u = numpy.array([node_list[index].get_z() for index in intersects[key]])
                data_w = numpy.array([node_list[index].get_x() * cos(tet) + node_list[index].get_y() * sin(tet) for index in intersects[key]])
            else:
                data_u = numpy.array([node_list[index].get_z() for index in intersects[key]])
                data_w = numpy.array([node_list[index].get_y() * cos(tet) + node_list[index].get_x() * sin(tet) for index in intersects[key]])

        temp_dict = {data_u[i]:data_w[i] for i in xrange(len(data_u))}
        data_u.sort()
        for i,coord in enumerate(data_u):
            data_w[i] = temp_dict[coord]

        inter_spline = scipy.interpolate.splrep(data_u, data_w)
        # U = numpy.linspace(data_u[0], data_u[-1], 100)
        # W = numpy.array([scipy.interpolate.splev(U[i],inter_spline) for i in xrange(len(U))])
        # if free == 0:
        #     V = [node_list[intersects[key][0]].get_y() for tg in xrange(100)]
        #     scat.scatter(U,V,W)
        # else:
        #     V= [node_list[intersects[key][0]].get_x() for tg in xrange(100)]
        #     scat.scatter(V,U,W)


        a, b, c, d = dict_bounds['plane'][key[0]][0][free], dict_bounds['plane'][key[0]][0][k], dict_bounds['plane'][key[0]][0][l], dict_bounds['plane'][key[0]][0][3]

        k_coefs = [-a/(b*(1-c*cos(tet)/(b*sin(tet)))), -c/(b*sin(tet)), -d/(b*(1-c*cos(tet)/(b*sin(tet))))]
        l_coefs = [a*cos(tet)/(sin(tet)*b*(1-c*cos(tet)/(b*sin(tet)))), (1+c*cos(tet)/(b*sin(tet)))/sin(tet), d*cos(tet)/(sin(tet)*b*(1-c*cos(tet)/(b*sin(tet))))]

        # if l == k:
        #     if sin(tet) != 0:
        #         k_coefs = [-a/(b*(1-c*cos(tet)/(b*sin(tet)))), -c/(b*sin(tet)), -d/(b*(1-c*cos(tet)/(b*sin(tet))))]
        #         m_coefs = [a*cos(tet)/(sin(tet)*b*(1-c*cos(tet)/(b*sin(tet)))), (1+c*cos(tet)/(b*sin(tet)))/sin(tet), d*cos(tet)/(sin(tet)*b*(1-c*cos(tet)/(b*sin(tet))))]
        #     else:
        #         k_coefs = [0,1,0]
        #         m_coefs = [0,0,0]
        # else:
        #     if cos(tet) != 0:
        #         m_coefs = [-a / (b * (1 - c * sin(tet) / (b * cos(tet)))), -c / (b * cos(tet)),
        #                    -d / (b * (1 - c * sin(tet) / (b * cos(tet))))]
        #         k_coefs = [a * sin(tet) / (cos(tet) * b * (1 - c * sin(tet) / (b * cos(tet)))),
        #                    (1 + c * sin(tet) / (b * cos(tet))) / cos(tet),
        #                    d * sin(tet) / (cos(tet) * b * (1 - c * sin(tet) / (b * cos(tet))))]
        #     else:
        #         m_coefs = [0, 1, 0]
        #         k_coefs = [0, 0, 0]
            # if cos(tet) != 0:
            #     m_coefs = [-a / (c * (1 - b * sin(tet) / (c * cos(tet)))), -b / (c * cos(tet)), -d / (c * (1 - b * sin(tet) / (c * cos(tet))))]
            #     k_coefs = [a * sin(tet) / (cos(tet) * c * (1 - b * sin(tet) / (c * cos(tet)))), (1 + b * sin(tet) / (c * cos(tet))) / cos(tet), d * sin(tet) / (cos(tet) * c * (1 - b * sin(tet) / (c * cos(tet))))]
            # else:
            #     m_coefs = [0,1,0]
            #     k_coefs = [0,0,0]
        if k < l:
            coefs = [k_coefs, l_coefs]
        else:
            coefs = [l_coefs, k_coefs]
        coefs.append(inter_spline)
        out[key] = coefs

    return out


def update_complex_intersect(node_list, out):

    for node in node_list:
        current_deg = node.get_deg_freedom()
        if current_deg == 1:
            current_bounds_type = node.get_bounds_type()
            if current_bounds_type[0] is 'plane' and current_bounds_type[1] is 'complex':
                current_bounds = tuple(node.get_bounds())
                for key in out:
                    if all(key[i] == current_bounds[i] for i in xrange(2)):
                        node.update_coefs(out[key])
    return node_list


def update_node_coords(node_list, var_coords, dict_bounds):
    "Update the coords of the node objects that are being optimized"

    for node in node_list:
        if node.get_freedom() is not None:
            node.update_coords(var_coords, dict_bounds)
    return node_list


def format_results(node_list):
    "Formats the results of the optimization into a list of new coords + corresponding index of the point"

    lenght = len(node_list)
    data = [[] for i in xrange(lenght)]
    counter = 0
    for i, node in enumerate(node_list):
        current_freedom = node.get_freedom()
        if current_freedom is not None:
            data[counter] = [node.get_x(), node.get_y(), node.get_z(), i]
            counter += 1
    data = data[:counter]
    return data


def build_neighbours(node_list, con_obj):
    "construct the list of neighbours of every node that is going to be optimized"

    faces = con_obj.get_faces()
    vertices = con_obj.get_vertices()
    for i,face in enumerate(faces):
        current_nodes = face + vertices[i]
        for node in current_nodes:
            node_list[node].add_neighbours([pot_neighbour for pot_neighbour in current_nodes if pot_neighbour is not node])
    return node_list


# Data extracting tools


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


# Quality Control


def quality_control(node_list, con_obj):
    "Tool used to estimate the quality of a mesh in term of orthogonality"

    faces = con_obj.get_faces()
    vertices = con_obj.get_vertices()
    max_lenght = len(faces)
    orth_deviation_distribution = numpy.zeros(max_lenght)
    counter = 0
    for i, face in enumerate(faces):
        vertice = vertices[i]
        if len(vertice) == 2:
            x0, y0, z0 = node_list[face[0]].get_x(), node_list[face[0]].get_y(), node_list[face[0]].get_z()
            x1, y1, z1 = node_list[face[1]].get_x(), node_list[face[1]].get_y(), node_list[face[1]].get_z()
            x2, y2, z2 = node_list[face[2]].get_x(), node_list[face[2]].get_y(), node_list[face[2]].get_z()
            x4, y4, z4 = node_list[vertice[0]].get_x(), node_list[vertice[0]].get_y(), node_list[vertice[0]].get_z()
            x5, y5, z5 = node_list[vertice[1]].get_x(), node_list[vertice[1]].get_y(), node_list[vertice[1]].get_z()
            dot_product = ((x5 - x4) * ((y0 - y1) * (z2 - z1) - (z0 - z1) * (y2 - y1)) - (y5 - y4) * (
                        (x0 - x1) * (z2 - z1) - (z0 - z1) * (x2 - x1)) + (z5 - z4) * (
                                       (x0 - x1) * (y2 - y1) - (y0 - y1) * (x2 - x1))) / (
                                      ((((x5 - x4) ** 2) + (y5 - y4) ** 2 + (z5 - z4) ** 2) ** (1. / 2)) * (
                                          ((y0 - y1) * (z2 - z1) - (z0 - z1) * (y2 - y1)) ** 2 + (
                                      -((x0 - x1) * (z2 - z1) - (z0 - z1) * (x2 - x1))) ** 2 + (
                                                      (x0 - x1) * (y2 - y1) - (y0 - y1) * (x2 - x1)) ** 2) ** (1. / 2))

            orth_deviation_distribution[counter] = 1 - abs(dot_product)
            counter += 1
    orth_deviation_distribution.resize(counter)
    return orth_deviation_distribution


def finding_targets(nodes_list,con_obj,limit, n_process):

    faces = con_obj.get_faces()
    lenght = len(faces)
    sub_size = int(lenght / nb_process) + 1
    vertices = con_obj.get_vertices()
    faces_tasks = numpy.array([faces[i:i + sub_size] for i in xrange(0, lenght, sub_size)])
    vertices_tasks = numpy.array([vertices[i:i + sub_size] for i in xrange(0, lenght, sub_size)])
    tasks = [[face, vertices_tasks[iv], nodes_list,limit] for iv, face in enumerate(faces_tasks)]
    output_tasks = Queue()
    workers = [Process(target=exctract_bad_elements, args=(task, output_tasks)) for task in tasks]
    for p in workers:
        p.start()
    # for p in processes:
    #     p.terminate()
    # for p in processes:
    #     p.join()
    bad_nodes = set([])
    for ind in xrange(len(tasks)):
        n_set = output_tasks.get()
        bad_nodes = bad_nodes | n_set

    return numpy.array(list(bad_nodes))


def exctract_bad_elements(t_in, t_out):
    "Return a list of nodes index that play a role in deviations higher than the given limit"

    faces, vertices, node_list, limit = t_in
    temp_bn = []
    # faces = con_obj.get_faces()
    # vertices = con_obj.get_vertices()
    # max_lenght = len(node_list)
    # bad_nodes = -1 * numpy.ones(max_lenght)
    # counter = 0
    for i, face in enumerate(faces):

        vertice = vertices[i]
        if len(vertice) == 2:
            x0, y0, z0 = node_list[face[0]].get_x(), node_list[face[0]].get_y(), node_list[face[0]].get_z()
            x1, y1, z1 = node_list[face[1]].get_x(), node_list[face[1]].get_y(), node_list[face[1]].get_z()
            x2, y2, z2 = node_list[face[2]].get_x(), node_list[face[2]].get_y(), node_list[face[2]].get_z()
            x4, y4, z4 = node_list[vertice[0]].get_x(), node_list[vertice[0]].get_y(), node_list[vertice[0]].get_z()
            x5, y5, z5 = node_list[vertice[1]].get_x(), node_list[vertice[1]].get_y(), node_list[vertice[1]].get_z()
            dot_product = ((x5 - x4) * ((y0 - y1) * (z2 - z1) - (z0 - z1) * (y2 - y1)) - (y5 - y4) * (
                    (x0 - x1) * (z2 - z1) - (z0 - z1) * (x2 - x1)) + (z5 - z4) * (
                                   (x0 - x1) * (y2 - y1) - (y0 - y1) * (x2 - x1))) / (
                                  ((((x5 - x4) ** 2) + (y5 - y4) ** 2 + (z5 - z4) ** 2) ** (1. / 2)) * (
                                  ((y0 - y1) * (z2 - z1) - (z0 - z1) * (y2 - y1)) ** 2 + (
                              -((x0 - x1) * (z2 - z1) - (z0 - z1) * (x2 - x1))) ** 2 + (
                                          (x0 - x1) * (y2 - y1) - (y0 - y1) * (x2 - x1)) ** 2) ** (1. / 2))

            dev = 1 - abs(dot_product)
            if dev > limit:
                for j in xrange(3):
                    temp_bn.append(face[j])
                    # if face[j] not in bad_nodes:
                    #     bad_nodes[counter] = face[j]
                    #     counter += 1
                for k in xrange(2):
                    temp_bn.append(vertice[k])
                    # if vertice[k] not in bad_nodes:
                    #     bad_nodes[counter] = vertice[k]
                    #     counter += 1
    t_out.put(set(temp_bn))
    # bad_nodes.resize(counter)
    # return bad_nodes


# Optimization


def initial_state(node_list):
    "Initial coordinate of the points to feed to the optimizing function"
    nb_points = len(node_list)
    init_state = numpy.zeros(3 * nb_points)
    last_index = 0
    for node in node_list:
        current_freedom = node.get_freedom()
        if current_freedom is not None:
            node.allocation_var_indexes(last_index)
            potential_import = [node.get_x(), node.get_y(), node.get_z()]
            for i in xrange(3):
                if current_freedom[i] == 1:
                    init_state[last_index] = potential_import[i]
                    last_index += 1
    init_state.resize(last_index)
    return init_state, node_list


def objective_function(N, *args):
    "Objective function: p-norm of the list of deviations"
    # function we want to minimize
    # N = list of variables
    # points = list of initial point coords
    # iteration = dictionary of connections
    # freedom = dictionary of degree of freedom of points
    # boundaries = dictionary of plane equations
    nodes, connect, dict_bounds = args
    global update_status
    update_status = False
    if update_status is False:
        nodes = update_node_coords(nodes, N, dict_bounds)
        update_status = True

    faces = connect.get_faces()
    vertices = connect.get_vertices()
    orth_deviation_distribution = 0
    for i, con in enumerate(faces):
        vertice = vertices[i]
        if len(vertices[i]) == 2:

            x0, y0, z0 = nodes[con[0]].get_x(), nodes[con[0]].get_y(), nodes[con[0]].get_z()
            x1, y1, z1 = nodes[con[1]].get_x(), nodes[con[1]].get_y(), nodes[con[1]].get_z()
            x2, y2, z2 = nodes[con[2]].get_x(), nodes[con[2]].get_y(), nodes[con[2]].get_z()
            x4, y4, z4 = nodes[vertice[0]].get_x(), nodes[vertice[0]].get_y(), nodes[vertice[0]].get_z()
            x5, y5, z5 = nodes[vertice[1]].get_x(), nodes[vertice[1]].get_y(), nodes[vertice[1]].get_z()
            dot_product = ((x5 - x4) * ((y0 - y1) * (z2 - z1) - (z0 - z1) * (y2 - y1)) - (y5 - y4) * (
                        (x0 - x1) * (z2 - z1) - (z0 - z1) * (x2 - x1)) + (z5 - z4) * (
                                       (x0 - x1) * (y2 - y1) - (y0 - y1) * (x2 - x1))) / (
                                      ((((x5 - x4) ** 2) + (y5 - y4) ** 2 + (z5 - z4) ** 2) ** (1. / 2)) * (
                                          ((y0 - y1) * (z2 - z1) - (z0 - z1) * (y2 - y1)) ** 2 + (
                                      -((x0 - x1) * (z2 - z1) - (z0 - z1) * (x2 - x1))) ** 2 + (
                                                      (x0 - x1) * (y2 - y1) - (y0 - y1) * (x2 - x1)) ** 2) ** (1. / 2))
            orth_deviation_distribution += 1 - abs(dot_product)

    return orth_deviation_distribution


def build_jacobian(N, *args):
    "Build the analytical partial derivatives of the objective function wrt to each variables"
    # Matrix with one 1 lign and len(N) columns with the derivatives of obj function wrt each variables
    # N = list of variables, 3*number_of_points variables
    # sequence = dictionary obtained with connection_mapping method
    # vertices_indexes = dictionary obtained with connection_mapping_opti method
    # elements_indexes = list of vertices indexes for every tet in the mesh
    nodes, connect, dict_bounds = args
    global update_status
    if update_status is False:
        nodes = update_node_coords(nodes, N, dict_bounds)
        update_status = True

    faces = connect.get_faces()
    vertices = connect.get_vertices()
    max_vars = len(N)
    jac = numpy.zeros(max_vars)
    for i,con in enumerate(faces):
        vertice = vertices[i]
        if len(vertice) == 2:
            x0, y0, z0 = nodes[con[0]].get_x(), nodes[con[0]].get_y(), nodes[con[0]].get_z()
            x1, y1, z1 = nodes[con[1]].get_x(), nodes[con[1]].get_y(), nodes[con[1]].get_z()
            x2, y2, z2 = nodes[con[2]].get_x(), nodes[con[2]].get_y(), nodes[con[2]].get_z()
            x4, y4, z4 = nodes[vertice[0]].get_x(), nodes[vertice[0]].get_y(), nodes[vertice[0]].get_z()
            x5, y5, z5 = nodes[vertice[1]].get_x(), nodes[vertice[1]].get_y(), nodes[vertice[1]].get_z()
            #Computing of local coeficient
            dot_product = ((x5 - x4) * ((y0 - y1) * (z2 - z1) - (z0 - z1) * (y2 - y1)) - (y5 - y4) * (
                    (x0 - x1) * (z2 - z1) - (z0 - z1) * (x2 - x1)) + (z5 - z4) * (
                                   (x0 - x1) * (y2 - y1) - (y0 - y1) * (x2 - x1))) / (
                                  ((((x5 - x4) ** 2) + (y5 - y4) ** 2 + (z5 - z4) ** 2) ** (1. / 2)) * (
                                  ((y0 - y1) * (z2 - z1) - (z0 - z1) * (y2 - y1)) ** 2 + (
                              -((x0 - x1) * (z2 - z1) - (z0 - z1) * (x2 - x1))) ** 2 + (
                                          (x0 - x1) * (y2 - y1) - (y0 - y1) * (x2 - x1)) ** 2) ** (1. / 2))
            sign = 1.0
            if dot_product > 0:
                sign = -1.0
            #Get Matrix of all partial derivatives (3*5)
            matrix_ders = numpy.zeros((5,3))
            if nodes[con[0]].get_freedom() is not None:
                matrix_ders[0][0] = (-0.5 * (-2 * y1 + 2 * y2) * ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) - 0.5 * (
                                        2 * z1 - 2 * z2) * (-(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1))) * (
                                             (-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                                                 (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                         (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                             (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                             ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                                 -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                         (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-1.5) + (
                                             (-y1 + y2) * (-z4 + z5) - (-y4 + y5) * (-z1 + z2)) * (
                                             (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                             ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                                 -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                         (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5)
                matrix_ders[0][1] = ((x1 - x2) * (-z4 + z5) + (-x4 + x5) * (-z1 + z2)) * (
                        (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5) + (
                                            -0.5 * (2 * x1 - 2 * x2) * ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) - 0.5 * (
                                            -2 * z1 + 2 * z2) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1))) * (
                                            (-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                                            (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                    (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-1.5)
                matrix_ders[0][2] = (-0.5 * (-2 * x1 + 2 * x2) * (-(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) - 0.5 * (
                        2 * y1 - 2 * y2) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1))) * (
                                            (-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                                            (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                    (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-1.5) + (
                                            -(x1 - x2) * (-y4 + y5) + (-x4 + x5) * (y1 - y2)) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5)
            if nodes[con[1]].get_freedom() is not None:
                matrix_ders[1][0] = ((y0 - y2) * (-z4 + z5) - (-y4 + y5) * (z0 - z2)) * (
                        (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5) + (
                                            -0.5 * (2 * y0 - 2 * y2) * ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) - 0.5 * (
                                            -2 * z0 + 2 * z2) * (-(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1))) * (
                                            (-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                                            (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                    (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-1.5)
                matrix_ders[1][1] = (-0.5 * (-2 * x0 + 2 * x2) * ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) - 0.5 * (
                        2 * z0 - 2 * z2) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1))) * (
                                            (-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                                            (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                    (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-1.5) + (
                                            (-x0 + x2) * (-z4 + z5) + (-x4 + x5) * (z0 - z2)) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5)
                matrix_ders[1][2] = (-(-x0 + x2) * (-y4 + y5) + (-x4 + x5) * (-y0 + y2)) * (
                        (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5) + (
                                            -0.5 * (2 * x0 - 2 * x2) * (-(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) - 0.5 * (
                                            -2 * y0 + 2 * y2) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1))) * (
                                            (-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                                            (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                    (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-1.5)
            if nodes[con[2]].get_freedom() is not None:
                matrix_ders[2][0] = (-0.5 * (-2 * y0 + 2 * y1) * ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) - 0.5 * (
                        2 * z0 - 2 * z1) * (-(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1))) * (
                                            (-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                                            (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                    (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-1.5) + (
                                            (-y0 + y1) * (-z4 + z5) - (-y4 + y5) * (-z0 + z1)) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5)
                matrix_ders[2][1] = ((x0 - x1) * (-z4 + z5) + (-x4 + x5) * (-z0 + z1)) * (
                        (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5) + (
                                            -0.5 * (2 * x0 - 2 * x1) * ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) - 0.5 * (
                                            -2 * z0 + 2 * z1) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1))) * (
                                            (-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                                            (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                    (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-1.5)
                matrix_ders[2][2] = (-0.5 * (-2 * x0 + 2 * x1) * (-(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) - 0.5 * (
                        2 * y0 - 2 * y1) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1))) * (
                                            (-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                                            (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                    (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-1.5) + (
                                            -(x0 - x1) * (-y4 + y5) + (-x4 + x5) * (y0 - y1)) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5)
            if nodes[vertice[0]].get_freedom() is not None:
                matrix_ders[3][0] = (-1.0 * x4 + 1.0 * x5) * ((-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                        (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                                      (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-1.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5) + (
                                            -(y0 - y1) * (-z1 + z2) + (-y1 + y2) * (z0 - z1)) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5)
                matrix_ders[3][1] = (-1.0 * y4 + 1.0 * y5) * ((-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                        (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                                      (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-1.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5) + (
                                            (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5)
                matrix_ders[3][2] = (-1.0 * z4 + 1.0 * z5) * ((-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                        (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                                      (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-1.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5) + (
                                            -(x0 - x1) * (-y1 + y2) + (-x1 + x2) * (y0 - y1)) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5)
            if nodes[vertice[1]].get_freedom() is not None:
                matrix_ders[4][0] = (1.0 * x4 - 1.0 * x5) * ((-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                        (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                                     (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-1.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5) + (
                                            (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5)
                matrix_ders[4][1] = (1.0 * y4 - 1.0 * y5) * ((-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                        (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                                     (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-1.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5) + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5)
                matrix_ders[4][2] = (1.0 * z4 - 1.0 * z5) * ((-x4 + x5) * ((y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) - (-y4 + y5) * (
                        (x0 - x1) * (-z1 + z2) - (-x1 + x2) * (z0 - z1)) + (-z4 + z5) * (
                                                                     (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1))) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-1.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5) + (
                                            (x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) * (
                                            (-x4 + x5) ** 2 + (-y4 + y5) ** 2 + (-z4 + z5) ** 2) ** (-0.5) * (
                                            ((x0 - x1) * (-y1 + y2) - (-x1 + x2) * (y0 - y1)) ** 2 + (
                                            -(x0 - x1) * (-z1 + z2) + (-x1 + x2) * (z0 - z1)) ** 2 + (
                                                    (y0 - y1) * (-z1 + z2) - (-y1 + y2) * (z0 - z1)) ** 2) ** (-0.5)

            current_points = [con[0], con[1], con[2], vertice[0], vertice[1]]
            for i,point in enumerate(current_points):
                current_freedom = nodes[point].get_freedom()
                if current_freedom is not None:
                    current_deg = nodes[point].get_deg_freedom()
                    current_coefs = nodes[point].get_coefs()
                    current_ders = [0., 0., 0.]
                    local_der_x = sign * matrix_ders[i][0]
                    local_der_y = sign * matrix_ders[i][1]
                    local_der_z = sign * matrix_ders[i][2]
                    if current_deg == 3:
                        current_ders[0] += local_der_x
                        current_ders[1] += local_der_y
                        current_ders[2] += local_der_z

                    elif current_deg == 2:
                        current_bound_types = nodes[point].get_bounds_type()
                        if 'complex' in current_bound_types:
                            current_bounds = nodes[point].get_bounds()
                            if dict_bounds['complex'][current_bounds[0]][2] == 0:
                                current_ders[1] += local_der_y + scipy.interpolate.bisplev(nodes[point].get_y(),
                                                                                           nodes[point].get_z(),
                                                                                           dict_bounds['complex'][
                                                                                               current_bounds[0]][0],
                                                                                           dx=1, dy=0) * local_der_x
                                current_ders[2] += local_der_z + scipy.interpolate.bisplev(nodes[point].get_y(),
                                                                                           nodes[point].get_z(),
                                                                                           dict_bounds['complex'][
                                                                                               current_bounds[0]][0],
                                                                                           dx=0, dy=1) * local_der_x
                            elif dict_bounds['complex'][current_bounds[0]][2] == 1:
                                current_ders[0] += local_der_x + scipy.interpolate.bisplev(nodes[point].get_x(),
                                                                                           nodes[point].get_z(),
                                                                                           dict_bounds['complex'][
                                                                                               current_bounds[0]][0],
                                                                                           dx=1, dy=0) * local_der_y
                                current_ders[2] += local_der_z + scipy.interpolate.bisplev(nodes[point].get_x(),
                                                                                           nodes[point].get_z(),
                                                                                           dict_bounds['complex'][
                                                                                               current_bounds[0]][0],
                                                                                           dx=0, dy=1) * local_der_y
                            else:
                                current_ders[0] += local_der_x + scipy.interpolate.bisplev(nodes[point].get_x(),
                                                                                           nodes[point].get_y(),
                                                                                           dict_bounds['complex'][
                                                                                               current_bounds[0]][0],
                                                                                           dx=1, dy=0) * local_der_z
                                current_ders[1] += local_der_y + scipy.interpolate.bisplev(nodes[point].get_x(),
                                                                                           nodes[point].get_z(),
                                                                                           dict_bounds['complex'][
                                                                                               current_bounds[0]][0],
                                                                                           dx=0, dy=1) * local_der_z
                        else:
                            if current_freedom[2] == 0:
                                current_ders[0] += local_der_x + current_coefs[0] * local_der_z
                                current_ders[1] += local_der_y + current_coefs[1] * local_der_z
                            elif current_freedom[1] == 0:
                                current_ders[0] += local_der_x + current_coefs[0] * local_der_y
                                current_ders[2] += local_der_z + current_coefs[1] * local_der_y
                            else:
                                current_ders[1] += local_der_y + current_coefs[0] * local_der_x
                                current_ders[2] += local_der_z + current_coefs[1] * local_der_x

                    else:
                        current_bound_types = nodes[point].get_bounds_type()
                        if 'complex' in current_bound_types:
                            current_coefs = nodes[point].get_coefs()
                            if 'plane' in current_bound_types:
                                #TODO modify the derivatives, not good projection  (Normally done)
                                if current_freedom[0] == 1:
                                    current_ders[0] += local_der_x + (current_coefs[0][0]-current_coefs[0][1]*scipy.interpolate.splev(nodes[point].get_x(),current_coefs[2], der=1)) * local_der_y + (current_coefs[1][0]-current_coefs[1][1]*scipy.interpolate.splev(nodes[point].get_x(),current_coefs[2], der=1)) * local_der_z
                                elif current_freedom[1] == 1:
                                    current_ders[1] += (current_coefs[0][0]-current_coefs[0][1]*scipy.interpolate.splev(nodes[point].get_y(),current_coefs[2], der=1)) * local_der_x + local_der_y + (current_coefs[1][0]-current_coefs[1][1]*scipy.interpolate.splev(nodes[point].get_y(),current_coefs[2], der=1)) * local_der_z
                                else:
                                    current_ders[2] += (current_coefs[0][0]-current_coefs[0][1]*scipy.interpolate.splev(nodes[point].get_z(),current_coefs[2], der=1)) * local_der_x + (current_coefs[1][0]-current_coefs[1][1]*scipy.interpolate.splev(nodes[point].get_z(),current_coefs[2], der=1)) * local_der_y + local_der_z
                            else:
                                #TODO Modify the derivatives, not good projection
                                if current_freedom[0] == 1:
                                    if current_bounds[0][2] == 1:
                                        loc_coefs = [scipy.interpolate.bisplev(nodes[point].get_x(), nodes[point].get_z(),current_bounds[0][0],dx=1, dy=0), scipy.interpolate.bisplev(nodes[point].get_x(), nodes[point].get_y(),current_bounds[1][0],dx=1, dy=0)]
                                    else:
                                        loc_coefs = [scipy.interpolate.bisplev(nodes[point].get_x(), nodes[point].get_z(), current_bounds[1][0], dx=1, dy=0), scipy.interpolate.bisplev(nodes[point].get_x(), nodes[point].get_y(), current_bounds[0][0], dx=1, dy=0)]
                                    current_ders[0] += local_der_x + loc_coefs[0] * local_der_y + loc_coefs[1] * local_der_z
                                elif current_freedom[1] == 1:
                                    if current_bounds[0][2] == 0:
                                        loc_coefs = [scipy.interpolate.bisplev(nodes[point].get_y(), nodes[point].get_z(),current_bounds[0][0],dx=1, dy=0), scipy.interpolate.bisplev(nodes[point].get_x(), nodes[point].get_y(),current_bounds[1][0],dx=0, dy=1)]
                                    else:
                                        loc_coefs = [scipy.interpolate.bisplev(nodes[point].get_y(), nodes[point].get_z(), current_bounds[1][0], dx=1, dy=0), scipy.interpolate.bisplev(nodes[point].get_x(), nodes[point].get_y(), current_bounds[0][0], dx=0, dy=1)]
                                    current_ders[0] += loc_coefs[0] * local_der_x + local_der_y + loc_coefs[1] * local_der_z
                                else:
                                    if current_bounds[0][2] == 0:
                                        loc_coefs = [scipy.interpolate.bisplev(nodes[point].get_x(), nodes[point].get_z(),current_bounds[0][0],dx=0, dy=1), scipy.interpolate.bisplev(nodes[point].get_y(), nodes[point].get_z(),current_bounds[1][0],dx=0, dy=1)]
                                    else:
                                        loc_coefs = [scipy.interpolate.bisplev(nodes[point].get_x(), nodes[point].get_z(), current_bounds[1][0], dx=0, dy=1), scipy.interpolate.bisplev(nodes[point].get_y(), nodes[point].get_z(), current_bounds[0][0], dx=0, dy=1)]
                                    current_ders[0] += loc_coefs[0] * local_der_x + local_der_y + loc_coefs[1] * local_der_z
                        else:
                            if current_freedom[0] == 1:
                                current_ders[0] += local_der_x + current_coefs[0][0] * local_der_y + current_coefs[1][
                                    0] * local_der_z
                            elif current_freedom[1] == 1:
                                current_ders[1] += local_der_y + current_coefs[0][0] * local_der_x + current_coefs[1][
                                    0] * local_der_z
                            else:
                                current_ders[2] += local_der_z + current_coefs[0][0] * local_der_x + current_coefs[1][
                                    0] * local_der_y
                    current_alloc = nodes[point].get_alloc_vars()
                    counter = 0
                    for i in xrange(3):
                        if current_freedom[i] == 1:
                            jac[current_alloc[counter]] += current_ders[i]
                            counter += 1

    return jac


def run_optimization(data, output):
    "Algorithm that creates the necessary objects and then launches the optimization, the output is added to the current queue of results"

    vars_indexes, dict_con, el_con, el_list, bound_list, points_list, maxiter, ftol, out = data
    cons = Con(dict_con, el_con, el_list,vars_indexes)
    nodes_list = build_node_objects(points_list,vars_indexes, bound_list)
    nodes_list = update_complex_intersect(nodes_list, out)
    source_vars, nodes_list = initial_state(nodes_list)
    global update_status
    update_status = True

    result = minimize(objective_function, source_vars,
                                     args=(nodes_list,cons, bound_list), method='SLSQP',
                                     jac=build_jacobian,
                                     options={'disp': False, 'ftol': ftol, 'maxiter': maxiter})

    nodes_list = update_node_coords(nodes_list, result.x, bound_list)
    new_points_list = format_results(nodes_list)
    if multiprocess:
        output.put(new_points_list)
    else:
        return [new_points_list]


# Input/Output Manipulation


def compile_results(point_list, results):
    "Rebuilds the list of points using the results from the previous chunks optimizations"

    for result in results:
        for point in result:
            point_list[point[3]][0] = point[0]
            point_list[point[3]][1] = point[1]
            point_list[point[3]][2] = point[2]
    return point_list


def pack_inputs(chunks, dict_con, el_con, el_list, bound_list, points_list, maxiter, ftol, out):
    "Compress the args needed in run_optimization into a list so that run_optimization can be called using only 1 arg"

    inputs = []
    for chunk in chunks:
        inputs.append([chunk, dict_con, el_con, el_list, bound_list, points_list, maxiter, ftol, out])
    return inputs


def calibrate_tol(el_list, points_list, chunks, dict_bounds, out):
    "Find the smallest tolerance in scipy.optimize.minimize so that the system converge in more that 1 iteration"
    max_ftol = -1
    for vars_indexes in chunks:
        ftol = 10
        nb_it = 1
        dict_con, el_con = connection_mapping(el_list)
        cons = Con(dict_con, el_con, el_list, vars_indexes)
        nodes_list = build_node_objects(points_list, vars_indexes, dict_bounds)
        nodes_list = update_complex_intersect(nodes_list, out)
        source_vars, nodes_list = initial_state(nodes_list)
        global update_status
        update_status = True
        while nb_it == 1:
            result = minimize(objective_function, source_vars,
                              args=(nodes_list, cons, dict_bounds), method='SLSQP',
                              jac=build_jacobian,
                              options={'disp': False, 'ftol': ftol, 'maxiter': 1})
            nb_it = result.nit
            ftol = ftol * 0.5
        if ftol > max_ftol:
            max_ftol = ftol

    return max_ftol


def unpack_bound_data(complex_data, plane_data):
    if complex_data is not None:
        bound_data = [complex_data]
    else:
        bound_data = []
    bulkplanes = numpy.loadtxt(plane_data)
    nb_it = len(bulkplanes)//4
    for it in xrange(nb_it):
        bound_data.append([bulkplanes[4*it+j] for j in xrange(4)])
    return bound_data


def format_bounds(complex_data, plane_data, points_list):
    bulk_bound_list = unpack_bound_data(complex_data, plane_data)
    form_bounds = {'plane':[],'complex':[]}
    for bound in bulk_bound_list:
        indexes = []
        if isinstance(bound[0], numpy.ndarray):
            plane_eq = get_plane_equation(bound[:-1])
            def_domain = get_vect_direction(bound)
            form_bounds['plane'].append([plane_eq, def_domain])
        else:
            if bound[0][-4:] == '.msh':
                ref_points, _, _, _, _ = meshio.read(bound[0])
                ref_plane_eq = get_plane_equation(bound[1][:-1])

                max = -1
                dir = None
                for ax, coef in enumerate(ref_plane_eq[:-1]):
                    if abs(coef) > max:
                        max = abs(coef)
                        dir = ax

                for ind, point in enumerate(ref_points):
                    current_eq = point[0]*ref_plane_eq[0] + point[1]*ref_plane_eq[1] + point[2]*ref_plane_eq[2] + ref_plane_eq[3]
                    if abs(current_eq) < epsilon:
                        indexes.append(ind)
                loc_index = [j for j in xrange(3) if j != dir]
                data_u = numpy.array([points_list[i][loc_index[0]] for i in indexes])
                data_v = numpy.array([points_list[i][loc_index[1]] for i in indexes])
                data_w = numpy.array([points_list[i][dir] for i in indexes])
                # smooth = len(data_u) + numpy.sqrt(2*len(data_u))
                fit = scipy.interpolate.bisplrep(data_u, data_v, data_w)
            else:
                ref_plane_eq = get_plane_equation(bound[1][:-1])
                max = -1
                dir = None
                for ax, coef in enumerate(ref_plane_eq[:-1]):
                    if abs(coef) > max:
                        max = abs(coef)
                        dir = ax
                loc_index = [j for j in xrange(3) if j != dir]
                bulk_data = numpy.loadtxt(bound[0])
                # x_res = len(set(list(bulk_data[:, 0])))
                # y_res = len(set(list(bulk_data[:, 1])))
                # data_u, data_v = numpy.meshgrid(numpy.linspace(0,10000,x_res),numpy.linspace(0,10000,y_res)) #42,23 --> 1000 ; 135, 73 --> 10000
                # data_u = numpy.array([point[loc_index[0]] for point in bulk_data])
                # data_v = numpy.array([point[loc_index[1]] for point in bulk_data])
                # data_w = numpy.array([point[dir] for point in bulk_data])
                fit = scipy.interpolate.bisplrep(bulk_data[:,0], bulk_data[:,1], bulk_data[:,2])

            form_bounds['complex'].append([fit, ref_plane_eq, dir])
    return form_bounds


#Chunking


def random_chunking(nb_points,sub_size):

    temp_list = numpy.arange(0,nb_points+1,1)
    numpy.random.shuffle(temp_list)
    temp_matrix = numpy.array([temp_list[i:i + sub_size] for i in xrange(0, nb_points, sub_size)])
    return temp_matrix


def chunking_bad_elements(bad_nodes, sub_size):
    nb_points = len(bad_nodes)
    numpy.random.shuffle(bad_nodes)
    temp_matrix = numpy.array([bad_nodes[i:i + sub_size] for i in xrange(0, nb_points, sub_size)])
    return temp_matrix


def chunking_disjoint(nodes_list, bad_nodes, sub_size, n):
    "Creates disjoint clusters of points so that those clusters can be optimized in parallel using multiprocessing"

    bad_nodes = list(bad_nodes)
    bulk_chunks = []
    counter = 0
    while len(bad_nodes) > 0:
        bulk_chunks += [[] for i in xrange(n)]
        infinite_loop = False
        while numpy.all(len(bulk_chunks[counter+i]) < sub_size for i in xrange(n)) and infinite_loop is False:
        # while len(bulk_chunks[counter]) < sub_size and len(bulk_chunks[counter+1]) < sub_size and len(bulk_chunks[counter+2]) < sub_size and infinite_loop is False:
            ind = 0
            infinite_loop = True
            for node in bad_nodes:
                current_neighbours = nodes_list[int(node)].get_neighbours()
                test = True
                for i in xrange(n):
                    if counter + i != counter + ind:
                        for element in bulk_chunks[counter + i]:
                            temp_neighbours = nodes_list[int(element)].get_neighbours()
                            if not temp_neighbours.isdisjoint(current_neighbours) or node in temp_neighbours:
                                test = False

                if test and len(bulk_chunks[counter+ind]) < sub_size:
                    bulk_chunks[counter+ind].append(int(node))
                    bad_nodes.remove(node)
                    infinite_loop = False

                if ind == n-1:
                    ind = 0
                else:
                    ind += 1
    #Rearrange chunk to maximize their lenght and minimize their number
        if sum([len(bulk_chunks[counter+i]) for i in xrange(n)])<= sub_size:
            for c in xrange(1,n):
                bulk_chunks[counter] += bulk_chunks[counter+c]
                bulk_chunks[counter+c] = []
        indexes = numpy.array([counter+i for i in xrange(n)])
        couples = combinations(indexes, 2)
        for couple in couples:
            if sum([len(bulk_chunks[index]) for index in couple]) < sub_size:
                bulk_chunks[couple[0]] += bulk_chunks[couple[1]]
                bulk_chunks[couple[1]] = []
        counter += n
        # if len(bulk_chunks[counter]) + len(bulk_chunks[counter+1]) + len(bulk_chunks[counter+2]) <= sub_size:
        #     bulk_chunks[counter] += bulk_chunks[counter + 1] + bulk_chunks[counter+2]
        #     bulk_chunks[counter + 1] = []
        #     bulk_chunks[counter + 2] = []
        # elif len(bulk_chunks[counter]) + len(bulk_chunks[counter+1]) <= sub_size:
        #     bulk_chunks[counter] += bulk_chunks[counter+1]
        #     bulk_chunks[counter+1] = []
        # elif len(bulk_chunks[counter]) + len(bulk_chunks[counter+2]) <= sub_size:
        #     bulk_chunks[counter] += bulk_chunks[counter+1]
        #     bulk_chunks[counter+2] = []
        # elif len(bulk_chunks[counter+1]) + len(bulk_chunks[counter+2]) <= sub_size:
        #     bulk_chunks[counter+1] += bulk_chunks[counter+2]
        #     bulk_chunks[counter+2] = []
        # counter += 3

    stop = False
    if len(bulk_chunks) < 6:
        stop = True
    while len(bulk_chunks) >= 6 and len(bulk_chunks[-3]) + len(bulk_chunks[-6]) < sub_size and stop is False:
        if len(bulk_chunks[-3]) + len(bulk_chunks[-6]) <= sub_size:
            bulk_chunks[-6] += bulk_chunks[-3]
            bulk_chunks = bulk_chunks[:-3]
        else:
            stop = True
    return bulk_chunks


def chunking_disjoint2(points, nodes_list, bad_nodes, sub_size, n):

    #Build M disjoint clusters (mathematical definition of disjoint not ours)
    n_cluster = int(len(bad_nodes)/float(sub_size))
    data = numpy.array([points[int(bn)] for bn in bad_nodes])
    clustering = KMeans(n_clusters= n_cluster, init='k-means++', n_init= 1, n_jobs= 5).fit(data)
    labels = clustering.labels_
    bulk_clusters = [[] for cn in xrange(n_cluster)]
    for ip in xrange(len(labels)):
        bulk_clusters[labels[ip]].append(int(bad_nodes[ip]))

    #Rearrange clusters in N sub groups of size K, with K= number of process, so that those sub groups can be optimised in parallel
    centroids = clustering.cluster_centers_
    parallel_clusters = []
    while len(bulk_clusters) > 0:
        temp_subgroup = [bulk_clusters[0]]
        last_neighbours = set([])
        for bc in bulk_clusters[0]:
            last_neighbours.union(nodes_list[bc].get_neighbours())
        temp_neighbours = [last_neighbours]
        dist_matrix = [numpy.dot(centroids[0],centroids[nc]) for nc in xrange(1,len(centroids))]
        bulk_clusters = bulk_clusters[1:]
        centroids = centroids[1:]
        while len(temp_subgroup) < n and len(dist_matrix) > 0:
            target = dist_matrix.index(min(dist_matrix))
            dist_matrix = dist_matrix[:target]+dist_matrix[target+1:]
            disjoint_test = True
            t_cluster = bulk_clusters[target]
            t_neighbours = set([])
            for tn in t_cluster:
                t_neighbours.union(nodes_list[tn].get_neighbours())
            for icl, cl in enumerate(temp_subgroup):
                if not t_neighbours.isdisjoint(temp_neighbours[icl]) or not t_neighbours.isdisjoint(set(cl)):
                    disjoint_test = False
            if disjoint_test:
                temp_subgroup.append(t_cluster)
                last_neighbours = set([])
                for bc in t_cluster:
                    last_neighbours.union(nodes_list[bc].get_neighbours())
                temp_neighbours.append(last_neighbours)
                bulk_clusters = [bulk_clusters[ib] for ib in xrange(len(bulk_clusters)) if ib != target]
                centroids = [centroids[ic] for ic in xrange(len(centroids)) if ic != target]
        while len(temp_subgroup) < n:
            temp_subgroup.append([])
        for cluster in temp_subgroup:
            parallel_clusters.append(cluster)
    return parallel_clusters


if __name__ == "__main__":

    #TODO to obtain better results on the TVZ, plane equation must be coupled with area of definition (some plane don't cut through the whole mesh)
    # pool = ProcessingPool(processes= nb_process)
    # no_multi = []
    # multi = []
    # n_nodes = []
    # hist_data = []
    evol_dev = []



    for index, filename in enumerate(filenames):
    #     chunk_size = chunk_sizes[index]
    #     fig = plt.figure()
    #     plist = [0.1, 0.2, 0.3, 0.4, 0.5]
    #     order2 = combinations(plist, 2)
    #     order3 = combinations(plist, 3)
    #     tasks = [order2, order3]

            # ##### Test complex topography #####
            # #Creating the topo
            # # creates a random topo of 100*100 points
            # nx = 100
            # ny = 100
            # dem1 = numpy.random.rand(nx, ny)
            #
            # sizex = 40
            # sizey = 40
            # x, y = numpy.mgrid[-sizex:sizex + 1, -sizey:sizey + 1]
            # g = numpy.exp(-0.333 * (x ** 2 / float(sizex) + y ** 2 / float(sizey)))
            # filter = g / g.sum()
            #
            # demSmooth = scipy.signal.convolve(dem1, filter, mode='valid')
            #
            # # rescales so elevation lies between 499 and 540
            # data = 499 + 41 * (demSmooth - demSmooth.min()) / (demSmooth.max() - demSmooth.min())
            # data_z = data.flatten()
            #
            # lenght = len(data_z)
            # data_x = numpy.zeros(lenght)
            # data_y = numpy.zeros(lenght)
            # new_x = -1000. / (len(data) - 1)
            # new_y = 0.
            # counter = 0
            # # adjusts the x,y range so that the surface is within 0 and 1
            # for i in xrange(len(data)):
            #     new_x += 1000. / (len(data) - 1)
            #     new_y = 0.
            #     for j in xrange(len(data[0])):
            #         data_x[counter] = new_x
            #         data_y[counter] = new_y
            #         new_y += 1000. / (len(data[0]) - 1)
            #         counter += 1

            #Formating boundary data
            # topo = scipy.interpolate.bisplrep(data_x, data_y, data_z, s=0)

            #TODO compute the bspline of the topo from the mesh so that it can be used as a boundary


            # X,Y = numpy.meshgrid(numpy.linspace(100,160965,100),numpy.linspace(100,86605,100))
            # XX = X.flatten()
            # YY = Y.flatten()
            # ZZ = numpy.array([scipy.interpolate.bisplev(XX[i],YY[i],dict_bounds['complex'][0][0]) for i in xrange(len(XX))])
            # Z = ZZ.reshape(X.shape)
            # fig = plt.figure()
            # ax = fig.gca(projection='3d')
            # ax.plot_wireframe(X, Y, Z)
            # plt.show()
            #TODO might need to compute complex intersection before optimizing clusters
            # If cluster append to have only 1 or 2 point on a complex intersection we can't compute the spline!!!!

            # init_scatter = [[list_points[p][0] for p in points_to_plot], [list_points[p][1] for p in points_to_plot],
            #                  [list_points[p][2] for p in points_to_plot]]

            # ##### End format data #####
        # Data Import
        # stdout.write('Reading File...')
    # list_points = numpy.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    #                    [0.8, 0, 0.5], [0.2, 1, 0.5], [0, 0.5, 0.2], [1, 0.5, 0.8], [0.8, 0.5, 0], [0.2, 0.5, 1],
    #                    [0.2, 0.8, 0.2]])
    #
    # list_tets = [[0, 4, 8, 14], [4, 5, 8, 14], [5, 1, 8, 14], [1, 0, 8, 14], [1, 5, 11, 14], [5, 6, 11, 14],
    #                  [6, 2, 11, 14], [2, 1, 11, 14], [2, 6, 9, 14], [6, 7, 9, 14], [7, 3, 9, 14], [3, 2, 9, 14],
    #                  [3, 7, 10, 14], [7, 4, 10, 14], [4, 0, 10, 14], [0, 3, 10, 14], [0, 1, 12, 14], [1, 2, 12, 14],
    #                  [2, 3, 12, 14], [3, 0, 12, 14], [4, 5, 13, 14], [5, 6, 13, 14], [6, 7, 13, 14], [7, 4, 13, 14]]
    #
    # segments = [[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7],[0,8],[1,8],[4,8],[5,8],[1,11],[2,11],[5,11],[6,11],[2,9],[3,9],[6,9],[7,9],[3,10],[0,10],[7,10],[4,10],[0,12],[1,12],[2,12],[3,12],[4,13],[5,13],[6,13],[7,13],[0,14],[1,14],[2,14],[3,14],[4,14],[5,14],[6,14],[7,14],[8,14],[9,14],[10,14],[11,14],[12,14],[13,14]]
    #
    # fig1 = plt.figure()
    # wr = fig1.gca(projection='3d')
    # for seg in segments:
    #     x = [list_points[i][0] for i in seg]
    #     y = [list_points[i][1] for i in seg]
    #     z = [list_points[i][2] for i in seg]
    #     wr.plot(x,y,z, c='b', linewidth= 1, alpha=0.8)
    #
    # cells = {'tetra': numpy.array(list_tets)}
    # mesh = meshio.Mesh(list_points,cells)
    # meshio.write(MESH_PATH+filename, mesh)
    #     mesh = meshio.read(MESH_PATH + filename)
    #     # meshio.write(MESH_PATH+filename[:-4]+'.vtk',mesh)
    #
    #     list_points = mesh.points
        # chunk_sizes = numpy.logspace(0,3,num=15)
        # clean_size = []
        # n_nodes.append(len(list_points))
        # counter = 1
        # for task in tasks:
        #     for perm in task:
        labels = ('init',)
        mesh = meshio.read(MESH_PATH + filename)
        list_points = mesh.points
        chunk_size = int(0.01*len(list_points))
        cells = mesh.cells
        list_tets = cells['tetra']
        # print len(list_tets)

        stdout.write('\rTreating Boundary Input...')
        bounds_dict = format_bounds(complex_bounds, BOUND_DAT_PATH + plane_bounds, list_points)
        # print bounds_dict
        neighbours, el_neighbours = connection_mapping(list_tets)
        temp_vars = numpy.arange(0,len(list_points)+1,1)

        temp_cons = Con(neighbours, el_neighbours, list_tets, temp_vars)
        # print len(temp_cons.get_faces())
        temp_nodes = build_node_objects(list_points, temp_vars, bounds_dict)
        # nodes_to_plot = []
        # for n,node in enumerate(temp_nodes):
            # if 'complex' in node.get_bounds_type():
                # nodes_to_plot.append(n)
        # init_scatter = [[temp_nodes[n].get_x() for n in nodes_to_plot],[temp_nodes[n].get_y() for n in nodes_to_plot],[temp_nodes[n].get_z() for n in nodes_to_plot]]
        # point_to_plot = []
        # log = []
        # for node in temp_nodes:
        #     if 'complex' in node.get_bounds_type():
        #         point_to_plot.append(list_points[node.get_index()])
        #         log.append([node.get_bounds(), node.get_index()])
        # point_to_plot = numpy.array(point_to_plot)
        # X, Y = numpy.meshgrid(numpy.arange(0., 160965., 100.), numpy.arange(0., 86605., 100.))
        # XX = X.flatten()
        # YY = Y.flatten()
        # ZZ = numpy.array([scipy.interpolate.bisplev(XX[i],YY[i],bounds_dict['complex'][0][0]) for i in xrange(len(XX))])
        # Z = ZZ.reshape(X.shape)
        # print len(point_to_plot)
        #
        # fig = plt.figure()
        # ax = fig.gca(projection='3d')
        # ax.plot_wireframe(X, Y, Z)
        # ax.scatter(point_to_plot[:,0],point_to_plot[:,1],point_to_plot[:,2], c='r')
        # plt.show()
        # fig = plt.figure()
        # scat = fig.gca(projection='3d')
        if toggle_complex:
            complex_inter = compute_complex_intersect(temp_nodes, bounds_dict)
            # print complex_inter
        else:
            complex_inter = {}
        stdout.write('\rEvaluating initial quality of the mesh...')
        init_dev = quality_control(temp_nodes,temp_cons)
        evol_dev.append(init_dev)
    ##########WORKING PART TO UNCOMMENT
        # limits = [-1]
        # limits = [perm[len(perm)-i-1] for i in xrange(len(perm))] #, 0.15, 0.05
        limits = [0.4,0.3,0.1]
        # chunk_size = 1000
        # td = datetime.now() - init_time
        # stdout.write('\rRunning time : ' + '{:02d}'.format(td.seconds // 3600) + ':' + '{:02d}'.format(
        #         (td.seconds // 60) % 60) + ':' + '{:02d}'.format(td.seconds % 60))
        # init_time = clock()
        init_time = datetime.now()
        for i,lim in enumerate(limits):
            labels += (str(lim),)
            if multiprocess:

                temp_nodes = build_neighbours(temp_nodes, temp_cons)
                stdout.write('\nSelecting bad nodes...')
                # nodes_to_optimize = exctract_bad_elements(temp_nodes, temp_cons, lim)
                nodes_to_optimize = finding_targets(temp_nodes,temp_cons,lim,nb_process)
                # print nodes_to_optimize
                # chunk_matrix = random_chunking(len(list_points),30)
                # chunk_matrix = chunking_bad_elements(nodes_to_optimize, chunk_size)
                stdout.write('\nClustering...')
                # chunk_matrix = chunking_disjoint(temp_nodes, nodes_to_optimize, chunk_size, nb_process)

                chunk_matrix = chunking_disjoint2(list_points, temp_nodes, nodes_to_optimize, chunk_size, nb_process)
                # print [len(chunk) for chunk in chunk_matrix]
                stdout.write('\nCalibrating tolerance...')
                cal_tol_chunks = []
                while len(cal_tol_chunks) < 5:
                    index = numpy.random.randint(0,len(chunk_matrix))
                    if len(chunk_matrix[index]) != 0:
                        cal_tol_chunks.append(chunk_matrix[index])
                ftol = calibrate_tol(list_tets, list_points, cal_tol_chunks, bounds_dict, complex_inter)
                # ftol = 0.001
                stdout.write('\rWorking on level of refinement ' + str(i+1) + '/' + str(len(limits)) +': deviation limit is ' + str(lim) + ' and tolerance is ' + str(ftol) + '\nPlease Wait...')

                #Multiprocess
                if len(chunk_matrix) % nb_process == 0:
                    nb_iteration = len(chunk_matrix) // nb_process
                else:
                    nb_iteration = len(chunk_matrix) // nb_process + 1

                last_optimize = 0
                for j in xrange(nb_iteration):
                    stdout.write('\rOptimizing iterations : ' + str(j+1) + '/' + str(nb_iteration))
                    if last_optimize + nb_process <= len(chunk_matrix):
                        chunk_to_optimize = chunk_matrix[last_optimize:last_optimize+nb_process]
                    else:
                        chunk_to_optimize = chunk_matrix[last_optimize:]

                    chunk_to_optimize = [chunk for chunk in chunk_to_optimize if len(chunk) > 0]

                    input_data = pack_inputs(chunk_to_optimize, neighbours, el_neighbours, list_tets, bounds_dict, list_points, iter, ftol, complex_inter)
                    output = Queue()
                    processes = [Process(target=run_optimization, args=(data, output)) for data in input_data]
                    for p in processes:
                        p.start()
                    # for p in processes:
                    #     p.terminate()
                    # for p in processes:
                    #     p.join()
                    resulting_coords = [output.get() for ind in xrange(len(input_data))]

                    # resulting_coords = pool.map(run_optimization, input_data)
                    list_points = compile_results(list_points, resulting_coords)
                    last_optimize += nb_process

            else:
                # temp_nodes = build_neighbours(temp_nodes, temp_cons)
                # nodes_to_optimize = finding_targets(temp_nodes,temp_cons,lim,nb_process)
                # # chunk_matrix = random_chunking(len(list_points),30)
                # # chunk_matrix = chunking_bad_elements(nodes_to_optimize, chunk_size)
                # chunk_matrix = chunking_disjoint(temp_nodes, nodes_to_optimize, chunk_size)
                # chunk_matrix = [chunk for chunk in chunk_matrix if len(chunk)>0]
                # # stdout.write('\nCalibrating tolerance...')
                # cal_tol_chunks = []
                # while len(cal_tol_chunks) < 5:
                #     index = numpy.random.randint(0,len(chunk_matrix))
                #     if len(chunk_matrix[index]) != 0:
                #         cal_tol_chunks.append(chunk_matrix[index])
                # ftol = calibrate_tol(list_tets, list_points, [nodes_to_optimize], bounds_dict, complex_inter)
                ftol = 0.00001
                input_data = pack_inputs([temp_vars], neighbours, el_neighbours, list_tets, bounds_dict, list_points,
                                         iter, ftol, complex_inter)

                resulting_coords = run_optimization(input_data[0], None)
                list_points = compile_results(list_points, resulting_coords)
            td = datetime.now() - init_time
            stdout.write('\rRunning time : ' + '{:02d}'.format(td.seconds // 3600) + ':' + '{:02d}'.format(
                    (td.seconds // 60) % 60) + ':' + '{:02d}'.format(td.seconds % 60) + '\n')

            temp_nodes = build_node_objects(list_points, temp_vars, bounds_dict)
            # temp_nodes = update_complex_intersect(temp_nodes,complex_inter)
            final_dev = quality_control(temp_nodes, temp_cons)
            evol_dev.append(final_dev)

        # histo = fig.add_subplot(2,10,counter)
        # histo.hist(evol_dev,label= labels)
        # histo.legend(loc='best')
        # counter +=1
        # final_scatter = [[list_points[p][0] for p in points_to_plot], [list_points[p][1] for p in points_to_plot], [list_points[p][2] for p in points_to_plot]]
        new_filename = MESH_PATH + filename[:-4] + '_opti.msh'
        # r_time = datetime.now() - ref_time
        # print r_time, len(list_points)
        # run_times.append(r_time)
        # print datetime.now()-ref_time, len(list_points)
        # hist_data.append(evol_dev)
        # new_filename2 = MESH_PATH + filename[:-4] + '_opti.vtk'
        ####
        # mesh.points = list_points
        meshio.write(new_filename, mesh)
        # meshio.write(new_filename2, mesh)
        # for seg in segments:
        #     x = [list_points[i][0] for i in seg]
        #     y = [list_points[i][1] for i in seg]
        #     z = [list_points[i][2] for i in seg]
        #     wr.plot(x,y,z, c='r', linewidth= 1, alpha=0.8)
    #########################################################################

        ####PLOT COMPLEX SURFACE
        # final_scatter = [[list_points[n][0] for n in nodes_to_plot],[list_points[n][1] for n in nodes_to_plot],[list_points[n][2] for n in nodes_to_plot]]
    for r, result in enumerate(evol_dev):
        log = SAVE_PATH + fname + ext
        with open(log, 'w+') as out:
            for val in result:
                out.write(str(val) + '\n')
            out.close()
    fig1 = plt.figure()
    histo = fig1.add_subplot(1,1,1)
    # # labels = ('Initial State',)
    # # labels += ('Optimized State',)
    histo.hist(evol_dev, label= labels)
    histo.legend(loc= 'best')
    plt.title('Optimisation result for Model3') #' + str(len(list_points)) + '
    plt.xlabel('Deviation from orthogonality')
    plt.ylabel('Number of faces')
    plt.show()
        # fig2 = plt.figure()
        # scat = fig2.gca(projection='3d')
        # scat.scatter(init_scatter[0],init_scatter[1],init_scatter[2], c= 'b', s=50)
        # scat.scatter(final_scatter[0], final_scatter[1],final_scatter[2], c='r', s=50)
        # topo_data = numpy.loadtxt(TOPO_PATH + 'tdata.txt', delimiter=' ')
        # nx = len(set(list(topo_data[:,0])))
        # ny = len(set(list(topo_data[:,1])))
        # X,Y = numpy.meshgrid(numpy.linspace(0,10000,nx),numpy.linspace(0,10000,ny))
        # Z = topo_data[:,2].reshape(X.shape)
        # scat.plot_surface(X,Y,Z, alpha= 0.5, color= 'g')


    # fig1 = plt.figure()
    # plt.plot([650,1733,6871,21618],multi, c='b', label= 'Multiprocessing enabled') #
    # plt.plot([650,1733,6871], multi, c='b')
    # plt.xlabel('Mesh size (number of nodes)')
    # plt.ylabel('Running time of the optimisation (s)')
    # plt.title('Impact of multiprocessing on running times')
    # plt.legend(loc='best')

    # fig2 = plt.figure()
    # for ind, hist in enumerate(hist_data):
    #     histo = fig2.add_subplot(2,4,ind+1)
    #     histo.hist(hist[1:])
    # plt.hist(evol_dev, label=['Initial', 'Optimized']) #, 'Optimization of bad elements', 'Fine Optimization'
    # plt.xlabel('Deviation from orthogonality (1-cos(theta))')
    # plt.ylabel('Frequency')
    # plt.title('Distribution of deviation before and after optimization, mesh name: ' + filename + ', number of vertices: ' + str(len(list_points)) + ', running time: ' + str(td))
    # plt.legend()
    # fig = plt.figure()
    # plt.plot(nb_processes,run_times)
    # plt.xlabel('Number of processors')
    # plt.ylabel('Running time (s)')
    # plt.title('Running time versus Number of processors used')
    # plt.show()


    #
    # regular grid covering the domain of the data
    # X, Y = numpy.meshgrid(numpy.arange(0., 160965., 100.), numpy.arange(0., 86605., 100.))
    # XX = X.flatten()
    # YY = Y.flatten()
    # Z = numpy.array([scipy.interpolate.bisplev(XX[i], YY[i], dict_bounds['complex'][0][0]) for i in xrange(len(XX))]).reshape(X.shape)
    #
    # plot points and fitted surface
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.plot_wireframe(X, Y, Z)
    # ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2, cmap= cm.coolwarm)
    # ax.scatter(init_scatter[0], init_scatter[1], init_scatter[2], c='r', s=50)
    # ax.scatter(final_scatter[0], final_scatter[1], final_scatter[2], c='b', s=50)
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # ax.set_zlabel('Z')
    # ax.axis('equal')
    # ax.axis('tight')
    #
    # fig2 = plt.figure()
    # histo = fig2.add_subplot(1, 1, 1)
    # histo.hist(evol_dev)
    # plt.show()

