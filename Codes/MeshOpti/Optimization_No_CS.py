import matplotlib.pyplot as plt
import meshio
from sys import stdout
from itertools import permutations
import numpy
from scipy.optimize import minimize
from datetime import datetime
from collections import OrderedDict
#from pathos.multiprocessing import ProcessingPool
from multiprocessing import Process
from fmq import Queue

# Input variables

# Mesh file name
filenames = ['model3.msh'] #, 'Model2.msh', 'model3.msh'

# Maximum number of iteration during optimization
iter = 1000



# Class Con

class Con:

    def __init__(self, dict_con, el_con, el_list,vars_list):
        self.faces_ = []
        self.vertices_ = []
        self.build_cons(dict_con, el_con, el_list, vars_list)

    def build_cons(self, dict_con, el_con, el_list, vars_list):
        vars_list = set(vars_list)
        size = len(dict_con)
        counter = 0
        self.faces_ = [[]] * size
        self.vertices_ = [[]] * size
        for con in dict_con:
            vertices = set(dict_con[con])
            if not vars_list.isdisjoint(vertices) or not vars_list.isdisjoint(con):
                self.vertices_[counter] = dict_con[con]
                face = [vert for vert in el_list[el_con[con][0]] if vert != dict_con[con][0]]
                face.sort()
                self.faces_[counter] = face
                counter += 1

    def get_faces(self): return self.faces_

    def get_vertices(self): return self.vertices_


# Class Nodes


class Node:

    def __init__(self, point_index, point_coord, node_bound_list, vars_list):
        self.global_index_ = point_index
        self.x_ = point_coord[0]
        self.y_ = point_coord[1]
        self.z_ = point_coord[2]
        self.alloc_vars_ = []
        self.neighbours_ = set([])
        self.bound_ = node_bound_list

        if point_index in vars_list:
            if len(self.bound_) < 3:
                self.determine_freedom()
            else:
                self.deg_freedom_ = 0
                self.freedom_ = None
        else:
            self.deg_freedom_ = 0
            self.freedom_ = None

        if self.freedom_ is not None:
            if self.deg_freedom_ != 3:
                self.determine_coefs()
            else:
                self.coefs_ = None

    def determine_freedom(self):
        self.deg_freedom_ = 3
        self.freedom_ = [1, 1, 1]
        nb_bound = len(self.bound_)
        if nb_bound == 1:
            i = 0
            while self.bound_[0][i] == 0:
                i += 1
            self.freedom_[i] = 0
            self.deg_freedom_ = 2
        elif nb_bound == 2:
            loc_axis = [0, 1, 2]
            i = 0
            while self.bound_[0][i] == 0:
                i += 1
            loc_axis = loc_axis[0:i] + loc_axis[i + 1:3]
            coefs = [self.bound_[0][j] * (-self.bound_[1][i] / self.bound_[0][i]) - self.bound_[1][j] for j
                     in range(4) if j != i]
            j = 0
            while coefs[j] == 0:
                j += 1
            self.freedom_[i] = 0
            self.freedom_[loc_axis[j]] = 0
            self.deg_freedom_ = 1

    def determine_coefs(self):

        index_1 = 0
        index_2 = 0
        for i in range(3):
            if self.freedom_[i] == 1:
                index_1 = i
            else:
                index_2 = i

        if self.deg_freedom_ == 1:
            k = 0
            while self.bound_[0][k] == 0:
                k += 1
            m = 3
            for j in range(3):
                if j != index_1 and j != k:
                    m = j
            l_coefs2 = [
                (self.bound_[1][k] * self.bound_[0][rang] / self.bound_[0][k] - self.bound_[1][rang]) / (
                        -self.bound_[1][k] * self.bound_[0][m] / self.bound_[0][k] + self.bound_[1][m]) for rang in
                [index_1, 3]]
            l_coefs1 = [
                -self.bound_[0][m] * l_coefs2[0] / self.bound_[0][k] - self.bound_[0][index_1] / self.bound_[0][k],
                -self.bound_[0][m] * l_coefs2[1] / self.bound_[0][k] - self.bound_[0][3] / self.bound_[0][k]]
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
            self.coefs_ = [-self.bound_[0][loc_vars[0]] / self.bound_[0][index_2],
                           -self.bound_[0][loc_vars[1]] / self.bound_[0][index_2],
                           -self.bound_[0][3] / self.bound_[0][index_2]]

    def allocation_var_indexes(self, current_lenght):

        if self.deg_freedom_ == 3:
            self.alloc_vars_ = [current_lenght, current_lenght + 1, current_lenght + 2]
        elif self.deg_freedom_ == 2:
            self.alloc_vars_ = [current_lenght, current_lenght + 1]
        else:
            self.alloc_vars_ = [current_lenght]

    def update_coords(self, current_vars):
        new_coords = [0, 0, 0]
        if self.deg_freedom_ == 3:
            new_coords[0] = current_vars[self.alloc_vars_[0]]
            new_coords[1] = current_vars[self.alloc_vars_[1]]
            new_coords[2] = current_vars[self.alloc_vars_[2]]
        elif self.deg_freedom_ == 2:
            loc_index = 0
            for i in range(3):
                if self.freedom_[i] == 1:
                    new_coords[i] = current_vars[self.alloc_vars_[loc_index]]
                    loc_index += 1
                else:
                    new_coords[i] = self.coefs_[0] * current_vars[self.alloc_vars_[0]] + self.coefs_[1] * current_vars[
                        self.alloc_vars_[1]] + self.coefs_[2]
        else:
            loc_index = 0
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

    def get_x(self):
        return self.x_

    def get_y(self):
        return self.y_

    def get_z(self):
        return self.z_

    def get_index(self):
        return self.global_index_

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


def build_node_objects(point_list, bound_list, vars_list):
    "Build a list of node objects from the list of points"

    nb_points = len(point_list)
    node_list = [0] * len(point_list)
    for i in xrange(nb_points):
        current_bounds = []
        for bound in bound_list:
            current_eq = abs(
                point_list[i][0] * bound[0] + point_list[i][1] * bound[1] + point_list[i][2] * bound[2] + bound[3])
            if current_eq < 0.0001:
                current_bounds.append(bound)
        node_list[i] = Node(i, point_list[i], current_bounds, vars_list)
    return node_list


def update_node_coords(node_list, var_coords):
    "Update the coords of the node objects that are being optimized"

    for node in node_list:
        if node.get_freedom() is not None:
            node.update_coords(var_coords)
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


def exctract_bad_elements(node_list, con_obj, limit):
    "Return a list of nodes index that play a role in deviations higher than the given limit"

    faces = con_obj.get_faces()
    vertices = con_obj.get_vertices()
    max_lenght = len(node_list)
    bad_nodes = -1 * numpy.ones(max_lenght)
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

            dev = 1 - abs(dot_product)
            if dev > limit:
                for j in xrange(3):
                    if face[j] not in bad_nodes:
                        bad_nodes[counter] = face[j]
                        counter += 1
                for k in xrange(2):
                    if vertice[k] not in bad_nodes:
                        bad_nodes[counter] = vertice[k]
                        counter += 1

    bad_nodes.resize(counter)
    return bad_nodes


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


def objective_function(N, nodes=None,connect=None):
    "Objective function: p-norm of the list of deviations"
    # function we want to minimize
    # N = list of variables
    # points = list of initial point coords
    # iteration = dictionary of connections
    # freedom = dictionary of degree of freedom of points
    # boundaries = dictionary of plane equations
    global update_status
    update_status = False
    if update_status is False:
        nodes = update_node_coords(nodes, N)
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


def build_jacobian(N, nodes=None,connect=None):
    "Build the analytical partial derivatives of the objective function wrt to each variables"
    # Matrix with one 1 lign and len(N) columns with the derivatives of obj function wrt each variables
    # N = list of variables, 3*number_of_points variables
    # sequence = dictionary obtained with connection_mapping method
    # vertices_indexes = dictionary obtained with connection_mapping_opti method
    # elements_indexes = list of vertices indexes for every tet in the mesh

    global update_status
    if update_status is False:
        nodes = update_node_coords(nodes, N)
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
                    current_bounds = nodes[point].get_coefs()
                    current_ders = [0., 0., 0.]
                    local_der_x = sign * matrix_ders[i][0]
                    local_der_y = sign * matrix_ders[i][1]
                    local_der_z = sign * matrix_ders[i][2]
                    if current_deg == 3:
                        current_ders[0] += local_der_x
                        current_ders[1] += local_der_y
                        current_ders[2] += local_der_z

                    elif current_deg == 2:
                        if current_freedom[2] == 0:
                            current_ders[0] += local_der_x + current_bounds[0] * local_der_z
                            current_ders[1] += local_der_y + current_bounds[1] * local_der_z
                        elif current_freedom[1] == 0:
                            current_ders[0] += local_der_x + current_bounds[0] * local_der_y
                            current_ders[2] += local_der_z + current_bounds[1] * local_der_y
                        else:
                            current_ders[1] += local_der_y + current_bounds[0] * local_der_x
                            current_ders[2] += local_der_z + current_bounds[1] * local_der_x

                    else:
                        if current_freedom[0] == 1:
                            current_ders[0] += local_der_x + current_bounds[0][0] * local_der_y + current_bounds[1][
                                0] * local_der_z
                        elif current_freedom[1] == 1:
                            current_ders[1] += local_der_y + current_bounds[0][0] * local_der_x + current_bounds[1][
                                0] * local_der_z
                        else:
                            current_ders[2] += local_der_z + current_bounds[0][0] * local_der_x + current_bounds[1][
                                0] * local_der_y
                    current_alloc = nodes[point].get_alloc_vars()
                    counter = 0
                    for i in xrange(3):
                        #Finding the maximum index for the derivatives introduced in the jacobian so that we can resize it properly
                        if current_freedom[i] == 1:
                            jac[current_alloc[counter]] += current_ders[i]
                            counter += 1

    return jac


def run_optimization(data, output):
    "Algorithm that creates the necessary objects and then launches the optimization, the output is added to the current queue of results"

    vars_indexes, dict_con, el_con, el_list, bound_list, points_list, maxiter, ftol = data
    cons = Con(dict_con, el_con, el_list,vars_indexes)
    nodes_list = build_node_objects(points_list, bound_list,vars_indexes)
    source_vars, nodes_list = initial_state(nodes_list)
    global update_status
    update_status = True

    result = minimize(objective_function, source_vars,
                                     args=(nodes_list,cons), method='SLSQP',
                                     jac=build_jacobian,
                                     options={'disp': False, 'ftol': ftol, 'maxiter': maxiter})

    nodes_list = update_node_coords(nodes_list, result.x)
    new_points_list = format_results(nodes_list)

    output.put(new_points_list)


def compile_results(point_list, results):
    "Rebuilds the list of points using the results from the previous chunks optimizations"

    for result in results:
        for point in result:
            point_list[point[3]][0] = point[0]
            point_list[point[3]][1] = point[1]
            point_list[point[3]][2] = point[2]
    return point_list


def pack_inputs(chunks, dict_con, el_con, el_list, bound_list, points_list, maxiter, ftol):
    "Compress the args needed in run_optimization into a list so that run_optimization can be called using only 1 arg"

    inputs = []
    for chunk in chunks:
        inputs.append([chunk, dict_con, el_con, el_list, bound_list, points_list, maxiter, ftol])
    return inputs


def calibrate_tol(el_list, bound_list, points_list, vars_indexes):

    ftol = 10
    nb_it = 1
    dict_con, el_con = connection_mapping(el_list)
    cons = Con(dict_con, el_con, el_list, vars_indexes)
    nodes_list = build_node_objects(points_list, bound_list, vars_indexes)
    source_vars, nodes_list = initial_state(nodes_list)
    global update_status
    update_status = True
    while nb_it == 1:
        result = minimize(objective_function, source_vars,
                          args=(nodes_list, cons), method='SLSQP',
                          jac=build_jacobian,
                          options={'disp': False, 'ftol': ftol, 'maxiter': 1})
        nb_it = result.nit
        ftol = ftol * 0.5
    return ftol


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


def chunking_disjoint(nodes_list, bad_nodes, sub_size):
    "Creates disjoint clusters of points so that those clusters can be optimized in parallel using multiprocessing"

    bad_nodes = list(bad_nodes)
    bulk_chunks = []
    counter = 0
    while len(bad_nodes) > 0:
        bulk_chunks += [[] for i in xrange(3)]
        infinite_loop = False
        while len(bulk_chunks[counter]) < sub_size and len(bulk_chunks[counter+1]) < sub_size and len(bulk_chunks[counter+2]) < sub_size and infinite_loop is False:
            ind = 0
            infinite_loop = True
            for node in bad_nodes:
                current_neighbours = nodes_list[int(node)].get_neighbours()
                test = True
                for i in xrange(3):
                    if counter + i != counter + ind:
                        for element in bulk_chunks[counter + i]:
                            temp_neighbours = nodes_list[int(element)].get_neighbours()
                            if not temp_neighbours.isdisjoint(current_neighbours) or node in temp_neighbours:
                                test = False

                if test and len(bulk_chunks[counter+ind]) < sub_size:
                    bulk_chunks[counter+ind].append(node)
                    bad_nodes.remove(node)
                    infinite_loop = False

                if ind == 2:
                    ind = 0
                else:
                    ind += 1
    #Rearrange chunk to maximize their lenght and minimize their number
        if len(bulk_chunks[counter]) + len(bulk_chunks[counter+1]) + len(bulk_chunks[counter+2]) <= sub_size:
            bulk_chunks[counter] += bulk_chunks[counter + 1] + bulk_chunks[counter+2]
            bulk_chunks[counter + 1] = []
            bulk_chunks[counter + 2] = []
        elif len(bulk_chunks[counter]) + len(bulk_chunks[counter+1]) <= sub_size:
            bulk_chunks[counter] += bulk_chunks[counter+1]
            bulk_chunks[counter+1] = []
        elif len(bulk_chunks[counter]) + len(bulk_chunks[counter+2]) <= sub_size:
            bulk_chunks[counter] += bulk_chunks[counter+1]
            bulk_chunks[counter+2] = []
        elif len(bulk_chunks[counter+1]) + len(bulk_chunks[counter+2]) <= sub_size:
            bulk_chunks[counter+1] += bulk_chunks[counter+2]
            bulk_chunks[counter+2] = []
        counter += 3

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


if __name__ == "__main__":

    # bounds = [[[1., 0., 0., 0.], [1., 0., 0., -10000.], [0., 1., 0., 0.], [0., 1., 0., -10000.], [0., 0., 1., 0.],
    #                    [0., 0., 1., 600.], [0., 0., 1., 1800.], [0., 0., 1., 4000.]],[[1., 0., 0., 0.], [1., 0., 0., -10000.], [0., 1., 0., 0.], [0., 1., 0., -10000.], [0., 0., 1., 0.],
    #                    [0., 0., 1., 600.], [0., 0., 1., 1800.], [0., 0., 1., 4000.], [1., 0., 0., -4900.], [1., 0., 0., -5100.], [0., 1., 0., -4900.], [0., 1., 0., -5100.]],[[1., 0., 0., 0.], [1., 0., 0., -10000.], [0., 1., 0., 0.], [0., 1., 0., -10000.], [0., 0., 1., 0.],
    #                    [0., 0., 1., 600.], [0., 0., 1., 1800.], [0., 0., 1., 4000.], [0.7071, 0., 0.7071, -1697.0563], [0.7071, 0., 0.7071, -1414.2136]]]
    bounds = [[[1., 0., 0., 0.], [1., 0., 0., -10000.], [0., 1., 0., 0.], [0., 1., 0., -10000.], [0., 0., 1., 0.],
                       [0., 0., 1., 600.], [0., 0., 1., 1800.], [0., 0., 1., 4000.], [0.7071, 0., 0.7071, -1697.0563], [0.7071, 0., 0.7071, -1414.2136]]]
    nb_process = 3
    # pool = ProcessingPool(processes= nb_process)
    fig = plt.figure()
    init_time = datetime.now()
    for bound,filename in enumerate(filenames):
        stdout.write('Reading File...')
        mesh = meshio.read(filename)
        list_points = mesh.points
        cells = mesh.cells
        # print len(list_points)
        list_tets = cells['tetra']
        list_bounds = bounds[bound]
        neighbours, el_neighbours = connection_mapping(list_tets)
        temp_vars = numpy.arange(0,len(list_points)+1,1)
        stdout.write('\rEvaluating initial quality of the mesh..')
        temp_cons = Con(neighbours, el_neighbours, list_tets, temp_vars)
        temp_nodes = build_node_objects(list_points, list_bounds, [])
        init_dev = quality_control(temp_nodes,temp_cons)
        evol_dev = [init_dev]

        # limits = [0.3, 0.2, 0.1, 0.05]
        limits = [0.3, 0.15, 0.05]
        # tols = [0.01, 0.01, 0.01, 0.001]
        chunk_size = 300

        td = datetime.now() - init_time
        stdout.write('\rRunning time : ' + '{:02d}'.format(td.seconds // 3600) + ':' + '{:02d}'.format(
            (td.seconds // 60) % 60) + ':' + '{:02d}'.format(td.seconds % 60))

        for i in xrange(3):

            temp_nodes = build_neighbours(temp_nodes, temp_cons)
            nodes_to_optimize = exctract_bad_elements(temp_nodes, temp_cons, limits[i])

            # chunk_matrix = random_chunking(len(list_points),30)
            # chunk_matrix = chunking_bad_elements(nodes_to_optimize, chunk_size)
            chunk_matrix = chunking_disjoint(temp_nodes, nodes_to_optimize, chunk_size)
            stdout.write('\nCalibrating tolerance...')
            ftol = calibrate_tol(list_tets, list_bounds, list_points, chunk_matrix[0])
            stdout.write('\rWorking on level of refinement ' + str(i+1) + '/3 : deviation limit is ' + str(limits[i]) + ' and tolerance is ' + str(ftol) + '\nPlease Wait...')
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
                input_data = pack_inputs(chunk_to_optimize, neighbours, el_neighbours, list_tets, list_bounds, list_points, iter, ftol)
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

            td = datetime.now() - init_time
            stdout.write('\rRunning time : ' + '{:02d}'.format(td.seconds // 3600) + ':' + '{:02d}'.format(
                (td.seconds // 60) % 60) + ':' + '{:02d}'.format(td.seconds % 60) + '\n')

            temp_nodes = build_node_objects(list_points, list_bounds, [])
            final_dev = quality_control(temp_nodes,temp_cons)
            evol_dev.append(final_dev)
        histo = fig.add_subplot(1,3,bound+1)
        histo.hist(evol_dev)
        new_filename = filename[:-4] + '_Opti.msh'
        meshio.write(new_filename, list_points, cells)
    plt.show()
