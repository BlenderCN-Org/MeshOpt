import meshio
import numpy
import json

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
                if model == 'm3':
                    f1 = [[2000.,0.,0.],[6000.,0.,-4000.],[6000.,10000.,-4000.]]
                    f2 = [[2400.,0.,0.],[6400.,0.,-4000.],[6400.,10000.,-4000.]]
                    eq1 = get_plane_equation(f1)
                    eq2 = get_plane_equation(f2)
                    d1 = eq1[0] * centroid[0] + eq1[1] * centroid[1] + eq1[2] * centroid[2] + eq1[3]
                    d2 = eq2[0] * centroid[0] + eq2[1] * centroid[1] + eq2[2] * centroid[2] + eq2[3]
                    if d1 <= 0 and d2 >= 0:
                        rt[i] = 4
                        found = True
                else:
                    if 4900. <= centroid[0] <= 5100. or 4900. <= centroid[1] <= 5100.:
                        rt[i] = 4
                        found = True
            if centroid[2] < -1800. and not found:
                rt[i] = 1
            elif -1800. <= centroid[2] < -600. and not found:
                rt[i] = 2
            elif -600. <= centroid[2] and not found:
                rt[i] = 3

    else:
        for i, el in enumerate(elements):
            dat = numpy.array([list_points[ind] for ind in el])
            centroid = [(max(dat[:,0])+min(dat[:,0]))/2, (max(dat[:,1])+min(dat[:,1]))/2, (max(dat[:,2])+min(dat[:,2]))/2]
            found = False
            if faults:
                if model == 'm3':
                    f1 = [[2000.,0.,0.],[6000.,0.,-4000.],[6000.,10000.,-4000.]]
                    f2 = [[2400.,0.,0.],[6400.,0.,-4000.],[6400.,10000.,-4000.]]
                    eq1 = get_plane_equation(f1)
                    eq2 = get_plane_equation(f2)
                    d1 = eq1[0] * centroid[0] + eq1[1] * centroid[1] + eq1[2] * centroid[2] + eq1[3]
                    d2 = eq2[0] * centroid[0] + eq2[1] * centroid[1] + eq2[2] * centroid[2] + eq2[3]
                    if d1 <= 0 and d2 >= 0:
                        rt[i] = 4
                        found = True
                else:
                    if 4900. <= centroid[0] <= 5100. or 4900. <= centroid[1] <= 5100.:
                        rt[i] = 4
                        found = True
            if centroid[2] < -1800. and not found:
                rt[i] = 1
            elif -1800. <= centroid[2] < -600.and not found:
                rt[i] = 2
            elif -600. <= centroid[2] and not found:
                rt[i] = 3

    return rt

PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Model1/Tet/cyl/'
meshname = 'otet_m1.vtk'
type = 'tet'
model = 'm3'
faults = False

mesh = meshio.read(PATH+meshname)
cells = mesh.cells
list_points = mesh.points
if type == 'tet':
    elements = cells['tetra']
else:
    elements = cells['hexahedron']

rt = numpy.zeros(len(elements))
ss = numpy.zeros(len(elements))
rt = fill_rocktypes()
with open(PATH+ meshname[:-4]+'.json', 'r') as f:
    data = json.load(f)
for source in data['source']:
    ss[source['cell']] = source['rate']

if type == 'tet':
    mesh.cell_data['tetra']['Sources'] = ss
    mesh.cell_data['tetra']['RT'] = rt
else:
    mesh.cell_data['hexahedron']['Sources'] = ss
    mesh.cell_data['hexahedron']['RT'] = rt
meshio.write(PATH+meshname, mesh)
