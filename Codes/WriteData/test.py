import numpy
import meshio

names = ['tet.msh','otet.msh','hex.msh']
for name in names:
    FULL_PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Results/isoT/isoT170_'+name
    mesh = meshio.read('/home/lmar626/Documents/Meshes/Basicmodels/Results/isoT/isoT170_'+name)
    meshio.write(FULL_PATH[:-4]+'.vtk',mesh)

# def get_cross_product(vector1, vector2):
#     "Returns the cross product of two vector to access the normal vector"
#
#     n_vect = [vector1[1] * vector2[2] - vector1[2] * vector2[1],
#               -1 * (vector1[0] * vector2[2] - vector1[2] * vector2[0]),
#               vector1[0] * vector2[1] - vector1[1] * vector2[0]]
#     return n_vect
#
#
# def get_plane_equation(point_list):
#     "Returns the coefs (a,b,c,d) of the equation of a given plane"
#
#     vect1 = [point_list[1][i] - point_list[0][i] for i in xrange(3)]
#     vect2 = [point_list[1][i] - point_list[2][i] for i in xrange(3)]
#     norm1 = numpy.linalg.norm(vect1)
#     norm2 = numpy.linalg.norm(vect2)
#     vect1 = [coord/norm1 for coord in vect1]
#     vect2 = [coord/norm2 for coord in vect2]
#     cross_product = get_cross_product(vect1, vect2)
#     coefs = [cross_product[0], cross_product[1], cross_product[2], - (
#             cross_product[0] * point_list[1][0] + cross_product[1] * point_list[1][1] + cross_product[2] *
#             point_list[1][2])]
#     # print coefs, sympy.Plane(sympy.Point3D(point_list[1]),sympy.Point3D(point_list[0]),sympy.Point3D(point_list[2])).equation()
#     return coefs
#
# def checkpos(point):
#     d1 = eq1[0]*point[0]+eq1[1]*point[1]+eq1[2]*point[2]+eq1[3]
#     d2 = eq2[0]*point[0]+eq2[1]*point[1]+eq2[2]*point[2]+eq2[3]
#     if d1 < 0 and d2 >0:
#         return True
#     else:
#         return False
#
#
# f1 = [[2000.,0.,0.],[6000.,0.,-4000.],[6000.,10000.,-4000.]]
# f2 = [[2400.,0.,0.],[6400.,0.,-4000.],[6400.,10000.,-4000.]]
# eq1 = get_plane_equation(f1)
# eq2 = get_plane_equation(f2)
#
# t = [[1000.,0.,0.],[2200.,0.,0.],[3000.,0.,0.]]
# f = [checkpos(p) for p in t]
#
# print f