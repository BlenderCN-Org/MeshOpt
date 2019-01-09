import numpy as np


def get_cross_product(vector1, vector2):
    "Returns the cross product of two vector to access the normal vector"

    n_vect = [vector1[1] * vector2[2] - vector1[2] * vector2[1],
              -1 * (vector1[0] * vector2[2] - vector1[2] * vector2[0]),
              vector1[0] * vector2[1] - vector1[1] * vector2[0]]
    return n_vect


def get_plane_equation(point_list):
    "Returns the coefs (a,b,c,d) of the equation of a given plane"

    vect1 = [point_list[1][i] - point_list[0][i] for i in range(3)]
    vect2 = [point_list[1][i] - point_list[2][i] for i in range(3)]
    norm1 = np.linalg.norm(vect1)
    norm2 = np.linalg.norm(vect2)
    vect1 = [coord/norm1 for coord in vect1]
    vect2 = [coord/norm2 for coord in vect2]
    cross_product = get_cross_product(vect1, vect2)
    coefs = [cross_product[0], cross_product[1], cross_product[2], - (
            cross_product[0] * point_list[1][0] + cross_product[1] * point_list[1][1] + cross_product[2] *
            point_list[1][2])]
    # print coefs, sympy.Plane(sympy.Point3D(point_list[1]),sympy.Point3D(point_list[0]),sympy.Point3D(point_list[2])).equation()
    return coefs


TOPO_PATH = 'C:\Users\lmar626\Documents\Documents\Meshes\TVZ_models\Topo'
DATA_PATH = 'C:\Users\lmar626\Documents\Documents\Meshes\TVZ_models'
toponame = '\data_topo1000000.txt'
litho = ['faults', 'res']

#Litho
l1 = [[0., 0., -1500.], [0., 15200., -1500.], [160965., 15200., -1500.], [160965., 0., -1500.]]
l2 = [[0., 16066., -3000.], [0., 30200., -3000.], [160965., 30200., -3000.], [160965., 16066., -3000.]]
l3 = [[0., 31066., -4500.], [0., 55531., -4500.], [160965., 55531., -4500.], [160965., 31066., -4500.]]
l4 = [[0., 56405., -3000.], [0., 70538., -3000.], [160965., 70538., -3000.], [160965., 56405., -3000.]]
l5 = [[0., 71405., -1500.], [0., 86605., -1500.], [160965., 86605., -1500.], [160965., 71405., -1500.]]
#Faults
f1 = [[0., 16732., -4500.], [160965., 16732., -4500.], [160965., 15000., -1500.], [0., 15000., -1500.]]
f2 = [[0., 16932., -4500.], [160965., 16932., -4500.], [160965., 15200., -1500.], [0., 15200., -1500.]]
f3 = [[0., 16732., -4500.], [160965., 16732., -4500.], [160965., 16932., -4500.], [0., 16932., -4500.]]

f4 = [[0., 31732., -6000.], [160965., 31732., -6000.], [160965., 30000., -3000.], [0., 30000., -3000.]]
f5 = [[0., 31932., -6000.], [160965., 31932., -6000.], [160965., 30200., -3000.], [0., 30200., -3000.]]
f6 = [[0., 31732., -6000.], [160965., 31732., -6000.], [160965., 31932., -6000.], [0., 31932., -6000.]]

f7 = [[0., 54872., -6000.], [160965., 54872., -6000.], [160965., 56605., -3000.], [0., 56605., -3000.]]
f8 = [[0., 54672., -6000.], [160965., 54672., -6000.], [160965., 56405., -3000.], [0., 56405., -3000.]]
f9 = [[0., 54672., -6000.], [160965., 54672., -6000.], [160965., 54872., -6000.], [0., 54872., -6000.]]

f10 = [[0., 69872., -4500.], [160965., 69872., -4500.], [160965., 71605., -1500.], [0., 71605., -1500.]]
f11 = [[0., 69672., -4500.], [160965., 69672., -4500.], [160965., 71405., -1500.], [0., 71405., -1500.]]
f12 = [[0., 69672., -4500.], [160965., 69672., -4500.], [160965., 69872., -4500.], [0., 69872., -4500.]]

x_res = 1501
y_res = 751
z_res = 201
X = np.linspace(0,160965,x_res)
Y = np.linspace(0,86605,y_res)
Z = np.linspace(-8000,2000,z_res)

bounds = [l1, l2, l3, l4, l5, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12]
planes_eqs = []
for bound in bounds:
    planes_eqs.append(get_plane_equation(bound[:-1]))
# print planes_eqs[-2]
# cond = '(0.86604*y0-0.50021*z0-62580>0)&&'
# bulk_topo = np.loadtxt(TOPO_PATH + toponame)
# pdx = np.linspace(0,160965,1362)
# pdy = np.linspace(0,86605,733)
#
# X_topo, Y_topo = np.meshgrid(pdx, pdy)
# XX = X_topo.flatten()
# YY = Y_topo.flatten()
#
# topo = [[XX[i],YY[i],bulk_topo[i][2]] for i in range(len(XX))]

# Bottom reservoir
# intr = []
# for x in X:
#     for y in Y:
#         test = False
#         if 0. <= y < 15200.:
#             eq = planes_eqs[0]
#             test = True
#         elif 16066. <= y < 30200.:
#             eq = planes_eqs[1]
#             test = True
#         elif 31066. <= y < 55531.:
#             eq = planes_eqs[2]
#             test = True
#         elif 56405. <= y < 70538.:
#             eq = planes_eqs[3]
#             test = True
#         elif 71605. <= y < 86605.:
#             eq = planes_eqs[4]
#             test = True
#         if test:
#             z = (-x * eq[0] - y * eq[1] - eq[3]) / eq[2]
#             intr.append([x,y,z])

# Bottom Fault
# b_fault = []
# for x in X:
#     for y in Y:
#         test = False
#         if 15000. <= y < 16732.:
#             eq = planes_eqs[5]
#             test = True
#         elif 16732. <= y < 16932.:
#             eq = planes_eqs[7]
#             test = True
#         elif 30000. <= y < 31732.:
#             eq = planes_eqs[8]
#             test = True
#         elif 31732. <= y < 31932.:
#             eq = planes_eqs[10]
#             test = True
#         elif 54672. <= y < 54872.:
#             eq = planes_eqs[13]
#             test = True
#         elif 54872. <= y < 56605.:
#             eq = planes_eqs[11]
#             test = True
#         elif 69672. <= y < 69872.:
#             eq = planes_eqs[16]
#             test = True
#         elif 69872. <= y < 71605.:
#             eq = planes_eqs[14]
#             test = True
#         if test:
#             z = (-x * eq[0] - y * eq[1] - eq[3]) / eq[2]
#             b_fault.append([x,y,z])
#
# # Top Fault
# t_fault = []
# for x in X:
#     for y in Y:
#         test = False
#         if 16066. <= y < 16932.:
#             eq = planes_eqs[6]
#             test = True
#         elif 15200. <= y < 16066.:
#             eq = planes_eqs[6]
#             test = True
#         elif 30200. <= y < 31066.:
#             eq = planes_eqs[9]
#             test = True
#         elif 31066. <= y < 31932.:
#             eq = planes_eqs[9]
#             test = True
#         elif 55531. <= y < 56405.:
#             eq = planes_eqs[12]
#             test = True
#         elif 54672. <= y < 55531.:
#             eq = planes_eqs[12]
#             test = True
#         elif 70538. <= y < 71405:
#             eq = planes_eqs[15]
#             test = True
#         elif 69672. <= y < 70538.:
#             eq = planes_eqs[15]
#             test = True
#         if test:
#             z = (-x * eq[0] - y * eq[1] - eq[3]) / eq[2]
#             t_fault.append([x,y,z])


# Domain Limits
box = []
for x in X:
    for y in Y:
        box.append([x,y,-8000.])
        box.append([x,y,2000.])
for z in Z:
    for x in X:
        box.append([x,0.,z])
        box.append([x,86605.,z])
    for y in Y:
        box.append([0., y, z])
        box.append([160965., y, z])


with open(DATA_PATH+'\TVZ_cloud_box.xyz', 'w') as o_box:
    for point in box:
        o_box.write(str(point[0]) + ' ' + str(point[1]) + ' ' + str(point[2]) + '\n')
    o_box.close()
# with open(DATA_PATH+'\TVZ_cloud_topo.xyz', 'w') as o_topo:
#     for point in topo:
#         o_topo.write(str(point[0]) + ' ' + str(point[1]) + ' ' + str(point[2]) + '\n')
#     o_topo.close()
# with open(DATA_PATH + '\TVZ_cloud_intr.xyz', 'w') as o_intr:
#     for point in intr:
#         o_intr.write(str(point[0]) + ' ' + str(point[1]) + ' ' + str(point[2]) + '\n')
#     o_intr.close()
# with open(DATA_PATH + '\TVZ_cloud_bfault.xyz', 'w') as o_bfault:
#     for point in b_fault:
#         o_bfault.write(str(point[0]) + ' ' + str(point[1]) + ' ' + str(point[2]) + '\n')
#     o_bfault.close()
# with open(DATA_PATH + '\TVZ_cloud_tfault.xyz', 'w') as o_tfault:
#     for point in t_fault:
#         o_tfault.write(str(point[0]) + ' ' + str(point[1]) + ' ' + str(point[2]) + '\n')
#     o_tfault.close()





