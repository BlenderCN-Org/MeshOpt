import numpy
from math import sin, cos, sqrt, atan2, radians, degrees, asin, tan
import scipy.interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import meshio


# def calculate_lenght(d_lat1, d_lon1, d_lat2, d_lon2):
#
#     lat1 = radians(d_lat1)
#     lon1 = radians(d_lon1)
#     lat2 = radians(d_lat2)
#     lon2 = radians(d_lon2)
#     R = 6371.0
#     dlon = lon2 - lon1
#     dlat = lat2 - lat1
#     a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
#     c = 2 * atan2(sqrt(a), sqrt(1 - a))
#     return R * c
#
#
# def apply_topo():
#     for point in list_points:
#         if abs(point[2]) < 0.001:
#             point[2] = fit(point[1], point[0])


if __name__ == "__main__":

    PATH = '/home/lmar626/Documents/Meshes/TestProblems/Scaling/'
    mname = 'model1.msh'
    tname = 'tdata.txt'

    mesh = meshio.read(PATH + mname)
    topo_data = numpy.loadtxt(PATH + tname)
    # real_x, real_y =X,Y = numpy.meshgrid(numpy.linspace(0,160965,42),numpy.linspace(0,86605,23))
    # XX= real_x.flatten()
    # YY=real_y.flatten()
    # with open(PATH+'TVZ_tdata.txt', 'w+') as out:
    #     for i in xrange(len(XX)):
    #         out.write(str(XX[i]) + ' ' + str(YY[i]) + ' ' + str(topo_data[i,2]) + '\n')
    #     out.close()

    fit = scipy.interpolate.bisplrep(topo_data[:, 0], topo_data[:, 1], topo_data[:, 2])

    points = mesh.points

    for point in points:
        if abs(point[2]) < 0.001:
            point[2] = scipy.interpolate.bisplev(point[0], point[1], fit)

    mesh.points = points
    meshio.write(PATH + '650_CS.msh', mesh)
    # X,Y = numpy.meshgrid(numpy.linspace(0,160965,200),numpy.linspace(0,86605,100))
    # XX = X.flatten()
    # YY = Y.flatten()
    # Z = numpy.array([scipy.interpolate.bisplev(XX[i],YY[i],fit) for i in xrange(len(XX))]).reshape(X.shape)
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.plot_surface(X,Y,Z)
    # plt.show()


    # MESH_PATH = '/home/mltn/Documents/Meshes/TVZ_models/640000/'
    # TOPO_PATH = '/home/mltn/Documents/Meshes/TVZ_models/Topo/'
    # toponame = 'data_topo.txt'
    # meshname = 'TVZ554000.msh'
    #
    # NW_vertice = [-37.629411, 176.108295] # Mount maunganui
    # NE_vertice = [-37.944351, 177.009664] # Estuaire Whakatane River
    # SE_vertice = [-39.208884, 176.108295] # West of Kaweka Forest Park
    # nb_points = 1000000
    #
    # NS_dist = calculate_lenght(NE_vertice[0], NE_vertice[1], SE_vertice[0], SE_vertice[1])
    # EW_dist = calculate_lenght(NE_vertice[0], NE_vertice[1], NW_vertice[0], NW_vertice[1])
    # dist_ratio = NS_dist / EW_dist
    # point_NS = int(numpy.sqrt(nb_points / dist_ratio))
    # point_EW = int(dist_ratio * point_NS)
    #
    # data = numpy.loadtxt(TOPO_PATH + toponame)
    #
    # # data_x = numpy.array([data[i][0] for i in xrange(len(data))])
    # # data_y = numpy.array([data[j][1] for j in xrange(len(data))])
    # pdx = numpy.linspace(0, EW_dist, point_EW)
    # pdy = numpy.linspace(0, NS_dist, point_NS)
    #
    # X, Y = numpy.meshgrid(pdx, pdy)
    # XX = X.flatten()
    # YY = Y.flatten()
    # data_z = numpy.array([data[k][2]/1000 for k in xrange(len(data))])
    #
    # with open(TOPO_PATH+'data_topo.xyz', 'w') as output:
    #
    #     for i in xrange(len(XX)):
    #         output.write(str(XX[i])+' '+str(YY[i])+' '+str(data_z[i])+'\n')

    # fit = scipy.interpolate.Rbf(XX, YY, data_z)
    # list_points, cells, point_data, cell_data, field_data = meshio.read(MESH_PATH + meshname)
    # print len(list_points)
    # apply_topo()
    # new_meshname = meshname[:-4] + '_wTopo.msh'
    # new_meshname1 = meshname[:-4] + '_wTopo.vtk'
    # meshio.write(MESH_PATH + new_meshname, list_points, cells, point_data= point_data, cell_data= cell_data, field_data= field_data)
    # meshio.write(MESH_PATH + new_meshname1, list_points, cells, point_data= point_data, cell_data= cell_data, field_data= field_data)
# Plot surface
#     U = numpy.linspace(0, 1000 * EW_dist, 371)
#     V = numpy.linspace(0, 1000 * NS_dist, 185)
#     U, V = numpy.meshgrid(U, V)
#     UU = U.flatten()
#     VV = V.flatten()
#     ZZ = numpy.array([fit(UU[i], VV[i]) for i in xrange(len(UU))])
#     Z = ZZ.reshape(U.shape)
#
#     fig = plt.figure()
#     ax = fig.gca(projection='3d')
#     ax.plot_surface(U, V, Z, rstride=1, cstride=1, cmap=cm.jet, linewidth=1, antialiased=True)
#     plt.show()