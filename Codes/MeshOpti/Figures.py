import numpy
import scipy.interpolate
import scipy.signal
import matplotlib.pyplot as plt
import plotly.tools as tls
from mpl_toolkits.mplot3d import Axes3D
import shapely
import meshio
from itertools import combinations
import Optimization_w_CS as ocs
import pygmsh

# mesh = meshio.read('/home/lmar626/Documents/Meshes/TestProblems/model1.msh')
# print len(mesh.points)

#introduction

# def centroid(p):
#     return [[(p[0][0]+p[1][0]+p[2][0])/3,(p[0][1]+p[1][1]+p[2][1])/3], [(p[0][0]+p[1][0]+p[3][0])/3,(p[0][1]+p[1][1]+p[3][1])/3]]
#
# p1 = [[0.,1.],[3.,1.],[3.,0.],[0.,2.]]
# p2 = [[1.,1.],[2.,1.],[1.5,0.],[1.5,2.]]
# segments = [[0,1],[0,2],[1,2],[0,3],[1,3]]
# c1 = centroid(p1)
# c2 = centroid(p2)
# fig = plt.figure()
# good = fig.add_subplot(1,2,2)
# for seg in segments:
#     x = [p2[i][0] for i in seg]
#     y = [p2[i][1] for i in seg]
#     good.plot(x,y,c='k')
# x = [c2[i][0] for i in xrange(2)]
# y = [c2[i][1] for i in xrange(2)]
# good.plot(x,y, c= 'r')
# good.scatter(x,y, c='r', s= 10)
# x = [p2[i][0] for i in xrange(len(p2))]
# y = [p2[i][1] for i in xrange(len(p2))]
# good.scatter(x,y, c= 'k', s=1.0)
#
# bad = fig.add_subplot(1,2,1)
# for seg in segments:
#     x = [p1[i][0] for i in seg]
#     y = [p1[i][1] for i in seg]
#     bad.plot(x,y,c='k')
# x = [c1[i][0] for i in xrange(2)]
# y = [c1[i][1] for i in xrange(2)]
# bad.plot(x,y,c='r')
# bad.scatter(x,y, c='r', s= 10)
# x = [p1[i][0] for i in xrange(len(p1))]
# y = [p1[i][1] for i in xrange(len(p1))]
# bad.scatter(x,y, c= 'k', s=1.0)
#
# plt.show()

#Spline

# x = numpy.linspace(0,5,5)
# y = numpy.array([numpy.random.randint(0,6) for i in xrange(len(x))])
# fit = scipy.interpolate.splrep(x,y)
# X = numpy.linspace(0,5,150)
# Y = [scipy.interpolate.splev(X[i],fit) for i in xrange(len(X))]
#
# fig = plt.figure()
# plt.scatter(x,y,c='r', s=10.)
# plt.plot(X,Y, c='b')
# plt.title('Graph showing how to fit a curve to a discrete cloud point using splines')
# plt.show()

#Clusters

# X,Y = numpy.meshgrid(numpy.arange(0,3,1),numpy.arange(0,3,1))
# XX = X.flatten()
# YY = Y.flatten()
# c = ['r', 'b', 'g']
# size = 20
# fig = plt.figure()
# random = fig.add_subplot(1,3,1)
# random.scatter([XX[0],XX[3],XX[7]],[YY[0],YY[3],YY[7]], c= c[0], s=size, label= 'Cluster #1')
# random.scatter([XX[8],XX[6],XX[5]],[YY[8],YY[6],YY[5]], c= c[1], s=size, label= 'Cluster #2')
# random.scatter([XX[1],XX[2],XX[4]],[YY[1],YY[2],YY[4]], c= c[2], s=size, label= 'Cluster #3')
# random.legend(loc = (0.,0.25))
# bad = fig.add_subplot(1,3,2)
# bad.scatter([XX[0],XX[3]],[YY[0],YY[3]], c= c[0], s=size, label= 'Cluster #1')
# bad.scatter([XX[8],XX[6],XX[5]],[YY[8],YY[6],YY[5]], c= c[1], s=size, label= 'Cluster #2')
# bad.scatter([XX[1],XX[2]],[YY[1],YY[2]], c= c[2], s=size, label= 'Cluster #3')
# bad.scatter([XX[7],XX[4]],[YY[7],YY[4]], c= 'k', s=size, label= 'Not to be optimized')
# bad.legend(loc = (0.,0.25))
# disjoint = fig.add_subplot(1,3,3)
# disjoint.scatter([XX[0],XX[3],XX[1]],[YY[0],YY[3],YY[1]], c= c[0], s=size, label= 'Cluster #1')
# disjoint.scatter([XX[2],XX[6],XX[4]],[YY[2],YY[6],YY[4]], c= c[1], s=size, label= 'Cluster #2')
# disjoint.scatter([XX[8],XX[7],XX[5]],[YY[8],YY[7],YY[5]], c= c[2], s=size, label= 'Cluster #3')
# disjoint.legend(loc = (0.,0.25))
# plt.show()

#Constraints

# p = [[0.,0.],[1.,0.],[1.,1.],[0.,1.],[0.5,0.5], [0.5,1.]]
# seg = [[0,1],[1,2],[0,3],[3,5],[5,2],[0,4],[1,4],[2,4],[3,4], [5,4]]
# fig = plt.figure()
# plt.scatter([p[0][0], p[4][0], p[5][0]], [p[0][1], p[4][1], p[5][1]], c= 'r', s= 20.)
# for s in seg:
#     x = [p[i][0] for i in s]
#     y = [p[i][1] for i in s]
#     plt.plot(x,y, c= 'k')
#
# plt.show()


#Complex intersection
# eq = [1.,0.,-1.,0.]
# PATH = '/home/lmar626/Documents/CartoOnTheGo/'
# topo_fname = 'topo.txt'
# topo_data = numpy.loadtxt(PATH + topo_fname, delimiter= ' ')
# inter = []
# low = []
# for point in topo_data:
#     c_eq = eq[0]*point[0]+eq[2]*point[2]+ eq[1]*point[1]
#     if abs(c_eq) < 0.005:
#         inter.append(point)
#     elif c_eq > 0.01:
#         low.append(point)
# inter = numpy.array(inter)
# low = numpy.array(low)
# nx = len(set(list(topo_data[:,0])))
# ny = len(set(list(topo_data[:,1])))
# tX,tY = numpy.meshgrid(numpy.linspace(0,1,nx),numpy.linspace(0,1,ny))
# tZ = topo_data[:,2].reshape(tX.shape)
#
#
# X,Y = numpy.meshgrid(numpy.linspace(0,1,150),numpy.linspace(0,1,150))
# XX = X.flatten()
# YY = Y.flatten()
# ZZ = numpy.array([XX[i] for i in xrange(len(XX))])
# Z = ZZ.reshape(X.shape)
#
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot_surface(X,Y,Z, alpha= 0.5, color='b')
# # ax.plot_surface(tX,tY,tZ, alpha= 0.5, color= 'g')
# ax.scatter(low[:,0],low[:,1],low[:,2], c='g', s=1.)
# ax.scatter(inter[:,0],inter[:,1],inter[:,2], c='r', s= 20)
# plt.show()

# list = [0.1,0.2,0.3,0.4,0.5]
# i= combinations(list, 3)
# t= combinations(list, 2)
# test1 = 0
# test2 = 0
# for i in xrange(2):
#     fig = plt.figure()
#     bla = fig.add_subplot(2,1,1)
# plt.show()







# PATH = '/home/lmar626/Documents/Meshes/TestProblems/Scaling/'
# topo_fname = 'tdata.txt'
#
#
#
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
# # rescales so elevation lies between 0 and 1
# Z = -500+1000*((demSmooth - demSmooth.min()) / (demSmooth.max() - demSmooth.min()))
# nx = Z.shape[0]
# ny = Z.shape[1]
# X,Y = numpy.meshgrid(numpy.linspace(0,10000,nx),numpy.linspace(0,10000,ny))
#
# XX = X.flatten()
# YY = Y.flatten()
# ZZ = Z.flatten()
# print len(XX)
# with open(PATH+'tdata.txt', 'w+') as out:
#     for i in xrange(len(XX)):
#         out.write(str(XX[i])+ ' ' + str(YY[i]) + ' ' + str(ZZ[i]) + '\n')
#     out.close()
# adjusts the x,y range so that the surface is within 0 and 1


# ZZ = r.flatten()

# X,Y = np.meshgrid(np.linspace(0,99,100),np.linspace(0,99,100))
# XX = X.flatten()
# YY = Y.flatten()



# PATH = '/home/lmar626/Documents/Meshes/TestProblems/'
# mname = 'TVZ10000.msh'
# tname = 'TVZ_tdata.txt'
#
# mesh = meshio.read(PATH + mname)
# topo_data = numpy.loadtxt(PATH + tname)
# # real_x, real_y =X,Y = numpy.meshgrid(numpy.linspace(0,160965,42),numpy.linspace(0,86605,23))
# # XX= real_x.flatten()
# # YY=real_y.flatten()
# # with open(PATH+'TVZ_tdata.txt', 'w+') as out:
# #     for i in xrange(len(XX)):
# #         out.write(str(XX[i]) + ' ' + str(YY[i]) + ' ' + str(topo_data[i,2]) + '\n')
# #     out.close()
#
# fit = scipy.interpolate.bisplrep(topo_data[:,0],topo_data[:,1],topo_data[:,2])
#
# points = mesh.points
#
# for point in points:
#     if abs(point[2]) < 0.001:
#         point[2] = scipy.interpolate.bisplev(point[0],point[1],fit)
#
# mesh.points = points
# meshio.write(PATH+'TVZ_CS10000.msh', mesh)
# # X,Y = numpy.meshgrid(numpy.linspace(0,160965,200),numpy.linspace(0,86605,100))
# # XX = X.flatten()
# # YY = Y.flatten()
# # Z = numpy.array([scipy.interpolate.bisplev(XX[i],YY[i],fit) for i in xrange(len(XX))]).reshape(X.shape)
# # fig = plt.figure()
# # ax = fig.gca(projection='3d')
# # ax.plot_surface(X,Y,Z)
# # plt.show()
#










# MESH_PATH = '/home/lmar626/Documents/Meshes/TestProblems/Scaling/'
# mesh_name = '650.msh'
# topo_data = numpy.loadtxt(PATH + topo_fname, delimiter= ' ')
# #
# #
#
# # #     for i in xrange(len(ZZ)):
# # #         out.write(str(topo_data[i][0]) + ' ' + str(topo_data[i][1]) + ' ' + str(ZZ[i]) + '\n')
# #
# fit = scipy.interpolate.bisplrep(topo_data[:,0],topo_data[:,1],topo_data[:,2])
# #
# mesh = meshio.read(MESH_PATH + mesh_name)
# list_point = mesh.points
# for point in list_point:
#     if abs(point[2])<0.001:
#         point[2] = scipy.interpolate.bisplev(point[0], point[1],fit)
# #
# mesh.point = list_point
# meshio.write(MESH_PATH+'650_CS.msh', mesh)
#
# # TPATH = '/home/lmar626/Documents/Meshes/TestProblems/Scaling/'
# # with open(TPATH + 'tdata.txt', 'w+') as out:
# #     for j in xrange(len(XX)):
# #         out.write(str(XX[j]) + ' ' + str(YY[j]) + ' ' + str(ZZ[j]) + '\n')
# #     out.close()
# #
# # fit = scipy.interpolate.bisplrep(X,Y,Z)
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# X, Y = numpy.meshgrid(numpy.linspace(0,10000,200),numpy.linspace(0,10000,200))
# XX, YY = X.flatten(), Y.flatten()
# Z = numpy.array([scipy.interpolate.bisplev(XX[i],YY[i],fit) for i in xrange(len(XX))]).reshape(X.shape)
# ax.plot_surface(X,Y,Z)
# plt.show()

# PATH = '/home/lmar626/Documents/Meshes/Articles/'
# filenames = ['MPr.txt']
# # labels = ['6871 nodes', '1733 nodes', '650 nodes']
# fig = plt.figure()
# ax = fig.add_subplot(1,1,1)
# for i, filename in enumerate(filenames):
#     data = numpy.loadtxt(PATH+filename, delimiter= ' ')
#     print data
#     speed_up = [data[0,0]/data[i,0] for i in xrange(len(data))]
#     # ratio = [data[j,1]/max(data[:,1]) for j in xrange(len(data))]
#     ax.scatter(data[:,1],speed_up,s=50, marker='P')
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# # ax.set_xbound(lower=0.0005,upper=2)
# ax.set_ybound(lower=0.8)
# plt.xlabel('Number of Processors', weight= 'bold', size= 'xx-large')
# plt.ylabel('Speed up', weight= 'bold', size= 'xx-large')
#
# for t in ax.get_xticklabels()+ax.get_yticklabels(): t.set_fontsize(15)
# plt.title('Speed up versus Number of processors used, using a cluster size of 1% of mesh size on a 6871 nodes mesh', weight='bold', size= 'xx-large')
# plt.legend(loc='best', fontsize = 'xx-large')
# plt.show()


#intr connection test
TPATH = '/home/lmar626/Documents/Meshes/Basicmodels/Model2/Tet/cyl/'
meshname = 'cyl_m2.msh'
mesh = meshio.read(TPATH+meshname)
cells = mesh.cells
# meshio.write('tet.exo', mesh)
cdata = mesh.cell_data['tetra']['gmsh:geometrical']
tets = cells['tetra']

n, en, = ocs.connection_mapping(tets)
# doublets = []
# for i,p1 in enumerate(mesh.points):
#     for j,p2 in enumerate(mesh.points):
#         if numpy.all(c1 == c2 for c1 in p1 for c2 in p2):
#             doublets.append((i,j))
# print len(i,j), len(mesh.points)
# print mesh.field_data
intr = {}
for con in en:
    if len(en[con]) == 2:
        if cdata[en[con][0]] != cdata[en[con][1]]:
            key = (cdata[en[con][0]],cdata[en[con][1]])
            if key in intr:
                intr[key] += 1
            else:
                intr[key] = 1
print intr
# meshio.write(TPATH+'test.vtk', mesh)

# obj = pygmsh.opencascade.Geometry()
