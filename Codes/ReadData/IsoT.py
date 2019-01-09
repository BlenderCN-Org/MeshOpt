import numpy as np
import meshio
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate

PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Model1/'
end = ['hex/', 'Tet/cyl/']
names = ['hex_m1.vtk', 'tet_m1.vtk', 'otet_m1.vtk']
types = ['hex', 'tet']
SAVE = '/home/lmar626/Documents/Meshes/Basicmodels/Results/isoT/'
sn = ['hex_T170.xyz','tet_T170.xyz','otet_T170.xyz']
tasks = [[0,1,1],[0,1,2]]

def extract_T(Temp, full_PATH, type):
    mesh = meshio.read(full_PATH)
    p = mesh.points
    if type == 'hex':
        el = mesh.cells['hexahedron']
        T = mesh.cell_data['hexahedron']['T_1']
    else:
        el = mesh.cells['tetra']
        T = mesh.cell_data['tetra']['T_1']

    coords = []
    for i,t in enumerate(T):
        if abs(t-Temp)<5.:
            cur_p = np.array([p[j] for j in el[i]])
            if type == 'hex':
                centroid = [(min(cur_p[:,k])+max(cur_p[:,k]))/2 for k in xrange(3)]
            else:
                centroid = [sum(cur_p[:,k])/4 for k in xrange(3)]
            coords.append(centroid)
    return np.array(coords)

for n in xrange(3):
    cloud = extract_T(170., PATH+end[tasks[0][n]]+names[tasks[1][n]], types[tasks[0][n]])
    # fit = scipy.interpolate.bisplrep(cloud[:,0],cloud[:,1],cloud[:,2])
    #
    # X,Y = np.meshgrid(np.linspace(0.,1000.,20),np.linspace(0.,1000.,20))
    # XX, YY = X.flatten(), Y.flatten()
    # Z = np.array([scipy.interpolate.bisplev(XX[m],YY[m],fit) for m in xrange(len(XX))])
    with open(SAVE+sn[tasks[1][n]], 'w+') as out:
        for cl in cloud:
            out.write(str(cl[0]) + ' ' + str(cl[1]) + ' ' + str(cl[2]) + '\n')
        out.close()
    fig = plt.figure()
    pl = fig.gca(projection= '3d')
    pl.scatter(cloud[:,0],cloud[:,1],cloud[:,2])
plt.show()