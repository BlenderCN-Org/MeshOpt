import numpy as np
import meshio
import quadpy

# def f(X):
#     return P - visc*q*X[0]/(dens*K)
def compute_volume(dat):
    if type == 'tet':
        vector1 = [dat[0][i]-dat[1][i] for i in xrange(3)]
        vector2 = [dat[0][j] - dat[2][j] for j in xrange(3)]
        vector3 = [dat[0][k] - dat[3][k] for k in xrange(3)]
        n_vect = [vector1[1] * vector2[2] - vector1[2] * vector2[1],
                  -1 * (vector1[0] * vector2[2] - vector1[2] * vector2[0]),
                  vector1[0] * vector2[1] - vector1[1] * vector2[0]]
        return abs(np.dot(n_vect,vector3))/6
    if type == 'hex':
        x_lenght = (max(dat[:,0])-min(dat[:,0]))
        y_lenght = (max(dat[:, 1]) - min(dat[:, 1]))
        z_lenght = (max(dat[:, 2]) - min(dat[:, 2]))
        return x_lenght*y_lenght*z_lenght

P = 5e5
T = 20.0
q = 0.5e-5
K = 1e-15
poro = 0.1
visc = 1.0017e-3
dens = 998.44

ROOT = '/home/lmar626/Documents/Meshes/Basicmodels/Simplest/'
PATHS = ['bigger/']
models = ['tet.vtk', 'tet_opti.vtk']
types = ['tet','hex']
tasks = [[0,0],[0,0]]

for m, model in enumerate(models):
    PATH = ROOT+PATHS[tasks[0][m]]
    type = types[tasks[1][m]]

    if type == 'tet':
        mesh = meshio.read(PATH+model)
        cell_data = mesh.cell_data
        els = mesh.cells['tetra']
        points = mesh.points
        P_th = []
        for el in els:
            dat = np.array([points[index] for index in el])
            P_th.append(quadpy.tetrahedron.integrate(lambda X: P - visc*q*X[0]/(dens*K),dat, quadpy.tetrahedron.Keast(10))/compute_volume(dat))
        errors = [(P_th[l]-cell_data['tetra']['P_1'][l])*0.00001 for l in xrange(len(P_th))]
        mesh.cell_data['tetra']['error'] = errors
        meshio.write(PATH+model, mesh)

    if type == 'hex':
        mesh = meshio.read(PATH+model)
        cell_data = mesh.cell_data
        els = mesh.cells['hexahedron']
        points = mesh.points
        P_th = []
        for el in els:
            dat = np.array([points[index] for index in el])
            x,y,z = [min(dat[:,0]),max(dat[:,0])],[min(dat[:,1]),max(dat[:,1])],[min(dat[:,2]),max(dat[:,2])]
            P_th.append(quadpy.hexahedron.integrate(lambda X: P - visc*q*X[0]/(dens*K),quadpy.hexahedron.cube_points(x,y,z),quadpy.hexahedron.Tyler(2))/compute_volume(dat))
        errors = [(P_th[l] - cell_data['hexahedron']['P_1'][l])*0.00001 for l in xrange(len(P_th))]
        mesh.cell_data['hexahedron']['error'] = errors
        meshio.write(PATH+model, mesh)