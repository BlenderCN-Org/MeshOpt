import meshio
import numpy as np

#hd max disp = 0.2; md max disp = 0.1; ld max disp = 0.05
PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Simplest/24tet/'
model = '24tet.msh'
mesh = meshio.read(PATH+model)

points = mesh.points
for point in points:
    test = 0
    pos = [0,0,0]
    for i,c in enumerate(point):
        if c == 0.25 or c == 0.75:
            test += 1
            pos[i] = 1
    if test == 2:
        rc = [0.01* np.random.randint(1, 5) for it in xrange(2)]
        si = [np.random.random() for ic in xrange(2)]
        ind = 0
        for j,p in enumerate(pos):
            if p == 1:
                sign = 1.0
                if si[ind]<0.5:
                    sign = -1.0
                point[j] = point[j] + sign*rc[ind]
                ind += 1
    elif test == 3:
        rc = [0.01 * np.random.randint(1, 5) for it in xrange(3)]
        si = [np.random.random() for ic in xrange(3)]
        ind = 0
        for j,p in enumerate(pos):
            if p == 1:
                sign = 1.0
                if si[ind]<0.5:
                    sign = -1.0
                point[j] = point[j] + sign*rc[ind]
                ind += 1
mesh.points = points
meshio.write(PATH+'ld.msh', mesh)