import bpy
import bmesh
from mathutils import Vector
from bpy_extras.object_utils import object_data_add
import numpy
from math import *

def calculate_lenght(d_lat1, d_lon1, d_lat2, d_lon2):

    lat1 = radians(d_lat1)
    lon1 = radians(d_lon1)
    lat2 = radians(d_lat2)
    lon2 = radians(d_lon2)
    R = 6371.0
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    return R * c


TOPO_PATH = '/home/mltn/Documents/Meshes/TVZ_models/Topo/'
toponame = 'data_topo10000.txt'

data = numpy.loadtxt(TOPO_PATH + toponame)
NW_vertice = [-37.629411, 176.108295]  # Mount maunganui
NE_vertice = [-37.944351, 177.009664]  # Estuaire Whakatane River
SE_vertice = [-39.208884, 176.108295]  # West of Kaweka Forest Park
nb_points = 10000

NS_dist = calculate_lenght(NE_vertice[0], NE_vertice[1], SE_vertice[0], SE_vertice[1])
EW_dist = calculate_lenght(NE_vertice[0], NE_vertice[1], NW_vertice[0], NW_vertice[1])
dist_ratio = NS_dist / EW_dist
point_NS = int(numpy.sqrt(nb_points / dist_ratio))
point_EW = int(dist_ratio * point_NS)
pdx = numpy.linspace(0, 1000 * EW_dist, point_EW)
pdy = numpy.linspace(0, 1000 * NS_dist, point_NS)

X, Y = numpy.meshgrid(pdx, pdy)
XX = X.flatten()
YY = Y.flatten()
ZZ = numpy.array([data[k][2] for k in range(len(data))])
verts = []
print ('Starting')
for i,el in enumerate(XX):
    verts.append([Vector((float(el),float(YY[i]),float(ZZ[i])))])
lmesh = bpy.data.meshes.new(name='Topo')
lmesh.verts.extend(verts)
lobj = bpy.object.new()
