import pygmsh
import meshio
import numpy

GMSH_PATH = '/home/lmar626/Documents/gmsh/gmsh'
MESH_PATH = '/home/lmar626/Documents/Meshes/TestProblems/'
ext = '.msh'

geo = pygmsh.built_in.Geometry()

basement_points = [[0,0,-1500],[0,15000,-1500],[0,16732,-4500],[0,16932,-4500],[0,16066,-3000],[0,30000,-3000],[0,31732,-6000],[0,31932,-6000],[0,31066,-4500],[0,55531,-4500],[0,54672,-6000],[0,54872,-6000],[0,56605,-3000],[0,70538,-3000],[0,69672,-4500],[0,69872,-4500],[0,71605,-1500],[0,86605,-1500],[0,86605,-8000],[0,0,-8000]]
top_points = [[0,0,-1500],[0,15000,-1500],[0,15200,-1500],[0,16066,-3000],[0,30000,-3000],[0,30200,-3000],[0,31066,-4500],[0,55531,-4500],[0,56405,-3000],[0,56605,-3000],[0,70538,-3000],[0,71405,-1500],[0,71605,-1500],[0,86605,-1500],[0,86605,0],[0,0,0]]
f1_points = [[0,15000,-1500],[0,15200,-1500],[0,16932,-4500],[0,16732,-4500]]
f2_points = [[0,30000,-3000],[0,30200,-3000],[0,31932,-6000],[0,31732,-6000]]
f3_points = [[0,71405,-1500],[0,71605,-1500],[0,69872,-4500],[0,69672,-4500]]
f4_points = [[0,56405,-3000],[0,56605,-3000],[0,54872,-6000],[0,54672,-6000]]

top = geo.add_polygon(top_points, lcar=1200)
basement = geo.add_polygon(basement_points,lcar= 1200)
f1 = geo.add_polygon(f1_points,lcar= 200)
f2 = geo.add_polygon(f2_points,lcar= 200)
f4 = geo.add_polygon(f3_points,lcar= 200)
f3 = geo.add_polygon(f4_points,lcar= 200)

axis = [160965, 0, 0]
vol_basement = geo.extrude(basement,translation_axis=axis, num_layers=200)
vol_top = geo.extrude(top,translation_axis=axis, num_layers=200)
vol_f1 = geo.extrude(f1,translation_axis=axis, num_layers=200)
vol_f2 = geo.extrude(f2,translation_axis=axis, num_layers=200)
vol_f3 = geo.extrude(f3,translation_axis=axis, num_layers=200)
vol_f4 = geo.extrude(f4,translation_axis=axis, num_layers=200)

# geo.add_physical_volume([vol_basement], label= 1)
# geo.add_physical_volume([vol_top], label= 2)
# geo.add_physical_volume([vol_f1,vol_f2,vol_f3,vol_f4], label= 3)

points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geo, gmsh_path= GMSH_PATH, geo_filename='test.geo')
# print len(points)
# meshname = 'TVZ' + str(len(points)) + ext
# meshio.write_points_cells(MESH_PATH + meshname, points, cells, point_data, cell_data, field_data)