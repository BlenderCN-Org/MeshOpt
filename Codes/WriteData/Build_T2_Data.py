from mulgrids import *
from t2data import *
import t2data_json as t2j
import json
import meshio
import numpy as np

# #Inputs
PATH = '/home/lmar626/Documents/Meshes/Basicmodels/Model1/Waiwera/CRU/'
filename = 'm1_hex.dat'
inc_name = 'm1_hex_inc.incon'
geoname = 'gm1_hex.dat'
T2_PATH = '~/bin/'

write_dat = False
write_json = True

# mesh = meshio.read(PATH+'/Model1/test.vtk')
# meshio.write(PATH+'grid.msh',mesh)



# Geometry

fault = False
faults_eq = []
origin = [0.,0.,0.]

xy_lres = 400
xy_mres = 200
xy_hres = 100

z_lres = 200
z_mres = 100
z_hres = 50

xy_gridlayout = [xy_lres]*8 + [xy_mres]*6 + [xy_hres]*12 + [xy_mres]*6 + [xy_lres]*8
z_gridlayout = [z_hres]*4 + [z_mres]*6 + [z_lres]*16

#Data

upflow_radius = 1200.
center = np.array([5000.,5000.])

total_heat = 1.

#Build Geometry file
geo = mulgrid().rectangular(xy_gridlayout, xy_gridlayout, z_gridlayout, origin= origin, atmos_type= 0, convention= 0)

if write_dat:
    geo.write(filename=PATH +geoname)
    # geo.write_mesh(filename=PATH + 'grid.vtk')

    # # Build TOUGH2 Grid file
    # t2 = t2grid().fromgeo(geo= geo)
    # t2.add_rocktype(rocktype('res', permeability= [1e-15, 1e-15, 1e-15], porosity= 0.2))
    # t2.add_rocktype(rocktype('cap', permeability= [1e-15, 1e-15, 0.1e-15], porosity= 0.1))
    # t2.add_rocktype(rocktype('surf', permeability= [1e-15, 1e-15, 1e-15], porosity= 0.2))
    # if fault:
    #     t2.add_rocktype(rocktype('fault', permeability=[1e-15, 1e-15, 3e-15], porosity=0.2))
    #
    # for block in t2.blocklist[t2.num_atmosphere_blocks:]:
    #     if block.centre[2] < -1800: block.rocktype = t2.rocktype['res']
    #     elif -1800 <= block.centre[2] < -600: block.rocktype = t2.rocktype['cap']
    #     else:
    #         # if fault:
    #         #     for f in faults_eq:
    #         #         pos_block = f[0][0] * block.centre[0] + f[0][1] * block.centre[1] + f[0][2] * block.centre[1]
    #         #         if f[0][3] <= pos_block <= f[1][3]: block.rocktype = t2.rocktype['fault']
    #         block.rocktype = t2.rocktype['surf']



    #Build TOUGH2 Data file
    dat = t2data()
    dat.title = 'Model1_hex grid_naturalstate'
    dat.filename = PATH + filename
    # dat.filename = PATH + 't2_dat.dat'
    # dat.meshfilename = PATH + 't2_gdat.dat'
    dat.grid = t2grid().fromgeo(geo)

    dat.parameter.update(
        {'max_timesteps': 9999,
         'tstop': 1.e15,
         'const_timestep': 1.e7,
         'print_interval': 50,
         'gravity': 9.81,
         'default_incons': [1.013e5, 20.]})
    dat.start = True

    dat.parameter['option'][1] = 1
    dat.parameter['option'][16] = 5

    # dat.relative_permeability = {'type': 3, 'parameters': [0.3, 0.1, 0., 0., 0.]}
    # dat.capillarity = {'type': 1, 'parameters': [0., 0., 1., 0., 0.]}

    res = rocktype('res', permeability= [1e-15, 1e-15, 1e-15], porosity= 0.2)
    surf = rocktype('surf', permeability= [1e-15, 1e-15, 1e-15], porosity= 0.2)
    cap = rocktype('cap', permeability= [0.1e-15, 0.1e-15, 0.1e-15], porosity= 0.1)
    dat.grid.add_rocktype(res)
    dat.grid.add_rocktype(cap)
    dat.grid.add_rocktype(surf)
    if fault:
        dat.grid.add_rocktype(rocktype('fault', permeability=[1e-15, 1e-15, 3e-15], porosity=0.2))

    for block in dat.grid.blocklist[1:]:
        if block.centre[2] < -1800.:
            block.rocktype = res
        elif -1800. <= block.centre[2] < -600.:
            block.rocktype = cap
        else:
            # if fault:
            #     for f in faults_eq:
            #         pos_block = f[0][0] * block.centre[0] + f[0][1] * block.centre[1] + f[0][2] * block.centre[1]
            #         if f[0][3] <= pos_block <= f[1][3]: block.rocktype = t2.rocktype['fault']
            block.rocktype = surf

    dat.multi['eos'] = 'EW'


    dat.parameter['gravity'] = 9.81


    bottom_layer = geo.layerlist[-1]
    columns = [col for col in geo.columnlist if np.linalg.norm(col.centre - center) <= upflow_radius]
    total_area = sum([col.area for col in columns])

    block_heat = total_heat / total_area

    for i,col in enumerate(columns):
        block_name = geo.block_name(bottom_layer.name, col.name)
        gen_name = str(i)
        if len(gen_name) == 1:
            gen_name = 's  0' + gen_name
        elif len(gen_name) == 2:
            gen_name = 's  ' + gen_name
        elif len(gen_name) == 3:
            gen_name = 's ' + gen_name
        else:
            gen_name = 's' + gen_name
        gen = t2generator(name= gen_name, block = block_name, type ='HEAT', gx = block_heat * col.area)
        dat.add_generator(gen)


    # dat.end_keyword = 'ENDCY'
    dat.write()

    # inc = dat.grid.incons()
    # for blk in dat.grid.blocklist:
    #     # print blk.centre
    #     if blk.centre is None:
    #         inc[blk.name].variable = [1e5, 20.]
    #     else:
    #         inc[blk.name].variable = [dat.parameter['default_incons'][0], dat.parameter['default_incons'][1]]
    #
    #
    # inc.write(PATH + inc_name)
    # dat = t2data(PATH+'t2_dat.dat', PATH+'t2_gdat.dat')

    # dat.run(simulator=T2_PATH+'AUTOUGH2_42D')

elif write_json:
    convert = t2j.t2data_export_json(PATH+filename)
    convert.write_exodus_json(geo)
    # data = convert.json(geo, 'model1')
    # with open(PATH + 'model1hex_dat.json', 'w') as output:
    #     json.dump(data, output, sort_keys=True, indent=4, separators=(', ', ': '))