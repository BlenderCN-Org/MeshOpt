

PATH = '/home/lmar626/Documents/Meshes/TestProblems/'


# bounds = []
# Z = [0.,-600.,-1800.,-4000.]
# for z in Z:
#     bounds.append([[0.,0.,z],[2000.-z,0.,z],[2000.-z,00.,z],[0.,00.,z]])
#     bounds.append([[2400.-z,0.,z],[00.,0.,z],[00.,00.,z],[2400. -z,00.,z]])
#
# bounds.append([[0.,0.,-4000.],[0.,00.,-4000.],[0.,00.,0.], [0.,0.,0.]])
# bounds.append([[00.,0.,-4000.],[00.,00.,-4000.],[00.,00.,0.], [00.,0.,0.]])
# bounds.append([[0.,0.,-4000.],[00.,0.,-4000.],[00.,0.,0.],[0.,0.,0.]])
# bounds.append([[0.,00.,-4000.],[00.,00.,-4000.],[00.,00.,0.],[0.,00.,0.]])
#
# bounds.append([[2000.,0.,0.],[6000.,0.,-4000.],[6000.,00.,-4000.],[2000.,00.,0.]])
# bounds.append([[2400.,0.,0.],[6400.,0.,-4000.],[6400.,00.,-4000.],[2400.,00.,0.]])



# bot = [[0., 0., -8000.], [160965., 0., -8000.], [160965., 86605., -8000.], [0., 86605., -8000.]]
# x_min = [[0., 0., -8000.], [0., 86605., -8000.], [0., 86605., 2000.], [0., 0., 2000.]]
# x_max = [[160965., 0., -8000.], [160965., 86605., -8000.], [160965., 86605., 2000.], [160965., 0., 2000.]]
# y_min = [[0., 0., -8000.], [160965., 0., -8000.], [160965., 0., 2000.], [0., 0., 2000.]]
# y_max = [[0., 86605., -8000.], [160965., 86605., -8000.], [160965., 86605., 2000.], [0., 86605., 2000.]]
# # Litho
# l1 = [[0., 0., -1500.], [160965., 0., -1500.], [160965., 15200., -1500.], [0., 15200., -1500.]]
# l2 = [[0., 16066., -3000.], [160965., 16066., -3000.], [160965., 30200., -3000.], [0., 30200., -3000.]]
# l3 = [[0., 31066., -4500.], [160965., 31066., -4500.], [160965., 55531., -4500.], [0., 55531., -4500.]]
# l4 = [[0., 56405., -3000.], [160965., 56405., -3000.], [160965., 70538., -3000.], [0., 70538., -3000.]]
# l5 = [[0., 71405., -1500.], [160965., 71405., -1500.], [160965., 86605., -1500.], [0., 86605., -1500.]]
#
# # Faults
# f1 = [[0., 16732., -4500.], [160965., 16732., -4500.], [160965., 15000., -1500.], [0., 15000., -1500.]]
# f2 = [[0., 16932., -4500.], [160965., 16932., -4500.], [160965., 15200., -1500.], [0., 15200., -1500.]]
# f3 = [[0., 16732., -4500.], [160965., 16732., -4500.], [160965., 16932., -4500.], [0., 16932., -4500.]]
#
# f4 = [[0., 31732., -6000.], [160965., 31732., -6000.], [160965., 30000., -3000.], [0., 30000., -3000.]]
# f5 = [[0., 31932., -6000.], [160965., 31932., -6000.], [160965., 30200., -3000.], [0., 30200., -3000.]]
# f6 = [[0., 31732., -6000.], [160965., 31732., -6000.], [160965., 31932., -6000.], [0., 31932., -6000.]]
#
# f7 = [[0., 54872., -6000.], [160965., 54872., -6000.], [160965., 56605., -3000.], [0., 56605., -3000.]]
# f8 = [[0., 54672., -6000.], [160965., 54672., -6000.], [160965., 56405., -3000.], [0., 56405., -3000.]]
# f9 = [[0., 54672., -6000.], [160965., 54672., -6000.], [160965., 54872., -6000.], [0., 54872., -6000.]]
#
# f10 = [[0., 69872., -4500.], [160965., 69872., -4500.], [160965., 71605., -1500.], [0., 71605., -1500.]]
# f11 = [[0., 69672., -4500.], [160965., 69672., -4500.], [160965., 71405., -1500.], [0., 71405., -1500.]]
# f12 = [[0., 69672., -4500.], [160965., 69672., -4500.], [160965., 69872., -4500.], [0., 69872., -4500.]]
#
bounds = [bot, x_min, x_max, y_min, y_max, l1, l2, l3, l4, l5, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12]

with open(PATH+'TVZ_bdata.txt', 'w+') as out:
    for bound in bounds:
        for point in bound:
            out.write(str(point[0]) + ' ' + str(point[1]) + ' ' + str(point[2]) + '\n')
    out.close()