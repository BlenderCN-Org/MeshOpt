import googlemaps
from time import sleep
import numpy
from math import sin, cos, sqrt, atan2, radians, degrees, asin, tan


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


def get_vector(point_1, point_2):

    lat1 = radians(point_1[0])
    lon1 = radians(point_1[1])
    lat2 = radians(point_2[0])
    lon2 = radians(point_2[1])
    var_x = sin(lon2-lon1) * cos(lat2)
    var_y = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon2-lon1)
    return (degrees(atan2(var_y, var_x))+360)%360


def get_intermediate_point(ang_dist, pos, point_1, point_2):

    lat1 = radians(point_1[0])
    lon1 = radians(point_1[1])
    lat2 = radians(point_2[0])
    lon2 = radians(point_2[1])
    a = sin((1-pos)*ang_dist) / sin(ang_dist)
    b = sin(pos * ang_dist) / sin(ang_dist)
    x = a * cos(lat1) * cos(lon1) + b * cos(lat2) * cos(lon2)
    y = a * cos(lat1) * sin(lon1) + b * cos(lat2) * sin(lon2)
    z = a * sin(lat1) + b * sin(lat2)
    rad_lat, rad_lon = atan2(z, numpy.sqrt(x**2 + y**2)), atan2(y, x)
    return [degrees(rad_lat), degrees(rad_lon)]


def get_end_point(start, bearing, dist):

    ang_dist = dist/6371.0
    s_lat = radians(start[0])
    s_lon = radians(start[1])
    f_lat = asin(sin(s_lat)*cos(ang_dist)+cos(s_lat)*sin(ang_dist)*cos(bearing))
    return [degrees(f_lat), degrees(s_lon + atan2(sin(bearing)*sin(ang_dist)*cos(s_lat),cos(ang_dist)-sin(s_lat)*sin(f_lat)))]


def random_chunking(temp_list, sub_size):
    nb_points = len(temp_list)
    temp_matrix = numpy.array([temp_list[i:i + sub_size] for i in xrange(0, nb_points, sub_size)])
    return temp_matrix


def build_matrix_request(point1, point2, point3, target_points):

    NS_dist = calculate_lenght(point1[0], point1[1], point3[0], point3[1])
    EW_dist = calculate_lenght(point1[0], point1[1], point2[0], point2[1])
    dist_ratio = NS_dist / EW_dist
    point_NS = int(numpy.sqrt(target_points / dist_ratio))
    point_EW = int(dist_ratio * point_NS)

    EW_bearing = get_vector(point1, point2)
    NS_sources = [get_intermediate_point(NS_dist/6371.0, float(i)/(point_NS-1), point1, point3) for i in xrange(point_NS)]
    grid = numpy.zeros((point_NS*point_EW,2))
    counter = 0
    for source in NS_sources:
        end_point = get_end_point(source, EW_bearing, EW_dist)
        for j in xrange(point_EW):
            grid[counter] = get_intermediate_point(EW_dist/6371.0, float(j)/(point_EW-1), source, end_point)
            counter += 1

    temp_mat = random_chunking(grid, 512)

    return random_chunking(temp_mat, 50)

if __name__ == "__main__":

    GMS_API_KEY = 'AIzaSyB3q6dqfR9gomJPgtwKTR_A1Se41Lhnzbc'
    ELEVATION_BASE_URL = 'https://maps.googleapis.com/maps/api/elevation/xml'
    TOPO_PATH = '/home/lmar626/Documents/'

    NW_vertice = [-42.5, 169.5] # Max lat; Min long
    NE_vertice = [-42.5, 171.5] # Max lat; Max long
    SE_vertice = [-44, 171.5] # Min lat; Max long
    nb_points = 100000

    NS_dist = calculate_lenght(NE_vertice[0], NE_vertice[1], SE_vertice[0], SE_vertice[1])
    EW_dist = calculate_lenght(NE_vertice[0], NE_vertice[1], NW_vertice[0], NW_vertice[1])
    dist_ratio = NS_dist / EW_dist
    point_EW = int(numpy.sqrt(nb_points / dist_ratio))
    point_NS = int(dist_ratio * point_EW)
    print point_EW, point_NS
    matrix_request = build_matrix_request(NE_vertice, NW_vertice, SE_vertice, nb_points)

    gmaps = googlemaps.Client(key=GMS_API_KEY)

    bulk_results = []
    for step in matrix_request:
        for request in step:
            bulk_results.append(gmaps.elevation(request))
        sleep(1)

    f = open(TOPO_PATH + 'box_100000.txt', 'w+')
    refine_results = []
    for res in bulk_results:
        for single_res in res:
            refine_results.append([single_res[u'location'][u'lat'],single_res[u'location'][u'lng'],single_res[u'elevation']])
            f.write(str(single_res[u'location'][u'lat']) + ' ' + str(single_res[u'location'][u'lng']) + ' ' +str(single_res[u'elevation']) + '\n')

    f.close()



