# import matplotlib.pyplot as plt
from math import *

walls_file = open("walls.txt","r")


class nn:
	def __init__(self,lat,lon):
		self.lat = lat
		self.lon = lon

def asRadians(degrees):
    return degrees * pi / 180

def getXYpos(relativeNullPoint, p):
    """ Calculates X and Y distances in meters.
    """
    deltaLatitude = float(p.lat) - relativeNullPoint.lat
    deltaLongitude = float(p.lon) - relativeNullPoint.lon
    latitudeCircumference = 40075160 * cos(asRadians(relativeNullPoint.lat))
    resultX = deltaLongitude * latitudeCircumference / 360
    resultY = deltaLatitude * 40008000 / 360
    return resultX, resultY



def plotAPloc(nullNode):
	AP_lat = [46.716524,46.716344,46.716985,46.716857,46.717188,46.716027,46.714034,46.715187,46.716451,46.717596]
	AP_lon = [11.652545,11.652767,11.654857,11.654715,11.657764,11.656766,11.660874,11.654817,11.65798,11.658105]

	for i in range(len(AP_lat)):
		ap_node = nn(AP_lat[i],AP_lon[i])
		temp_x,temp_y = getXYpos(nullNode,ap_node)
		plt.plot(temp_x,temp_y,'ro')


nullNode = nn(46.7129922666101,11.6499116535865)

for line in walls_file:
    coord = line.strip().split(",")
    # print coord[0],coord[1],coord[3],coord[4],coord[6],coord[7],coord[9],coord[10]
    plt.plot([coord[0],coord[6]],[coord[1],coord[7]],'b')

plotAPloc(nullNode)
plt.show()
