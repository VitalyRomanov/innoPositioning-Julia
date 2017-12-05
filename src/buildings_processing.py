import overpy
from decimal import *
from math import *
import matplotlib.pyplot as plt
api = overpy.Overpass()

walls_file = open("walls_ways.txt","w+")

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
    # p.lat = resultY
    # p.lon = resultX
    return resultX, resultY


def plotAPloc(nullNode):
	AP_lat = [55.75336]
	AP_lon = [48.74138]
	# AP_lat = [46.716524,46.716344,46.716985,46.716857,46.717188,46.716027,46.714034,46.715187,46.716451,46.717596]
	# AP_lon = [11.652545,11.652767,11.654857,11.654715,11.657764,11.656766,11.660874,11.654817,11.65798,11.658105]

	for i in range(len(AP_lat)):
		ap_node = nn(AP_lat[i],AP_lon[i])
		temp_x,temp_y = getXYpos(nullNode,ap_node)
		print("AP = ", temp_x," | y=", temp_y)
		plt.plot(temp_x,temp_y,'ro')

def writeAreaBorder(nullNode):
	lat = [46.7129922666101,46.7186377333899]
	lon = [11.6499116535865,11.6635073464135]
	b1 = nn(55.75335,48.73943)
	b2 = nn(55.75469,48.73982)
	b3 = nn(55.75441,48.74222)
	b4 = nn(55.75301,48.74163)
	# b1 = nn(46.7129922666101,11.6499116535865)
	# b2 = nn(46.7186377333899,11.6499116535865)
	# b3 = nn(46.7186377333899,11.6635073464135)
	# b4 = nn(46.7129922666101,11.6635073464135)

	temp_x,temp_y = getXYpos(nullNode,b1)
	walls_file.write(repr(temp_x)+","+repr(temp_y)+","+repr(0.0))
	walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(300.0))
	temp_x,temp_y = getXYpos(nullNode,b2)
	walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(300.0))
	walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(0.0)+"\n")

	walls_file.write(repr(temp_x)+","+repr(temp_y)+","+repr(0.0))
	walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(300.0))
	temp_x,temp_y = getXYpos(nullNode,b3)
	walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(300.0))
	walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(0.0)+"\n")

	walls_file.write(repr(temp_x)+","+repr(temp_y)+","+repr(0.0))
	walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(300.0))
	temp_x,temp_y = getXYpos(nullNode,b4)
	walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(300.0))
	walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(0.0)+"\n")

	walls_file.write(repr(temp_x)+","+repr(temp_y)+","+repr(0.0))
	walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(300.0))
	temp_x,temp_y = getXYpos(nullNode,b1)
	walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(300.0))
	walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(0.0)+"\n")



nullNode = nn(55.75335,48.73943)

# writeAreaBorder(nullNode)

print ("Receiving....")

result = api.query("(way(55.75335,48.73943,55.75441,48.74222)[\"building\"~\".\"];node(w););out;")
# result = api.query("(way(46.7129922666101,11.6499116535865,46.7186377333899,11.6635073464135)[\"building\"~\".\"];node(w););out;")

print ("Processing...")
value =0
ind_way = 0
for way in result.ways:
	x = []
	y = []
	print("index way = ", ind_way)
	if (ind_way != 5):
		for ind_node,node in enumerate(way.nodes):
			print("int node = ", ind_node, "value = ", value)
			temp_x,temp_y = getXYpos(nullNode,node)
			# if (value != 59):
			# plt.plot(temp_x, temp_y, 'b')
			# plt.show()
			# value +=1
			x.append(temp_x)
			y.append(temp_y)
			print("x=",x)
			print("y=",y)
			# walls_file.write(repr(temp_x) + " "+repr(temp_y)+ " ")
			value += 1
			if(ind_node):
				walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(2.0))
				walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(0)+"\n")
			walls_file.write(repr(temp_x)+","+repr(temp_y)+","+repr(0))
			walls_file.write(","+repr(temp_x)+","+repr(temp_y)+","+repr(2.0))
		# x.append(x[0])
		# y.append(y[0])
		walls_file.write(","+repr(x[0])+","+repr(y[0])+","+repr(2.0))
		walls_file.write(","+repr(x[0])+","+repr(y[0])+","+repr(0)+"\n")
		# walls_file.write("\n")
		plt.plot(x,y,'b')
		# plt.show()
		print("first loop")
	ind_way+=1
print("plotAPloc(nullNode)")
plotAPloc(nullNode)
print("plt.show()")
plt.show()
print("walls_file.close()")
walls_file.close()
print("finish")

