import matplotlib.pyplot as plt
import sys
import os

def getDataFromFile(datafile):
	with open(datafile) as f:
		first = True
		xaxisLabel = "" 
		yaxisLabel = ""
		xvals = []
		yvals = []
		for line in f.readlines():
			if first:
				try:
					s =  line.split(",")
					xaxisLabel = s[0]
					yaxisLabel = s[1]
					first = False
				except:
					print "Please include your axis information in the following way: xaxisname,yaxisname"
					print "Failed on the following line:"
					print line
					break
			else:
				try:
					data = line.split("\n")[0].split(",")
					xvals.append( float(data[0]) )
					yvals.append( float(data[1]) )
				except:
					print "Required format not met: x,y"
					print "Failed on the following line:"
					print line
					break

		return (xaxisLabel,yaxisLabel,xvals,yvals)

def drawPlots(folder):
	names = os.listdir(folder)
	l = len(names)
	s = int(l**0.5)

	i = 1
	for n in names:
		d = getDataFromFile(folder+"/" + n)
		#plt.subplot(s,s,i)
		plt.subplot(1,l,i)
		plt.plot(d[2],d[3], "b")
		plt.xlabel(d[0])
		plt.ylabel(d[1])
		plt.title(n)
		i += 1

	#plt.subplots_adjust(wspace=0.4, hspace=0.02, top=.9, bottom=0.02, left=0.02, right=0.98)
	plt.show()

drawPlots("plotting")