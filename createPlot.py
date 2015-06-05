import matplotlib.pyplot as plt
import sys
import os

colors = { 0:"r" ,1:'g', 2:'k', 3: 'b', 4:'c', 5:'m', 6:"y", 7:"y"}

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
					data = (line.split("/n")[0])[:-1].split(" ")[:-1]

					for d in range(len(data)):
						s = data[d].split(",")
						if len(xvals) <= d:
							xvals.append([float(s[0])])
							yvals.append([float(s[1])])
						else:
							xvals[d].append(float(s[0]) )
							yvals[d].append(float(s[1]) )  

				except Exception, e:
					print "Required format not met"					
					print str(e)
					print s, datafile
					break

		return (xaxisLabel,yaxisLabel,xvals,yvals)

def drawPlots(folder):
	names = os.listdir(folder)
	l = len(names)

	i = 1
	for n in names:
		d = getDataFromFile(folder+"/" + n)
		print d[2]
		plt.subplot(1,l,i)
		for j in range(len(d[2])):
			plt.plot(d[2][j],d[3][j], colors[j], label='line' + str(j))

		plt.xlabel(d[0])
		plt.ylabel(d[1])
		plt.title(n)
		i += 1

	#plt.subplots_adjust(wspace=0.4, hspace=0.02, top=.9, bottom=0.02, left=0.02, right=0.98)
	#plt.tight_layout()
	plt.show()

drawPlots("plotting")