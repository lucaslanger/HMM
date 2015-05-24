import matplotlib.pyplot as plt
import sys

datafile = sys.argv[1]

with open(datafile) as f:
	first = True
	xaxisLabel = "" 
	yaxisLabel = ""
	xaxis = []
	yaxis = []
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
				xaxis.append( float(data[0]) )
				yaxis.append( float(data[1]) )
			except:
				print "Required format not met: x,y"
				print "Failed on the following line:"
				print line
				break


plt.plot(xaxis, yaxis)
plt.xlabel(xaxisLabel)
plt.ylabel(yaxisLabel)
plt.grid(True)
plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.title(datafile)

plt.show()
