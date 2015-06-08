import matplotlib.pyplot as plt
import sys
import os
import math

colors = { 0:"g" ,1:'r', 2:'k', 3: 'b', 4:'c', 5:'m', 6:"y", 7:"y"}

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
	#names = os.listdir(folder)
	names = ['Query_Errors_Base', 'Query_Errors_Naive', 'Non-Comm_Query_Error', "QError_Abs_Base_vs_Naive",'QError_Reg_Base_vs_Naive','Non-Comm_Matrix_Error','True_Hankel_vs_Emperical', 'True_Ax_vs_Emperical_Ax',  '(Ax)^2_v.s A(x^2)']
	l = len(names)

	i = 1
	s = math.ceil(l**0.5)
	L = ''
	for n in names:
		d = getDataFromFile(folder+"/" + n)
		plt.subplot(s,s,i)
		for j in range(len(d[2])):
			if j == 0:
				L = 'Abs'
			else:
				L = 'Non-Abs'
			plt.plot(d[2][j],d[3][j], colors[j], label=L)

		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
        #   ncol=2, mode="expand", borderaxespad=0.)

		plt.xlabel(d[0])
		plt.ylabel(d[1])
		plt.title(n)
		i += 1

	plt.subplots_adjust( hspace=0.50 )
	#plt.tight_layout()
	plt.show()

drawPlots("plotting")