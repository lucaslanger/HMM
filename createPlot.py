import matplotlib.pyplot as plt
import sys
import os
import math

colors = {}
for i in range(10):
	r = 55 + 20*i
	g = 50
	b = 55 + 20*i
	colors[i] = '#%02x%02x%02x' % (r,g,b)

colorsOther = {0:'g', 1:'r', 2:'b', 3:'y'}

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

def takeLogOfArray(a):
	return [math.log(i) for i in a]

def drawPlots(folder, names):	
	l = len(names)

	i = 1
	s = math.ceil(l**0.5)
	L = ''
	for name in names:
		n = name[0]
		try:
			d = getDataFromFile(folder+"/" + n)
			plt.subplot(s,s,i)
			for j in range(len(d[2])):
				xdata = d[2][j]
				ydata = d[3][j]		
				if name[1] == 'log':
					xdata = takeLogOfArray(xdata)
				if name[2] == 'log':
					ydata = takeLogOfArray(ydata)
				
				plt.plot(xdata,ydata, colorsOther[j%4])

			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

			plt.xlabel(d[0])
			plt.ylabel(d[1])
			plt.title(n)
			i += 1
		except Exception,e:
			print str(e)
			print "Trouble plotting ",n  

	plt.subplots_adjust( hspace=0.88 )
	#plt.tight_layout()
	plt.show()

datafile = sys.argv[1]

temp = ['Query_Errors_Base', 'Query_Errors_Naive', 'Comm_Query_Error', "QError_Base_vs_Naive",'Comm_Matrix_Error','True_H_vs_Emp', 'True_Ax_vs_Emp',  '(Ax)^2_v.s A(x^2)', 'ConditionalError','ConditionalEmp','ConditionalTrue', 'Base_Errors']
initialTests = [(i,'normal','normal') for i in temp]

modelBased = [('BaseComp_Area', 'log','normal'), ("MinError_Dif_Bases", 'log','log'), ("ArgMin_Dif_Bases",'log','normal') , ("Multiple_Trials_ModelError",'normal','normal')]

#drawPlots(datafile, initialTests)
drawPlots(datafile, modelBased)
