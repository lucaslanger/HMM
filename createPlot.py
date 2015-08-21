import matplotlib.pyplot as plt
import sys
import os
import math


numColors = 10
colorsAnalysis = {}
for i in range(numColors):
	r = 0
	g = 100
	b = 0
	colorsAnalysis[i] = '#%02x%02x%02x' % (r,g,b)

colorsTests = {0:'g', 1:'r', 2:'b', 3:'y'}


def generateFont(size):
	font = {'family' : 'serif',
		'color'  : 'black',
		'weight' : 'normal',
		'size'   : size,
		}
	return font

def getDataFromFile(datafile):
	print "Plotting " + datafile
	with open(datafile) as f:
		lineNum = 0
		internalComment = ""
		title = ""
		xaxisLabel = "" 
		yaxisLabel = ""
		xvals = []
		yvals = []
		for line in f.readlines():
			if lineNum == 0:
				try:
					s =  line.split(",")
					xaxisLabel = s[0]
					yaxisLabel = s[1]
				except:
					print "Please include your axis information in the following way: xaxisname,yaxisname"
					print "Failed on the following line:"
					print line
					break
			elif lineNum == 1:
				internalComment = line
			elif lineNum == 2:
				title = line
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
			lineNum = lineNum + 1

		return (xaxisLabel,yaxisLabel,xvals,yvals, title, internalComment)

def takeLogOfArray(a):
	return [math.log(i) for i in a]

def drawPlots(folder, names, colors, alpha_scaling,  verticalPlotSize=-1, horizontalPlotSize=-1):

	l = len(names)	
	if horizontalPlotSize == -1:
		horizontalPlotSize = math.ceil(l**0.5)	
		verticalPlotSize = math.ceil(l**0.5)	
	i = 1
	L = ''
	for name in names:
		n = name[0]
		try:
			d = getDataFromFile(folder+"/" + n)
			plt.subplot(horizontalPlotSize,verticalPlotSize,i)
			title = d[4]
			internalComment = d[5]

			for j in range(len(d[2])):
				xdata = d[2][j]
				ydata = d[3][j]
						
				if name[1] == 'log':
					xdata = takeLogOfArray(xdata)
				if name[2] == 'log':
					ydata = takeLogOfArray(ydata)

				l = len("SingularValues")
				if alpha_scaling and n[:l] != "SingularValues":
					a=(1.0*(j+1))/len(d[2])

				else:
					a = 1.0				
				plt.plot(xdata,ydata, colors[j%len(colors)], alpha=a)

			#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

			plt.xlabel(d[0], fontdict= generateFont(22) )
			plt.ylabel(d[1], fontdict= generateFont(22) )
			plt.title(title, fontdict= generateFont(28) )	

			xcoord = xdata[int(len(xdata)*5/8 )]
			ycoord = ydata[0]*1.5
			print xcoord
			print ycoord
			print internalComment			

			plt.text(xcoord, ycoord, internalComment, fontdict=generateFont(16))

			i += 1
		except Exception,e:
			print str(e)
			print "Trouble plotting ",n  

	plt.subplots_adjust( hspace=0.88 )
	#plt.tight_layout()
	plt.show()

datafile = sys.argv[1]
t = sys.argv[2]
#-v
temp = ['Query_Errors_Base', 'Query_Errors_Naive', 'Comm_Qerror', "QError_Base_vs_Naive",'True_H_vs_Emp', 'True_Ax_vs_Emp',  '(Ax)^2_v.s A(x^2)', 'ConditionalError','ConditionalEmp','ConditionalTrue', 'Base_Errors']
validityTests = [(i,'normal','normal') for i in temp]
#

#-a
modelBased = [('BaseComp_Area', 'log','normal'), ("MinError_Dif_Bases", 'log','normal'),  ("Difference Plot_FIXEDMS",'log','normal'), ("ArgMin_Dif_Bases",'log','normal') , ("Multiple_Trials_ModelError",'normal','normal') , ("Difference Plot",'log','normal')]
#

#-k
baseComparisonKeyPredictions = [("Datasize:256000,64:16,ST:0.0",'normal','normal')]
#[ ('KeyFindingErrorTesting_MaxK:250','normal','normal')]
#[("Datasize:10000PacMan",'normal','normal')]

				#
				#('KeyFindingErrorTraining_MaxK:100','normal','normal')]
				#[("Datasize:10000PacMan",'normal','normal')]

				
#fixedSizeModelBaseComparison = [10000]
#baseComparisonKeyPredictions = []
#for s in fixedSizeModelBaseComparison:
#	baseComparisonKeyPredictions.append(('BaseImprovementOverModelSizesDatasize:' + str(s),'normal','normal'))
#baseComparisonKeyPredictions.append(('SingularValues','normal','normal'))
#

#-mo
multipleObservations = []
dataSizes = [100,1000]
loopPairs = ['16:32','17:27']
for d in dataSizes:
	for l in loopPairs:
		multipleObservations.append( ('errorModelSizesBase,' + str(d) + "," + l,'normal','normal' ) )
		multipleObservations.append(('SingularValues,' + str(d) + "," + l,'normal','normal'))
#

multipleObservations = [('errorModelSizesBase,1000,32:16', 'normal','normal')]

if t=='-v':	
	drawPlots(datafile, validityTests, colorsTests, False)
elif t=='-a':
	drawPlots(datafile, modelBased, colorsAnalysis, True)
elif t=='-k':
	drawPlots(datafile, baseComparisonKeyPredictions, colorsAnalysis, True)	#2,5 to split
elif t=='-mo':
	drawPlots(datafile, multipleObservations, colorsAnalysis, True, 1, 1)
else:
	print "invalid format"
