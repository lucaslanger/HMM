from HMM import *

prior = [0.6, 0.3, 0.1]   # happy,sad,angry

T = [
        [0.8,  0.6, 0.2],
        [0.15, 0.2, 0  ],
        [0.05, 0.2, 0.8]
    ]

O = [(-500,3),(300,10),(1000,25)]

def testSampleGeneration():
    global prior,T,O
    print "Generated Samples: \n\n", generate_samples(prior,T,O)
    
def testStateAssignment():
    global prior,T,O
    counts = [0,0,0]
    for i in range(1000):
        counts[assign_state([0.5,0.3,0.2])] += 1
        
    print "Counts of assignments ", counts
    
def testNormPDF():
    for i in range(10):
        print normPDF(0, 1, 1*i)
        
def testEM():
    EM( generate_samples(prior,T,O)[1][0], debug=True )  #returns 2-d array, hence the 0. The 1 index grabs the observation seqeuence and not the state sequence
    
testEM()
        