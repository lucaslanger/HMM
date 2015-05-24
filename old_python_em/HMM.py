#This model has observed states of purchases with the hidden states being whether the customer is {1:Happy, 2:Sad, 3:Angry}

import numpy
import math
from numpy import matrix
from numpy import linalg
from numpy import random



#input prior probability, transition matrix 3x3, Observation gaussians 3x2
def generate_samples(prior, T, O, number_of_samples = 1, duration = 500):       
    sample_seq = []
    observation_seq = []
    for i in range(number_of_samples):
        sample = []
        observation = []
        
        state = assign_state(prior)
        sample.append(state)
            
        for j in range(duration):
            state = assign_state( T[state] )
            sample.append(state)
            observation.append( generate_observation(state,O) )
            
        sample_seq.append(sample)
        observation_seq.append(observation)
            
    return sample_seq, observation_seq
        
def generate_observation(state, O):
    return numpy.random.normal( O[state][0], O[state][1] )

def assign_state(p_a):
    r = numpy.random.random()
    if r < p_a[0] :
        return 0
    elif r < p_a[0] + p_a[1]:
        return 1
    else:
        return 2

def EM( observations , num_states=3 , iterations = 10, debug = False):
    assignments = {0:[], 1:[], 2:[] }
    stats = {0:[0,0],1:[0,0],2:[0,0]}
    
    for o in observations:                            #Initial random guess
        a = assign_state([0.333,0.333,0.334])
        assignments[a].append(o)
    
    for i in range(iterations):
        
        #Maximize
        for key in assignments:
            l = len( assignments[key] )
            if l!= 0:
                stats[key][0] = sum( assignments[key] ) / l                                                    #Mean
                stats[key][1] = computeSD(assignments[key], stats[key][0] ) #Standard Deviation
            else:
                print "No items in group ", key
        #Expect         Note: The standard deviation is important!!! can't just use distance from the mean because even though it might be closer to some point, another point might have larger sigma --> more likely to be part of this group  
        assignments = {0:[], 1:[], 2:[] }
        for o in observations:
            m = max( [normPDF( stats[i][0],stats[i][1],o ) for i in stats] )    # compute max probability density 
            
            a = None
            for i in stats:
                 if normPDF( stats[i][0],stats[i][1],o ) == m:                  # find which group (in 0,1,2) produces this maximal values m
                    a = i
                    assignments[a].append(o)
                    break
                 
            if a == None:
                print "Disaster in assignment"
                
        if debug:
            print "Models:" , stats
            print "Assignment counts (size of groups)" , [len(assignments[s]) for s in assignments]
    
    return stats, assignments

def computeSD(array, mean): #array, arraylength, mean
    l = len(array) 
    return (sum( [ (array[i] - mean)**2 for i in range(l) ] )/ l)**0.5

def normPDF(mean, sd, datapoint):
    var = float(sd)**2
    pi = 3.1415926
    denom = (2*pi*var)**.5
    num = math.exp(-(float(datapoint)-float(mean))**2/(2*var))
    return num/denom
