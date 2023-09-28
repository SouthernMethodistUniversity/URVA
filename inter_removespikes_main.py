# This script is designed to remove the
# spikes in the curvature plots.
#
# The removal of spikes is based on the following:
# 1. a fix threshold ( if the value is above this threshold, cut it. )
# 2. Get a distribution of difference in value between two steps,
#    check into the highest 80% or 90%
#   2.b check peripheral gradient, if the gradient of the point in question
#       is too different from peripheral gradient, cut this point.
# 3. cut out the two ends first
#


from inter_removespikes_utils import *
import sys
from fileio import expList2File

print ("Info: Entering RemoveSpikes interface...")

if len(sys.argv) != 6:
    print ("Error: Number of parameters wrong for \'RemoveSpikes\' interface...  ")

filename = sys.argv[1]
cuthigh = float( sys.argv[2]  )
percentage = float( sys.argv[3] )
gradientratio = float( sys.argv[4]  )
outfilename = sys.argv[5]


#filename = "data.txt"
#cuthigh = 20.0
#percentage = 0.85
#gradientratio = 1.2
#outfilename = "output.dat"



#print ("filename: "+filename)
#print ("[Parameter] cuthigh: "+str(cuthigh))
#print ("[Parameter] percentage: "+str(percentage))
#print ("[Parameter] gradientratio: "+str(gradientratio))



s,curv = readdat(filename)

# check whether we have duplicate points
todelete = checkdupli(s)

if len(todelete) != 0:
    s1 = []
    curv1 = []

    for i in range(len(s)):
        if i in todelete:
            pass
        else:
            s1.append( s[i] )
            curv1.append( curv[i] )

    s = s1
    curv = curv1
    # CAUTION

# calculate the difference of each point with regard to two points nearby
difs = calcdiff(curv)
#print difs

#for i in difs:
#    print i


# First filter -> possible spikes are 1 in the "labels" list

labels = checkpercent( s,curv,difs,percentage )


# Second filter -> drop some wrong spike points by checking the gradient

labels = gradientfilter(s,curv,labels,gradientratio)


# Third filter -> cut down the points above a threshold

labels = directcut(curv, labels, cuthigh)





outvalues=[]

for i in range(len(labels)):
    if labels[i] == 0:
        #print s[i], "   ", curv[i]
        outvalues.append( [ s[i],  curv[i] ] )

expList2File( outvalues, outfilename )
