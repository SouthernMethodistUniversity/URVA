# This script is designed to do
# fully automatic smoothening function
# for the curves including:
# 1) Scalar curvature plot
# 2) components of scalar curvature plot
# 3) components of path direction
#
# The smoothening consists of 2 major steps:
# 1. Using exponential function to merge two ends
# 2. Using cubit spline fitting to correct central part
#
# Input original curve data
# 1. Old curve whose central region is problematic and needs correct
# 2. New curve who started from the central region, however, needs extra work to
#   2.a Make it merge into the old curve at two ends
#   2.b Correct the central region
#
# Note: This script drops most error handles...
#       This script is expected to incorporate into pURVA.

import sys
from scipy import optimize
import math
from scipy.interpolate import interp1d
import numpy as np
from inter_autosmooth_utils import *
from fileio import expList2File

print ("Info: Entering AutoSmooth interface...")

if len(sys.argv) != 8:
    print ("Error: Number of inputs wrong in autosmooth interface...  ")
    sys.exit()

RL = []
oldcurvf = sys.argv[1]
newcurvf = sys.argv[2]
dstep = float( sys.argv[3] )
RL.append( int(sys.argv[4]) )
RL.append( int(sys.argv[5]) )
continuethresh = float( sys.argv[6] )
outfilename = sys.argv[7]
#


# USER SETTING-UP ##########
#oldcurvf = 'grmpme/oldcurv.dat'
#newcurvf = 'grmpme/newcurv.dat'
#dstep = 0.03
#RL = [3,3]

#continuethresh = 0.5 # Increase this value WHEN NECESSARY

#############################






# ADVANCED SETTING-UP
#continuethresh = 0.4
prederrthresh = 0.1 # Way No. 1

#RL = [1,1]          # Way No. 2 # Chosen as the default approach




# Read in data
olddat=[]
with open(oldcurvf) as oldf:
    for line in oldf:
        if len(line) > 1:
            #print line
            s = float( line.split()[0] )
            v = float( line.split()[1] )
            olddat.append( [s,v] )

newdat=[]
with open(newcurvf) as newf:
    for line in newf:
        if len(line) > 1:
            s = float( line.split()[0] )
            v = float( line.split()[1] )
            newdat.append( [s,v] )


newfor=[]
chkfor=[]
newrev=[]
chkrev=[]

# split forward and reverse
for i in range(len(newdat)):
    if newdat[i][0] >= 0.0:
        newfor.append( newdat[i] )
        chkfor.append(0)
        #print newdat[i]


for i in range(len(newdat)):
    if newdat[i][0] <= 0.0:
        newrev.append( newdat[i] )
        chkrev.append(0)
#      print newrev

# Make sure that new dat and old dat have the same grid point
for i in range(len(newfor)):
    for j in range(len(olddat)):
        if newfor[i][0] == olddat[j][0]:
            chkfor[i] = 1


for i in range(len(newrev)):
    for j in range(len(olddat)):
        if newrev[i][0] == olddat[j][0]:
            chkrev[i] = 1


if sum(chkfor) != len(chkfor):
    print ("Error: Grid point not match...")
    sys.exit()
if sum(chkrev) != len(chkrev):
    print ("Error: Grid point not match...")
    sys.exit()


# Get the label for TS point
TSindex = -1
for i in range(len(olddat)):
    if olddat[i][0] == 0.0:
        TSindex = i
#print TSindex


# Detect problematic region in old curve (maybe optional)

# Scale the new curve starting from the center and expand to two ends
#       Check 2 things:
#          1. The value
#          2. The gradient
# This step should be an iterative procedure, which move step by step
# from transition state point to the end(s)
#



forx=[]
forsign=[]
forratio=[]
forkappa=[]
for i in range(len(newfor)-1):
    #print "i",i
    j = i + 1 # Current point
    p = i     # Point before current point
    sj = newfor[j][0]
    sp = newfor[p][0]
    vj = newfor[j][1]
    vp = newfor[p][1]
    #print sp, sj
    #print vp, vj
    #
    #
    # p  -> previous point
    # j  -> current point
    # q  -> past point
    #
    #
    oj = olddat[ TSindex + j ][1] # Current point value in old curve
    q = TSindex + j + 1 # The point past currect point in old curve
    sq = olddat[q][0]
    vq = olddat[q][1]
    #print sq
    Snew = sj
    Vnew = vj
    Vold = oj
    #
    #
    # "diffe" here is a user-defined function
    # x -> variable will be changed by optimizer
    # Snew, Vnew, Vold and dstep -> parameters that are fixed during optimization
    #
    #
    res = optimize.minimize_scalar(diffe,method="Bounded",bounds=(1e-2,1e10),args=(Snew,Vnew,Vold,dstep))
    unknownx =  res.x    # The "x" after optimization

    forx.append(unknownx)
    if abs(vj) > abs(oj):
        forsign.append(-1)
    else:
        forsign.append(1)

    # check derivative then
    #print unknownx
    #print sj, unknownx
    #print sp,sj,sq,'-',vp,vj,"-",oj,vq,"-",dstep
    #forratio.append(

    # " unknownx " is used here

    tmpresult =  graddiff(unknownx,sp,sj,sq,vp,vj,oj,vq,dstep)

    #print tmpresult
    forratio.append( tmpresult[0] )
    forkappa.append( tmpresult[1] ) # This is curvature value

#print forratio
#print forkappa

# Re-order newrev list
tmp=[]
for i in range(len(newrev)):
    j = len(newrev) - 1 - i
    tmp.append( newrev[j] )
newrev = tmp


revratio=[]
revkappa=[]
revx=[]
revsign=[]
for i in range(len(newrev)-1):
    j = i + 1
    p = i
    sj = newrev[j][0]
    sp = newrev[p][0]
    vj = newrev[j][1]
    vp = newrev[p][1]

    oj = olddat[ TSindex - j ][1]
    q = TSindex - j - 1
    sq = olddat[q][0]
    vq = olddat[q][1]

#    print sq,sj,sp

    Snew = sj
    Vnew = vj
    Vold = oj
    res = optimize.minimize_scalar(diffe,method="Bounded",bounds=(1e-2,1e10),args=(Snew,Vnew,Vold,dstep))
    unknownx = res.x

    revx.append(unknownx)
    if abs(vj) > abs(oj):
        revsign.append(-1)
    else:
        revsign.append(1)

#   print unknownx
    #revratio.append( graddiff(unknownx,sp,sj,sq,vp,vj,oj,vq,dstep) )
    tmpresult = graddiff(unknownx,sp,sj,sq,vp,vj,oj,vq,dstep)
    revratio.append( tmpresult[0] )
    revkappa.append( tmpresult[1] )



#print revratio

# Checking gradient and then determine when to cut
## input: forratio, revratio

#print forratio
#print revratio

#detercut(forratio)


# detercut is based on the gradient ratio which is not so robust
 #fn = detercut(forratio,continuethresh)
 #rn = detercut(revratio,continuethresh)





# So we need to use curvature instead, which is implemented as detercutB
fn = detercutB(forkappa,continuethresh)
rn = detercutB(revkappa,continuethresh)





print(("Info: In forward direction, "+str(fn)+" points need correction."))
print(("Info :In reverse direction, "+str(rn)+" points need correction."))


#########################################################################
#########################################################################
# Next step is to smooth back the central region.
# Another iterative procedure is expected.
#

# calculate the scaling factor and scaled values
# input: forx, revx, forsign, revsign

#print forx,revx,forsign,revsign


# Get the last point information
forx = forx[ fn-1 ]
forsign = forsign[ fn-1 ]
revx = revx[ rn-1 ]
revsign = revsign[ rn-1 ]

#print forx,forsign
#print revx,revsign

print(("Info: The scaling factor for forward direction is: exp("+str(forsign)+"*abs(\"s\")/"+str(dstep)+"/"+str(round(forx,2))+")"))
print(("Info: The scaling factor for reverse direction is: exp("+str(revsign)+"*abs(\"s\")/"+str(dstep)+"/"+str(round(revx,2))+")"))


#


forsv = []  # s value
forkv = []  # scaling factor
forvv = []  # old value
forvscale = [] # scaled value

revsv = []
revkv = []
revvv = []
revvscale = []

for i in range( fn ):
    forsv.append( newfor[i+1][0] )
    forvv.append( newfor[i+1][1] )
    k = math.exp( forsign* abs(newfor[i+1][0])/dstep/forx )
    forkv.append( k )
    forvscale.append( k* newfor[i+1][1] )

for i in range( rn ):
    revsv.append( newrev[i+1][0] )
    revvv.append( newrev[i+1][1] )
    k = math.exp( revsign* abs(newrev[i+1][0])/dstep/revx )
    revkv.append( k )
    revvscale.append( k* newrev[i+1][1] )


#print revvscale

formark = -1
revmark = -1

for i in range(len(olddat)):
    if olddat[i][0] == forsv[-1]:
        formark = 0
        continue
    if formark == 0:
        forp1s = olddat[i][0]
        forp1v = olddat[i][1]
        formark = formark + 1
        continue
    if formark == 1:
        forp2s = olddat[i][0]
        forp2v = olddat[i][1]
        formark = formark + 1
        continue

# Get the information for 2 extra points
# These 2 extra points are used to help the cubic spline fitting

#print forp1s,forp2s


for i in range(len(olddat)):

    j = len(olddat) - i - 1
    if olddat[j][0] == revsv[-1]:
        revmark = 0
        continue
    if revmark == 0:
        revp1s = olddat[j][0]
        revp1v = olddat[j][1]
        revmark = revmark + 1
        continue
    if revmark == 1:
        revp2s = olddat[j][0]
        revp2v = olddat[j][1]
        revmark = revmark + 1
        continue

# Get the information for 2 extra points
# These 2 extra points are used to help the cubic spline fitting


#print revp1s,revp2s

#
# 'forsf' is the data list for spline fitting
#

forsf = []
# Put 2 extra points in reverse direction into the list
forsf.append( [ revp2s, revp2v ] )
forsf.append( [ revp1s, revp1v ] )

#print revsv

# Put scaled values in reverse direction into the list
for i in range(len(revsv)):
    j = len(revsv) - i - 1
    forsf.append( [ revsv[j], revvscale[j] ] )  # bug!


#print forsf

# Put s = 0.0 into the list
forsf.append([ 0.0, 0.0] )

# Put scaled values in forward direction into the list
for i in range(len(forsv)):
    forsf.append( [ forsv[i], forvscale[i] ] ) # bug !

# Put 2 extra points in forward direction into the list
forsf.append( [ forp1s, forp1v ] )
forsf.append( [ forp2s, forp2v ] )

#print forsf # This is corrected data before spline fitting (4 extra original data included)


# Till now, the data for spline fitting is ready
################################################

##########################
# | | | | | | | | | | | |
# > > > > > > > > > > > >

marks = []
sfmark= []
for i in range(len(forsf)):
    marks.append(1)
    sfmark.append(0)

marks[0] =0
marks[1] =0
marks[2] =0
marks[-1]=0
marks[-2]=0
marks[-3]=0


# < < < < < < < < < < < <
# | | | | | | | | | | | |
#########################








#print sum(marks)


# Spline fitting procedure
# This is an iterative process
# input:
#  1. data for spline fitting -> forsf
#  2. label -> marks
#

#def spline1dlist(lis1):

#def get1point(data,marks):

#def get2points(data,marks):
#    # This function returns the two points to be predicted by spline fitting

#def checkerr(f,pts):
#    # check the error between predicted value and train value

#def adjustmark(marks,l,r):
#    # this function could make 1 into 0

#def updateval( err,thresh,data,pts,f ):

## SPLINE FITTING
#def itersplif(data,marks):

#print forsf

#################################################################
#forsf = itersplif(forsf,marks) # THIS IS NOT WORKING PROPERLY!!!
#################################################################



# Another way to automatic spline fitting this region is to
# use a quadratic function try to fit N-3, N-5, N-7, ... points
# until the residue error is below a threshold ...



##########################
# | | | | | | | | | | | |
# > > > > > > > > > > > >

quadgood = 0
expand = 0

while quadgood == 0:
    for i in range(len(forsf)):
        if forsf[i][0] == 0.0:
            TSindex = i
    #
    expand = expand + 1
    quadgood = 1


# < < < < < < < < < < < <
# | | | | | | | | | | | |
#########################

##########################
# | | | | | | | | | | | |
# > > > > > > > > > > > >
#  This part is useless #
# < < < < < < < < < < < <
# | | | | | | | | | | | |
#########################






#
# Most simple way to do the automatic spline fitting
# is do tell the program how many points need to be
# ruled out for spline fitting on the left side and
# and the right side.
# The user need to give this number !
#

#def semispline( forsf, RL ):


forsf = semispline( forsf,RL)


#print forsf

####for i in range(len(forsf)):
####    print round(forsf[i][0],4),"        ",round(forsf[i][1],4)

expList2File(forsf,outfilename)
