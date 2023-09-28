# This is an interface of " curcorbra " for pURVA
# This module could be imported by pURVA directly 
# However, it relies on two dependent files:
#  1. interf-curvcor-util.py  
#  2. interf-curvcor-util2.py
#
# INPUT FROM pURVA:
#  1. List of objects
#  2. Number of points for correction on the left
#  3. Number of points for correction on the right
#
# OUTPUT to pURVA:
#  A list that contains s values along with corrected curvature value
#
#
# CALLED IN pURVA:
#   corcurv(oBlist,Ln,Rn)
#


from calc import solWils
from fileio import stop 
from util import LT2Sqr
import gc
import numpy as np

def corcurv(oBlist,Ln,Rn):
     
     print ("Info: Entering CurvCor interface...")	
     # 
     for i in range(len(oBlist)):
	 if oBlist[i].s == 0.0:
	     #print "zero "
             zeroindex = i
     points = []
     for i in range(Ln):
         j = Ln - i
	 #print zeroindex - j
         points.append( zeroindex - j )
     #print zeroindex
     points.append( zeroindex )
     for i in range(Rn):
	 j = i + 1
	 points.append( zeroindex + j )

     #

     corkappas = []
     #print "points",points
     #points = [563]
     #print tripleprocess(zeroindex,oBlist)
     for i in range(len(points)):
         sval,kappaval =  tripleprocess( points[i],oBlist )
         corkappas.append([sval,kappaval])

     #print corkappas
     return corkappas


def extmcgh(pindex,oBlist):
    # Extract mass, coordinates, gradients and hessian
    # for a specific point

    Vmass = oBlist[pindex].Vmass
    Vcoor = oBlist[pindex].Vcoor
    Vgrad = oBlist[pindex].Vgrad
    Vhess = oBlist[pindex].Vhess

    return (Vmass, Vcoor, Vgrad, Vhess )


def tripleprocess(pindex,oBlist):
    a = pindex - 1 
    b = pindex 
    c = pindex + 1
    
    #print oBlist[b].Vhess

    Vmass, VcoorA, VgradA, VhessA = extmcgh( a, oBlist )
    Vmass, VcoorB, VgradB, VhessB = extmcgh( b, oBlist )
    Vmass, VcoorC, VgradC, VhessC = extmcgh( c, oBlist )

    # Get the normal modes
    eig,nm,cdisp=solWils(VcoorB,VhessB,Vmass) 

    # Flatten Vcoor
    VcoorA = [item for sublist in VcoorA for item in sublist]
    VcoorB = [item for sublist in VcoorB for item in sublist]
    VcoorC = [item for sublist in VcoorC for item in sublist]

    # Re-calculate step size of IRC
    s = 0.0
    for i in range(len(VcoorB)):
	    j = i / 3
	    tmp = ((VcoorA[i] - VcoorB[i])**2 ) * Vmass[j]
	    s = s + tmp
    dsn = s**0.5

    s = 0.0 # fixed here
    for i in range(len(VcoorA)):
            j = i / 3 
            tmp = ((VcoorB[i] - VcoorC[i])**2 ) * Vmass[j]
	    s = s + tmp
    dsp = s**0.5

    #print dsp,dsn
    ## Complain if dsn or dsp goes wrong 
    if dsn == 0.0 or dsp == 0.0:
       stop("Error: dsn/dsp in interface curvcor... ")
    if (dsn/dsp)>= 5.0 or (dsp/dsn)>= 5.0:
       stop("Error: dsn/dsp huge difference in interface curvcor...")

    # Calculation 
    NAtom = len(Vmass)
    ffM1 =  LT2Sqr( 3*NAtom, VhessB )
    ffMp =  LT2Sqr( 3*NAtom, VhessC )
    ffMn =  LT2Sqr( 3*NAtom, VhessA )

    del VhessA
    gc.collect()
    del VhessC
    gc.collect()


    for i in range(3*NAtom):
        for j in range(3*NAtom):
            im = i / 3
	    jm = j / 3
            ffMp[i][j] = ffMp[i][j] / ( (Vmass[im])**0.5 * (Vmass[jm]**0.5) )
            ffMn[i][j] = ffMn[i][j] / ( (Vmass[im])**0.5 * (Vmass[jm]**0.5) )
            ffM1[i][j] = ffM1[i][j] / ( (Vmass[im])**0.5 * (Vmass[jm]**0.5) )

    #
    ffMp = np.matrix( ffMp )
    ffMn = np.matrix( ffMn )
    #print dsn, dsp
    mF1 = (ffMp - ffMn ) / (dsn+dsp)
    del ffMn
    del ffMp
    gc.collect()

   
    ffM1 = np.matrix(ffM1)
    Veta = np.array( nm[0]  )

    tmp = (Veta * ffM1 )
    tmp = np.asarray( tmp )
    a = np.dot(tmp, Veta)
    del tmp
    gc.collect()

    I = np.identity( 3*NAtom )
    I = I * 2 * a

    Left = np.matrix(I) - ffM1
    Leftinv = np.linalg.inv( Left )

    del Left
    gc.collect()

    tmp = mF1.dot(Veta)
    del mF1
    gc.collect()
    tmp = np.asarray(tmp)[0]

    a = np.dot(Veta,tmp)
    
    av = Veta * a

    Right = tmp - av

    del tmp
    del av
    gc.collect()

    cur = Leftinv.dot(Right)
    tt = cur.tolist()

    result = np.linalg.norm(cur)
    sval = oBlist[b].s
    
    return (sval,result)
   





    #print "hello...."








