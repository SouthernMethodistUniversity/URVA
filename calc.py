import numpy as np
import scipy
from sympy import *
import math
import sys
from util import *
import copy
import gc
import os
from fileio import ProgressBar
from fileio import stop
from fileio import rdkeywrds
from cologne import DMODriv
from cologne import ABscIn




## Calculate vibrational frequencies and normal modes for each point on the path
## And also reorder them using DMO

def calfrq_dmo(Brlist,excludelist,tfcurvcpl,tfcoriolis,inpfile,freqfile):
    # Brlist is usually oBlist
    # excludelist contains the information of discontinuities

    #-------------------------------------
    # Is this overlap threshold too high ?
    Sthresh0 = 0.990
    Slowest0 = 0.890
    # # # # # # # #
    Np0 = 30
    NMax0 = 200 #500
    Cut0 = 0
    Skip0 = 0
    #-------------------------------------

    # Read DMO parameters from user input file
    dmosetlist= [0]

    Skip = rdkeywrds(inpfile,"DMO","Skip")
    if Skip == "toberead":
        Skip = Skip0
    if Skip == 1:
        SkipA = rdkeywrds(inpfile,"DMO","SkipA")
        SkipB = rdkeywrds(inpfile,"DMO","SkipB")


    Cut = rdkeywrds(inpfile,"DMO","Cut")
    if Cut == "toberead":
        Cut = Cut0
    #print "Cut = ",Cut
    if Cut == 1:
        CutA = rdkeywrds(inpfile,"DMO","CutA")
        CutB = rdkeywrds(inpfile,"DMO","CutB")


    Sthresh = rdkeywrds(inpfile,"DMO","Sthresh")
    if Sthresh == 'toberead':
        Sthresh = Sthresh0
    if Sthresh < 0.5 or Sthresh > 1.0:
        stop("Error: Please provide a reasonable value of Sthresh in DMO")


    Slowest = rdkeywrds(inpfile,"DMO","Slowest")
    if Slowest == 'toberead':
        Slowest = Slowest0
    if Slowest >= Sthresh:
        stop("Error: Slowest should be smaller than Sthresh")
    if Slowest < 0.5 or Slowest > 1.0:
        stop("Error: Please provide a reasonable value of Slowest in DMO")


    Np = rdkeywrds(inpfile,"DMO","Np")
    if Np == 'toberead':
        Np = Np0
    if Np < 5:
        stop("Error: Please provide a reasonable value of Np in DMO")

    NMax = rdkeywrds(inpfile,"DMO","NMax")
    if NMax == 'toberead':
        NMax = NMax0
    if NMax <= Np:
        stop("Error: Please provide a reasonable value of NMax in DMO")





    #--------------------------------------------------------------------
    #
    # Here we need to add one more option in which only a specific range of s
    # has vibrational frequencies calcuated and ordered
    # input: cutA, cutB

    skiplist = []
    if Skip ==1:
        skiplist.append(0)  # TS point skipped
        if SkipA == "toberead":
            SkipA = 0
        if SkipB == "toberead" :
            SkipB = 0
    else:
        SkipA = 0
        SkipB = 0
    if SkipA >= 2:
        print("Warning: SkipA is over 2")
    if SkipB >= 2:
        print("Warning: SkipB is over 2")


    if Cut == 1:
        if CutA == "toberead":
            CutA = Brlist[0].s  - 0.01
        if CutB == "toberead":
            CutB = Brlist[-1].s + 0.01
        cutA =  CutA
        cutB =  CutB
        #print "cutA",cutA, "cutB",cutB
    else:
        cutA = Brlist[0].s  - 0.01 # Full length
        cutB = Brlist[-1].s + 0.01 # Full length


    numPn = len(Brlist)

    klib = []
    k = "no"  # index for TS point
    for i in range(numPn):
        if Brlist[i].s == 0.0:
            k = i
            klib = [ k-2, k-1, k, k+1, k+2 ]

    # skip
    for i in range(SkipB):
        skiplist.append( i+1 )
    for i in range(SkipA):
        skiplist.append( -1*(i+1) )
    for i in range(len(skiplist)):
        skiplist[i] = skiplist[i] + k




    pool = [] # store the index of points in Brlist

#   numPn = len(Brlist)

    for i in range(numPn):
        if i not in skiplist:
            if Brlist[i].s > cutA and Brlist[i].s < cutB:   #
                pool.append( i )

    if len(pool) < 3:
        stop("Error: Too few points on reaction path used for DMO. ")


#    klib = []
#    k = "no"  # index for TS point
#    for i in range(numPn):
#       if Brlist[i].s == 0.0:
#          k = i
#          klib = [k]#[ k-2, k-1, k, k+1, k+2 ]


    # First point as the starting step
    if pool[0] in klib:            # k == pool[0]:  # corrected for klib
        # if the first point is TS point
        e0,nm0,cdisp0 = solWilsProj( Brlist[pool[0]].Vcoor,Brlist[pool[0]].Vhess,Brlist[pool[0]].Vmass,Brlist[pool[0]].Vnm_0 )
    else:
        e0,nm0,cdisp0 = solWilsProj( Brlist[pool[0]].Vcoor,Brlist[pool[0]].Vhess,Brlist[pool[0]].Vmass,Brlist[pool[0]].Veta )
    #prtVib(e0) # <- first point in Browsing file

    # paused on May 15th, 2017


    if os.path.isfile(freqfile):
        stop("Error: this file already exists- "+freqfile)

    f1h = open(freqfile,"w")
    nfreq = len(e0)
    f1h.write("   s ")
    for l in range(nfreq):
        f1h.write("vib"+str(l+1)+" ")
    f1h.write("\n")





    # Loop over the rest points
    nvib = len(nm0)
    natom3 = len(nm0[0])
    #print nvib,natom3

    bar = ProgressBar(total = (len(pool)-1) )

    for i in range(len(pool)-1):

        bar.move()
        bar.log('We have arrived at: ' + str(i ))

        #j = pool[i] + 1
        j = pool[i+1]
        if i == 0:
            nm1 = copy.deepcopy( nm0 )
        else:
            nm1 = copy.deepcopy( nm2 )

        if j in klib:   # j == k: # k is for TS  # corrected for klib
            e2,nm2,cdisp2 = solWilsProj( Brlist[j].Vcoor,Brlist[j].Vhess,Brlist[j].Vmass,Brlist[j].Vnm_0 )
        else:
            e2,nm2,cdisp2 = solWilsProj( Brlist[j].Vcoor,Brlist[j].Vhess,Brlist[j].Vmass,Brlist[j].Veta  )

        #prtVib(e2)  # <- debug use

        # Stop before DMO if discontinuity is encountered
        ## at this moment, we just stop the program. Later, a solution will be figured out and tested...

        if len(excludelist) != 0:
            for m in range(len(excludelist)):
                if excludelist[m][1] == j:
                    #stop("Error: DMO module crashed due to discontinuity...")
                    pass


        nm2,e2 = DMODriv(nm1,nm2,e2,natom3,nvib)


        # Check overlap after 1st DMO run
        ov_test =  overlap_dmo(nm1,nm2)  # the min overlap value is returned


        #print Brlist[j].s,  round(ov_test,4)
        #if ov_test < (nvib-1):   # This part seems useless.....
                                  # However, later I realize that this could help to identify the local difficulty


        # Identify local difficulty via overlap value
        # try to solve local difficulty with reduced stepsize
        if ov_test < Sthresh:
            print(( "Warning: Local difficulty encountered at "+str(Brlist[j].s) ))

            e2,nm2 = solvLocDif(nm1, pool[i] ,j,Np,Sthresh,NMax,Brlist,Slowest) # modified




        # DEBUG print out freq
        frqdata = prtViblist(e2)

#       if os.path.isfile(freqfile):
#          stop("Error: this file already exists- "+freqfile)

#        f1h = open(freqfile,"w")
#        nfreq = len(e2)
#        f1h.write("   s ")
#       for l in range(nfreq):
#           f1h.write("vib"+str(l+1)+" ")
#        f1h.write("\n")

        f1h.write( str(Brlist[j].s)+" ")
        for f in frqdata:
            f1h.write( str(f)+" "   )
        f1h.write("\n")
    #print ("Output: Generalized vibrational frequencies written to "+freqfile)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




        # Calculate Curvature coupling coefficient
        # Input:
        #  1. normal mode
        #  2. Sq Mass
        #  3. Kappa
        #sys.stdout.write( "Bmunu "+str(Brlist[j].s)+" ")
        if tfcurvcpl == 1:
            pass
            #calcurvcpl(nm2,Brlist[j].Vmass,Brlist[j].Vkappa)


        # Calculate Coriolis coupling coefficient
        # Input:
        #   1. normal mode
        #   2. Mass
        #   3. Kappa
        #   4. Eta
        #   5. Hessian (un-mass-weighted)
        #   6. Cartesian coordinates
        #
        # Possible issue:
        #   1. mass-weighting of normal modes
        #   2. inconsistency of sign of eta for forward and reverse direction  (in CalPdP)
        #
        if tfcoriolis == 1:
            pass
            #calcoriolis( nm2,Brlist[j].Vmass,Brlist[j].Vkappa,Brlist[j].Veta,Brlist[j].Vhess,Brlist[j].Vcoor,Brlist[j].s  )





#
    print(("Output: Generalized vibrational frequencies written to "+freqfile))
#
        # DEBUG print out freq
#       frqdata = prtViblist(e2)

#       sys.stdout.write("freq "+ str(Brlist[j].s)+" ")
#       for f in frqdata:
#           sys.stdout.write( str(f)+" "   )
#       sys.stdout.write("\n")





def calcurvcpl(nm,Vmass,Vkappa):
    nvib = len(nm)
    natom = len(Vmass)
    natm3 = natom * 3

    SqMass = []
    Vsqmass = []

    for i in range(natom):
        for j in range(3):
            Vsqmass.append( Vmass[i] )

    for i in range(natom*3):
        SqMass.append( [] )
    for i in range(natom*3):
        for j in range(natom*3):
            SqMass[i].append( 0 )

    for k in range(3*natom):
        SqMass[k][k] = Vsqmass[k]

    Lkappa = [Vkappa]
    # Un-mass weight kappa
    for i in range(natm3):
        Lkappa[0][i] = Lkappa[0][i] / math.sqrt( SqMass[i][i] )

    # Un-mass weight normal mode
    nm_2 = []
    for i in range(nvib):
        nm_2.append( [] )
    for i in range(nvib):
        for j in range(natm3):
            #nm_2[i].append( 0 )
            tmp = nm[i][j]
            tmp = tmp / math.sqrt( Vsqmass[j] )
            nm_2[i].append( tmp )


    for i in range(nvib):
        length = np.linalg.norm( nm_2[i] )
        for j in range(natm3):
            nm_2[i][j] = nm_2[i][j] / length


    Amp1 = AnBInA( natm3,nvib,1,nm_2,SqMass,Lkappa,True )
    print(Amp1)
    pass







def solvLocDif(aNM,p,q,Np,Sthresh,NMax,Brlist,Slowest):
    # We also need to check overlap in this part

    # If the overlap check fails, increase Np, until it hits NMax

    Np_orig = Np
    NMax_orig = NMax
    aNM_orig = copy.deepcopy( aNM )

    # debug mode print option
    debug = True

    if (p+1) != q:
        #stop("Error: Two consecutive points are expected in solvLocDif...")
        print(("Warning: point(s) skipped between "+str(Brlist[p].s)+" and "+str(Brlist[q].s)))


    while( Np <= NMax ):
    #print (Np,NMax,Sthresh)

    # Three situations should be considered
    # 1. Run once, then done -> noerror = 1
    # 2. Run several times, then done -> noerror = 1
    # 3. Exceed NMax, still fail

        noerror = -1

        for i in range(Np):
            j = i + 1
            bE,bNM,bCAR=solWilsProj( \
                                     ABscIn(j,Np,Brlist[p].Vcoor,Brlist[q].Vcoor),\
                                     ABscIn(j,Np,Brlist[p].Vhess,Brlist[q].Vhess),\
                                     Brlist[p].Vmass,\
                                     ABscIn(j,Np,Brlist[p].Veta, Brlist[q].Veta) \
                                     # might need improvement for TS
                                     )

            nvib = len(aNM)
            natom3 = len(aNM[0])

            bNM,bE = DMODriv(aNM,bNM,bE,natom3,nvib)

            #aNM = bNM ---> replaced by deepcopy below


            # Check overlap
            tmpL = []
            for k in range(nvib):
                tmpL.append( np.dot(aNM[k],bNM[k]) )
            minv = min(tmpL)


            aNM = copy.deepcopy( bNM )


            if minv < Sthresh:
                #print ("We need to increase Np.")
                Np = Np + 20
                noerror = 0

                aNM = copy.deepcopy( aNM_orig ) # added later...

                if debug:
                    print(("break at "+str(i)+" Np:"+str(Np-20)+" minv:"+str(minv)))

                break # hopefully this could break the for loop

        # End of for loop ------------------------------------------

        if Np > NMax:
            if Sthresh >= Slowest:
                Np = Np_orig
                aNM = copy.deepcopy( aNM_orig )
                Sthresh = Sthresh - 0.02
                if debug:
                    print(("Warning: Np exceeds NMax, current overlap min: "+str(minv)+" Sthresh= "+str(Sthresh)  ))
            else:
                if debug:
                    print(("Error: Np exceeds NMax, current overlap min: "+str(minv)+" Sthresh= "+str(Sthresh)  ))
                stop("Error: Failed in resolve of local difficulty in DMO...")



        if noerror == -1: # Get through situation
            #if debug:
            print(("Np: "+str(Np)+" current overlap min: "+str(minv)))
            Np = NMax + 1 # In order to exit this while loop to avoid continuing the while loop
    # End of while loop -- -- -- -- --

    return (bE,bNM)








def overlap_dmo(l1,l2):
    nv1 = len(l1)
    nat1= len(l1[0])

    result=[]
    for i in range(nv1):
        p = np.dot( l1[i], l2[i] )
        result.append(p)
    #total = result #sum(result)
    minv = min(result)

#    # special debug
#    i = 5
#    j = 6
#    k = 7
#    Q= []
#    for i in [5,6,7]:
#           #print i
#           q = []
#           pass
#           for j in [5,6,7]:
#                p = np.dot( l1[i], l2[j] )
#               q.append(p)

#           Q.append(q)
#            print q
    return minv





## Check whether the IRC path is fully continuous or aborted somewhere in the middle
def chkcont(Brlist,io):
    numPn = len(Brlist)
    Thresh = 0.985

    loglist = []

    kpool = []
    for i in range(numPn):
        if Brlist[i].s == 0.0:
            kpool = [i-2, i-1, i, i+1, i+2]   # Modified for gold catalysis system



    if io == 1:
        stop("Error: io = 1 has not been tested in chkcont")

    if io == 0:
        for i in range(numPn-1):
            j = i + 1

            e1 = Brlist[i].Veta
            e2 = Brlist[j].Veta

            if  i in kpool:   #Brlist[i].s == 0.0:
                e1 = Brlist[i].Vnm_0
            if  j in kpool:   #Brlist[j].s == 0.0:
                e2 = Brlist[j].Vnm_0

            p = np.dot(e1,e2)
            #if p < 0:
            #   p = abs(p)
            #print Brlist[i].s, Brlist[j].s,p   # Uncomment this line for debug use

            if abs(p) < Thresh:
                row = [ i, j, p ]
                loglist.append( row )

    ndis = len(loglist)
    if ndis != 0:
        print ("Warning: Discontinuity found in the path...")
        print(("         Number of discontinuities: "+str(ndis)   ))
        for j in range(ndis):
            print(("         discontinuity found between "+str( Brlist[loglist[j][0]].s) +" "+ str( Brlist[loglist[j][1]].s) +" S = "+str(loglist[j][2])     ))
    else:
        print ("Info: No discontinuity found in the path.")

    return loglist




## Calculate path direction and curvature decomposition into internal coordinates
# caldircurcomp(oldnew, Brlist, numPn, parB, parID, "eta-q_n.csv", "kappa-q_n.csv" )
def caldircurcomp( oldnew, Brlist, numPn, parB ):

    nparm = len( parB[0] )

    # data structure
    #  [
    #  [s, q_n1, q_n2, ...],
    #  [                  ],
    #  ]
    #
    #

    etadata = []
    kappadata = []


    for i in range(numPn): # Loop over all points
        etarow = []
        kapparow = []
        etarow.append( Brlist[i].s )
        kapparow.append( Brlist[i].s )
        for j in range(nparm): # Loop over all parameters
            u1 = calU(i,j,Brlist,parB)
            if oldnew == 1:
                eta0 = [item for sublist in Brlist[i].Veta for item in sublist]
                kappa0 = [item for sublist in Brlist[i].Vkappa for item in sublist]
            elif oldnew == 0:
                eta0 = Brlist[i].Veta
                kappa0 = Brlist[i].Vkappa
            e1 = np.dot( u1, eta0 )
            if Brlist[i].s <= 0.0:
                e1 = e1 * -1
            etarow.append( e1 )
            e2 = np.dot( u1, kappa0 ) * -1  # multiply -1 # 05/17/2017
            kapparow.append( e2 )

        #print kapparow
        etadata.append( etarow )
        kappadata.append( kapparow )

    pass
    return (etadata,kappadata)



## Interface with alm to calculate the B matrix of ring coordinates
def calRingExB(na,coords,prm,basepath):
    deg2rad = 0.0174533
    ang2bohr = 1.889725989
    bohr2ang = 0.529177249
    fh1 = open("TMPLM0","w")
    fh1.write(str(na)+"\n") # XYZ file
    fh1.write("title here\n")
    for i in range(na):
        fh1.write("C " +" "+ str((float(coords[i][0]))*bohr2ang ) +" "\
                           + str((float(coords[i][1]))*bohr2ang ) +" "\
                           + str((float(coords[i][2]))*bohr2ang )  + "\n")  # Coordinates are in angstroms now
    fh1.write("\n")
    fh1.close()

    fh2 = open("TMPLM1","w")
    fh2.write("Title\n\n")
    fh2.write(" $contrl\n")
    fh2.write("  qcprog=\"xyz\" \n ")
    fh2.write("  IFMATLAB=.True.\n")
    fh2.write(" $end\n\n")
    fh2.write(" $qcdata\n")
    fh2.write("  fchk=\"TMPLM0\"   \n")
    fh2.write(" $end\n\n" )

    fh2.write(" $ring nring="+str(na)+" $end\n ")
    for i in range(na):
        fh2.write(str(i+1)+" "  )
    fh2.write("\n\n")

    fh2.write(" $LocMod  $End  \n")
    fh2.write(" 10  : noname\n")
    fh2.write("  "+str(prm[0])+" "+str(prm[1])+"\n\n"    )
    fh2.close()
    os.system(basepath+"lm90.test.exe -b < TMPLM1 > TMPLM1.log ")
    os.remove("TMPLM0")
    os.remove("TMPLM1")
    os.remove("TMPLM1.log")

    cmd1 = "grep \"b=\\[\" job\.m  -A 1 |grep -v b"
    tmp = os.popen(cmd1).read()
    #print tmp
    os.remove("job.m")
    tmp = tmp.split()
    value =[]
    for i in tmp:
        value.append( float(i) )
    #print value

    return value




## Interface with alm to calculate ring coordinates
def calRingEx(na,coords,prm,basepath):
    deg2rad = 0.0174533
    ang2bohr = 1.889725989
    bohr2ang = 1.0/ang2bohr
    fh1 = open("TMPLM0","w")
    fh1.write(str(na)+"\n") # XYZ file
    fh1.write("title here\n")
    for i in range(na):
        fh1.write("C " +" "+ str((float(coords[i][0]))*bohr2ang ) +" "\
                           + str((float(coords[i][1]))*bohr2ang ) +" "\
                           + str((float(coords[i][2]))*bohr2ang )  + "\n")  # Coordinates are in angstroms now
    fh1.write("\n")
    fh1.close()

    fh2 = open("TMPLM1","w")
    fh2.write("Title\n\n")
    fh2.write(" $contrl\n")
    fh2.write("  qcprog=\"xyz\" \n ")
    fh2.write(" $end\n\n")
    fh2.write(" $qcdata\n")
    fh2.write("  fchk=\"TMPLM0\"   \n")
    fh2.write(" $end\n\n" )

    fh2.write(" $ring nring="+str(na)+" $end\n ")
    for i in range(na):
        fh2.write(str(i+1)+" "  )
    fh2.write("\n\n")

    fh2.write(" $LocMod  $End  \n")
    fh2.write(" 10  : noname\n")
    fh2.write("  "+str(prm[0])+" "+str(prm[1])+"\n\n"    )
    fh2.close()
    os.system(basepath+"lm90.test.exe -b < TMPLM1 > TMPLM1.log ")
    os.remove("TMPLM0")
    os.remove("TMPLM1")

    cmd1=" grep \"q\_n        Name\" TMPLM1\.log  -A 2 |grep noname|awk -F' ' '{print $6}'     "
    tmp = os.popen(cmd1).read()
    os.remove("TMPLM1.log")
    value =  float(tmp)

    if prm[1] == 0 or prm[1] == 1 or prm[1] == 3:
        value = value * ang2bohr
    else:
        value = value * deg2rad

    #print value
    return value

## Interface with alm to calculate B matrix of out-of-plane angle
def calOopExB(c1,c2,c3,c4,basepath):
    deg2rad = 0.0174533
    bohr2ang = 0.529177249
    c1 = bohr2ang*c1
    c2 = bohr2ang*c2
    c3 = bohr2ang*c3
    c4 = bohr2ang*c4
    fh1 = open("TMPLM0","w")
    fh1.write(str(4)+"\n") # XYZ file first
    fh1.write("title here\n")
    fh1.write("C " +" "+ str(float(c1[0])) +" "+ str(float(c1[1])) +" "+ str(float(c1[2]))  + "\n" )
    fh1.write("C " +" "+ str(float(c2[0])) +" "+ str(float(c2[1])) +" "+ str(float(c2[2]))  + "\n" )
    fh1.write("C " +" "+ str(float(c3[0])) +" "+ str(float(c3[1])) +" "+ str(float(c3[2]))  + "\n" )
    fh1.write("C " +" "+ str(float(c4[0])) +" "+ str(float(c4[1])) +" "+ str(float(c4[2]))  + "\n" )
    fh1.write("\n")
    fh1.close()

    fh2 = open("TMPLM1","w")
    fh2.write("Title\n\n")
    fh2.write(" $contrl\n")
    fh2.write("  qcprog=\"xyz\" \n ")
    fh2.write("  IFMATLAB=.True. \n")
    fh2.write(" $end\n\n")
    fh2.write(" $qcdata\n")
    fh2.write("  fchk=\"TMPLM0\"   \n")
    fh2.write(" $end\n\n" )
    fh2.write(" $LocMod  $End  \n")
    fh2.write(" -1 2 3 4 \n\n")
    fh2.close()
    os.system(basepath+"lm90.test.exe -b < TMPLM1 > TMPLM1.log ")
    os.remove("TMPLM0")
    os.remove("TMPLM1")
    os.remove("TMPLM1.log")

    cmd1 = "grep \"b=\\[\" job\.m  -A 1 |grep -v b"
    tmp = os.popen(cmd1).read()
    #print tmp
    os.remove("job.m")
    tmp = tmp.split()
    value =[]
    for i in tmp:
        value.append( float(i) )

    return value




## Interface with alm to calculate out-of-plane angle
def calOopEx(c1,c2,c3,c4,basepath):
    deg2rad = 0.0174533
    fh1 = open("TMPLM0","w")
    fh1.write(str(4)+"\n") # XYZ file first
    fh1.write("title here\n")
    fh1.write("C " +" "+ str(float(c1[0])) +" "+ str(float(c1[1])) +" "+ str(float(c1[2]))  + "\n" )
    fh1.write("C " +" "+ str(float(c2[0])) +" "+ str(float(c2[1])) +" "+ str(float(c2[2]))  + "\n" )
    fh1.write("C " +" "+ str(float(c3[0])) +" "+ str(float(c3[1])) +" "+ str(float(c3[2]))  + "\n" )
    fh1.write("C " +" "+ str(float(c4[0])) +" "+ str(float(c4[1])) +" "+ str(float(c4[2]))  + "\n" )
    fh1.write("\n")
    fh1.close()

    fh2 = open("TMPLM1","w")
    fh2.write("Title\n\n")
    fh2.write(" $contrl\n")
    fh2.write("  qcprog=\"xyz\" \n ")
    fh2.write(" $end\n\n")
    fh2.write(" $qcdata\n")
    fh2.write("  fchk=\"TMPLM0\"   \n")
    fh2.write(" $end\n\n" )
    fh2.write(" $LocMod  $End  \n")
    fh2.write(" -1 2 3 4 \n\n")
    fh2.close()
    os.system(basepath+"lm90.test.exe -b < TMPLM1 > TMPLM1.log ")
    os.remove("TMPLM0")
    os.remove("TMPLM1")

    cmd1=" grep \"q\_n        Name\" TMPLM1\.log  -A 2 |grep ang|awk -F' ' '{print $6}'     "

    tmp = os.popen(cmd1).read()
    os.remove("TMPLM1.log")
    value =  float(tmp)*deg2rad
    return value

## Interface with alm to calcuate B matrix of pyramidalization angle
def calPyrExB(c1,c2,c3,c4,basepath):
    deg2rad = 0.0174533
    bohr2ang= 0.529177249

    fh1 = open("TMPLM0","w")
    fh1.write(str(4)+"\n")
    fh1.write("title here\n")
    c1 =bohr2ang*c1 #
    c2= bohr2ang*c2 #
    c3= bohr2ang*c3 #
    c4= bohr2ang*c4 #
    fh1.write("C " +" "+ str(float(c1[0])) +" "+ str(float(c1[1])) +" "+ str(float(c1[2]))  + "\n" )
    fh1.write("C " +" "+ str(float(c2[0])) +" "+ str(float(c2[1])) +" "+ str(float(c2[2]))  + "\n" )
    fh1.write("C " +" "+ str(float(c3[0])) +" "+ str(float(c3[1])) +" "+ str(float(c3[2]))  + "\n" )
    fh1.write("C " +" "+ str(float(c4[0])) +" "+ str(float(c4[1])) +" "+ str(float(c4[2]))  + "\n" )
    fh1.write("\n")
    fh1.close()

    fh2 = open("TMPLM1","w")
    fh2.write("Title\n\n")
    fh2.write(" $contrl\n")
    fh2.write("  qcprog=\"xyz\" \n ")
    fh2.write("  IFMATLAB= .True. \n ")  # that file is called "job.m"
    fh2.write(" $end\n\n")
    fh2.write(" $qcdata\n")
    fh2.write("  fchk=\"TMPLM0\"   \n")
    fh2.write(" $end\n\n" )
    fh2.write(" $LocMod  $End  \n")
    fh2.write(" 1 -2 -3 -4 \n\n")
    fh2.close()
    os.system(basepath+"lm90.test.exe -b < TMPLM1 > TMPLM1.log ")
    os.remove("TMPLM0")
    os.remove("TMPLM1")
    os.remove("TMPLM1.log")

    cmd1 = "grep \"b=\\[\" job\.m  -A 1 |grep -v b"
    tmp = os.popen(cmd1).read()
    #print tmp
    os.remove("job.m")
    tmp = tmp.split()
    value =[]
    for i in tmp:
        value.append( float(i) )

    return value


## Interface with alm to calculate pyramidalization angle
def calPyrEx(c1,c2,c3,c4,basepath):
    deg2rad = 0.0174533

    fh1 = open("TMPLM0","w")
    fh1.write(str(4)+"\n")
    fh1.write("title here\n")
    fh1.write("C " +" "+ str(float(c1[0])) +" "+ str(float(c1[1])) +" "+ str(float(c1[2]))  + "\n" )
    fh1.write("C " +" "+ str(float(c2[0])) +" "+ str(float(c2[1])) +" "+ str(float(c2[2]))  + "\n" )
    fh1.write("C " +" "+ str(float(c3[0])) +" "+ str(float(c3[1])) +" "+ str(float(c3[2]))  + "\n" )
    fh1.write("C " +" "+ str(float(c4[0])) +" "+ str(float(c4[1])) +" "+ str(float(c4[2]))  + "\n" )
    fh1.write("\n")
    fh1.close()

    fh2 = open("TMPLM1","w")
    fh2.write("Title\n\n")
    fh2.write(" $contrl\n")
    fh2.write("  qcprog=\"xyz\" \n ")
    fh2.write(" $end\n\n")
    fh2.write(" $qcdata\n")
    fh2.write("  fchk=\"TMPLM0\"   \n")
    fh2.write(" $end\n\n" )
    fh2.write(" $LocMod  $End  \n")
    fh2.write(" 1 -2 -3 -4 \n\n")
    fh2.close()
    os.system(basepath+"lm90.test.exe -b < TMPLM1 > TMPLM1.log ")
    os.remove("TMPLM0")
    os.remove("TMPLM1")

    cmd1=" grep \"q\_n        Name\" TMPLM1\.log  -A 2 |grep ang|awk -F' ' '{print $6}'     "

    tmp = os.popen(cmd1).read()
    os.remove("TMPLM1.log")
    value =  float(tmp)*deg2rad
    return  value



## calc. the value of parameter
### works for both new browsing file and old browsing file
def calparVal(numberofpoints,Brlist,parID,basepath):
    parVal=[]
    bar = ProgressBar(total = numberofpoints)
    for i in range(numberofpoints):
    #sys.stdout.write(".")
    #sys.stdout.flush()
        bar.move()
        bar.log('We have arrived at: ' + str(i + 1))
#
        parVal.append([])
        for  j in range(len(parID)):
            if parID[j][0]=="g": # ring coordinates
            # 0 -> g
            # 1 -> number of atoms
            # 2 -> [atom labels]
            # 3 -> [parameters]

                na = parID[j][1] # number of atoms
                coords=[]
                atmspec = parID[j][2]
                prm = parID[j][3]
                for n in range(na):
                    coor = np.array(Brlist[i].Vcoor[ atmspec[n] - 1 ])
                    coords.append(coor)
                #print len(Brlist[i].Vcoor)
                #print coords
                tmpV = calRingEx(na,coords,prm,basepath)
                tmp1=[]
                tmp1.append("g")
                tmp1.append(parID[j][-1])
                tmp1.append(tmpV)
                parVal[i].append( tmp1 )



            if parID[j][0]=="o": # out-of-plane angle
                r1=parID[j][1] - 1
                r2=parID[j][2] - 1
                r3=parID[j][3] - 1
                r4=parID[j][4] - 1
                c1=np.array(Brlist[i].Vcoor[r1])
                c2=np.array(Brlist[i].Vcoor[r2])
                c3=np.array(Brlist[i].Vcoor[r3])
                c4=np.array(Brlist[i].Vcoor[r4])
                tmpV = calOopEx(c1,c2,c3,c4,basepath)
                tmp1=[]
                tmp1.append("o")
                tmp1.append(parID[j][5])
                tmp1.append(tmpV)
                parVal[i].append( tmp1 )


            if parID[j][0]=="y": # pyramidalization angle
                r1=parID[j][1] - 1
                r2=parID[j][2] - 1
                r3=parID[j][3] - 1
                r4=parID[j][4] - 1
                c1=np.array(Brlist[i].Vcoor[r1])
                c2=np.array(Brlist[i].Vcoor[r2])
                c3=np.array(Brlist[i].Vcoor[r3])
                c4=np.array(Brlist[i].Vcoor[r4])
                tmpV = calPyrEx(c1,c2,c3,c4,basepath)
                tmp1=[]
                tmp1.append("y")
                tmp1.append(parID[j][5])
                tmp1.append(tmpV)
                parVal[i].append( tmp1  )

            if parID[j][0]=="r": # bond
                r1=parID[j][1] - 1
                r2=parID[j][2] - 1
                c1=np.array(Brlist[i].Vcoor[r1])
                c2=np.array(Brlist[i].Vcoor[r2])
                #print len(Brlist[i].Vcoor)
                bdr=np.linalg.norm(c1-c2)
                tmp1=[]
                tmp1.append("r")
                tmp1.append(parID[j][3])
                tmp1.append(bdr)
                parVal[i].append(tmp1)
                #print bdr
                #print c1
            if parID[j][0]=="a": # angle
                a1=parID[j][1] - 1
                a2=parID[j][2] - 1
                a3=parID[j][3] - 1
                c1=np.array(Brlist[i].Vcoor[a1])
                c2=np.array(Brlist[i].Vcoor[a2])
                c3=np.array(Brlist[i].Vcoor[a3])
                d12=np.linalg.norm(c1-c2)
                d32=np.linalg.norm(c3-c2)
                tmp2=np.inner(c1-c2,c3-c2)
                cphi=tmp2/(d12*d32)
                if abs(cphi) > 1:
                    cphi = 1.0 * np.sign(cphi)
                bda = math.acos(cphi)
                tmp1=[]
                tmp1.append("a")
                tmp1.append(parID[j][4])
                tmp1.append(bda)
                parVal[i].append(tmp1)
            if parID[j][0]=="d": # dihedral
                d1=parID[j][1] - 1
                d2=parID[j][2] - 1
                d3=parID[j][3] - 1
                d4=parID[j][4] - 1
                c1=np.array(Brlist[i].Vcoor[d1])
                c2=np.array(Brlist[i].Vcoor[d2])
                c3=np.array(Brlist[i].Vcoor[d3])
                c4=np.array(Brlist[i].Vcoor[d4])
                ba = c1 - c2
                bc = c3 - c2
                cd = c4 - c3
                u1 = np.cross(ba,bc)
                u2 = np.cross(cd,bc)
                u3 = np.cross(ba,cd)
                r1=np.linalg.norm(u1)
                r2=np.linalg.norm(u2)
                if ((r1<1.0e-5) or (r2<1.0e-5)):
                    tor=0
                r1 = 1/r1
                r2 = 1/r2
                u1 = u1*r1
                u2 = u2*r2
                tmp = np.inner(u1,u2)
                if abs(tmp) > 1:
                    tmp = 1.0 * np.sign(tmp)
                tor = math.acos(tmp)
                tmp=np.inner(u3,bc)
                tor = tor * np.sign(tmp)
                tmp1=[]
                tmp1.append("d")
                tmp1.append(parID[j][5])
                tmp1.append(tor)
                parVal[i].append(tmp1)
    #sys.stdout.write("\n")
    return parVal
###############################

### calculate dihedral
def torsion(c1,c2,c3,c4):
                #c1=np.array(Brlist[i].Vcoor[d1])
                #c2=np.array(Brlist[i].Vcoor[d2])
                #c3=np.array(Brlist[i].Vcoor[d3])
                #c4=np.array(Brlist[i].Vcoor[d4])
    ba = c1 - c2
    bc = c3 - c2
    cd = c4 - c3
    u1 = np.cross(ba,bc)
    u2 = np.cross(cd,bc)
    u3 = np.cross(ba,cd)
    r1=np.linalg.norm(u1)
    r2=np.linalg.norm(u2)
    if ((r1<1.0e-5) or (r2<1.0e-5)):
        tor=0
    r1 = 1/r1
    r2 = 1/r2
    u1 = u1*r1
    u2 = u2*r2
    tmp = np.inner(u1,u2)
    if abs(tmp) > 1:
        tmp = 1.0 * np.sign(tmp)
    tor = math.acos(tmp)

    return tor
###############

## nmGenBD
def nmGenBD(c1,c2,c3,c4):
    "Numerical method to calculate B matrix for dihedral"
    q0 = torsion(c1,c2,c3,c4)
    delta = 0.01
    bv=[]

    c11=copy.deepcopy(c1)
    c11[0] = c11[0] + delta
    bv.append( ( torsion(c11,c2,c3,c4)-q0 )/delta )
    c12=copy.deepcopy(c1)
    c12[1] = c12[1] + delta
    bv.append( ( torsion(c12,c2,c3,c4)-q0 )/delta )
    c13=copy.deepcopy(c1)
    c13[2] = c13[2] + delta
    bv.append( ( torsion(c13,c2,c3,c4)-q0 )/delta )

    c21=copy.deepcopy(c2)
    c21[0] = c21[0] + delta
    bv.append( ( torsion(c1,c21,c3,c4)-q0 )/delta )
    c22=copy.deepcopy(c2)
    c22[1] = c22[1] + delta
    bv.append( ( torsion(c1,c22,c3,c4)-q0 )/delta )
    c23=copy.deepcopy(c3)
    c23[2] = c23[2] + delta
    bv.append( ( torsion(c1,c23,c3,c4)-q0 )/delta )

    c31=copy.deepcopy(c3)
    c31[0] = c31[0] + delta
    bv.append( ( torsion(c1,c2,c31,c4)-q0 )/delta )
    c32=copy.deepcopy(c3)
    c32[1] = c32[1] + delta
    bv.append( ( torsion(c1,c2,c32,c4)-q0 )/delta )
    c33=copy.deepcopy(c3)
    c33[2] = c33[2] + delta
    bv.append( ( torsion(c1,c3,c33,c4)-q0 )/delta )

    c41=copy.deepcopy(c4)
    c41[0] = c41[0] + delta
    bv.append( ( torsion(c1,c2,c3,c41)-q0 )/delta )
    c42=copy.deepcopy(c4)
    c42[1] = c42[1] + delta
    bv.append( ( torsion(c1,c2,c3,c42)-q0 )/delta )
    c43=copy.deepcopy(c4)
    c43[2] = c43[2] + delta
    bv.append( ( torsion(c1,c2,c3,c43)-q0 )/delta )

    return bv


## Wenli's way to calculate b matrix for bond angles
### Give the same result as analytical way does
def GenBA(c1,c2,c3):
    "Calculate Wilson b matrix for bond angles in Wenli's way"
    bv=[]
    c1=np.array(c1)
    c2=np.array(c2)
    c3=np.array(c3)
    tmp1 = c1-c2
    tmp2 = c3-c2
    d12 = np.linalg.norm(tmp1)
    d32 = np.linalg.norm(tmp2)
    tmp1a = tmp1 * (1/d12)
    tmp2a = tmp2 * (1/d32)
    cphi = np.dot(tmp1a,tmp2a)
    sphi = (1-cphi*cphi)**0.5

    b1 = tmp2a * -1
    b1a= tmp1a * cphi + b1
    cf=1/(d12*sphi)
    b1b= b1a * cf

    b3 = tmp1a * -1
    b3a= tmp2a * cphi + b3
    cf=1/(d32*sphi)
    b3b= b3a * cf

    b2 = b1b * -1
    b2b= b2 - b3b

    for i in b1b:
        bv.append(i)
    for i in b2b:
        bv.append(i)
    for i in b3b:
        bv.append(i)
    return bv

#######################################


# B matrix based on Torsion displacement
def GenBD2(c1,c2,c3,c4):
    xa = c1[0]
    ya = c1[1]
    za = c1[2]
    xb = c2[0]
    yb = c2[1]
    zb = c2[2]
    xc = c3[0]
    yc = c3[1]
    zc = c3[2]
    xd = c4[0]
    yd = c4[1]
    zd = c4[2]


    bvAD = TorsDispB(  xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd  )
    bv = []

    for i in range(3):
        bv.append( bvAD[i] )
    for i in range(6):
        bv.append( 0.0)
    for i in range(3):
        j = i + 3
        bv.append( bvAD[j] )

    return bv


## Wenli's way to calculate b matrix for dihedrals
def GenBD(c1,c2,c3,c4):
    "Calculate wilson b matrix for dihedral angle in wenli's way"
    bv=[]
    c1=np.array(c1)
    c2=np.array(c2)
    c3=np.array(c3)
    c4=np.array(c4)

    tmp1 = c2 - c1
    tmp2 = c3 - c2
    tmp3 = c4 - c3
    d21 = np.linalg.norm(tmp1)
    d32 = np.linalg.norm(tmp2)
    d43 = np.linalg.norm(tmp3)

    tmp1a= tmp1*(1/d21)
    tmp2a= tmp2*(1/d32)
    tmp3a= tmp3*(1/d43)

    cph1 = -1 * np.dot(tmp1a, tmp2a)
    sph1 = (1 - cph1 * cph1)**0.5
    cph2 = -1 * np.dot(tmp3a, tmp2a)
    sph2 = (1 - cph2 * cph2)**0.5

    tmp4 = np.cross(tmp1a,tmp2a)
    tmp5 = np.cross(tmp2a,tmp3a)

    cf= -1/(d21*sph1*sph1)
    b1= tmp4*cf  ###
    cf=  1/(d43*sph2*sph2)
    b4= tmp5*cf  ####
    cf= (d32-d21*cph1)/(d32*d21*sph1*sph1)
    b2= tmp4*cf
    cf = cph2/(d32*sph2*sph2)
    b3 = tmp5*cf
    b2a=b2 - b3   #####
    cf= -1
    b3= b1*cf
    b3a=b3 - b2
    b3b=b3a- b4    #####

    for i in b1:
        bv.append(i)
    for i in b2a:
        bv.append(i)
    for i in b3b:
        bv.append(i)
    for i in b4:
        bv.append(i)


    return bv


##################################################


## B matrix calculation engine
### work for both new browsing file and old browsing file
def calparB(numberofpoints,Brlist,parID,basepath):
    "Calculate the B matrix for parameters"
    x1,x2,x3,x4 = symbols('x1 x2 x3 x4')
    y1,y2,y3,y4 = symbols('y1 y2 y3 y4')
    z1,z2,z3,z4 = symbols('z1 z2 z3 z4')

    bondr = ((x1-x2)**2 +  \
             (y1-y2)**2 +  \
             (z1-z2)**2)**0.5
    #drdx1 = diff(bondr,x1)
    bonda = acos(  \
                 ((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3))/ \
                 (sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* \
                  sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2))     )
    #print diff(bonda,x1)
    parB=[]

    bar = ProgressBar(total = numberofpoints)

    for i in range(numberofpoints):
        #sys.stdout.write(".")
        #sys.stdout.flush()
        bar.move()
        bar.log('We have arrived at: ' + str(i + 1))
        parB.append([])
        for j in range(len(parID)):
            if parID[j][0]=="g": # ring coordinates
            # 0 -> g
            # 1 -> number of atoms
            # 2 -> [atom labels]
            # 3 -> [parameters]
                na = parID[j][1]
                coords=[]
                atmspec = parID[j][2]
                prm = parID[j][3]
                for n in range(na):
                    coor = np.array(Brlist[i].Vcoor[ atmspec[n] - 1 ])
                    coords.append(coor)
                bv = calRingExB(na,coords,prm,basepath) # basepath added
                tmp1=[]
                tmp1.append("g")
                tmp1.append( atmspec )
                tmp1.append( bv )
                #print tmp1
                parB[i].append(tmp1)   # "How to expand the compact B matrix into full length B matrix?"  implemented or not?


            if parID[j][0]=="o": # out-of-plane angle
                r1=parID[j][1] - 1
                r2=parID[j][2] - 1
                r3=parID[j][3] - 1
                r4=parID[j][4] - 1
                c1=np.array(Brlist[i].Vcoor[r1])
                c2=np.array(Brlist[i].Vcoor[r2])
                c3=np.array(Brlist[i].Vcoor[r3])
                c4=np.array(Brlist[i].Vcoor[r4])
                bv = calOopExB( c1,c2,c3,c4 ,basepath)
                tmp1=[]
                tmp1.append("o")
                tmp1.append(r1+1)
                tmp1.append(r2+1)
                tmp1.append(r3+1)
                tmp1.append(r4+1)
                for k in bv:
                    tmp1.append(k)
                parB[i].append(tmp1)
                #print tmp1

            if parID[j][0]=="y": # pyramidalization angle
                r1=parID[j][1] - 1
                r2=parID[j][2] - 1
                r3=parID[j][3] - 1
                r4=parID[j][4] - 1
                c1=np.array(Brlist[i].Vcoor[r1])
                c2=np.array(Brlist[i].Vcoor[r2])
                c3=np.array(Brlist[i].Vcoor[r3])
                c4=np.array(Brlist[i].Vcoor[r4])
                bv = calPyrExB( c1,c2,c3,c4,basepath)
                tmp1=[]
                tmp1.append("y")
                tmp1.append(r1+1)
                tmp1.append(r2+1)
                tmp1.append(r3+1)
                tmp1.append(r4+1)
                for k in bv:
                    tmp1.append(k)
                parB[i].append(tmp1)

            if parID[j][0]=="r": # bond
                r1=parID[j][1] - 1
                r2=parID[j][2] - 1
                c1=np.array(Brlist[i].Vcoor[r1])
                c2=np.array(Brlist[i].Vcoor[r2])
                tmp=c2-c1
                d=1/np.linalg.norm(tmp)
                b31=-1*d*tmp
                b32=   d*tmp
                tmp1=[]
                tmp1.append("r")
                tmp1.append(r1+1)
                tmp1.append(r2+1)
                for k in b31:
                    tmp1.append(k)
                for k in b32:
                    tmp1.append(k)
                parB[i].append(tmp1)
                #For validation
                #b11=drdx1.subs([  (x1,c1[0]),(y1,c1[1]),(z1,c1[2]),(x2,c2[0]),(y2,c2[1]),(z2,c2[2])       ])
                #print b11
                #print b31[0]
                #print ""
            if parID[j][0]=="a": # angle
                a1=parID[j][1] - 1
                a2=parID[j][2] - 1
                a3=parID[j][3] - 1
                c1=np.array(Brlist[i].Vcoor[a1])
                c2=np.array(Brlist[i].Vcoor[a2])
                c3=np.array(Brlist[i].Vcoor[a3])
                tmp1=[]
                tmp1.append("a")
                tmp1.append(a1+1)
                tmp1.append(a2+1)
                tmp1.append(a3+1)
                b11=(diff(bonda,x1)).subs([ (x1,c1[0]),(y1,c1[1]),(z1,c1[2]),\
                                            (x2,c2[0]),(y2,c2[1]),(z2,c2[2]),\
                                            (x3,c3[0]),(y3,c3[1]),(z3,c3[2])  ])
                b21=(diff(bonda,y1)).subs([ (x1,c1[0]),(y1,c1[1]),(z1,c1[2]),\
                                            (x2,c2[0]),(y2,c2[1]),(z2,c2[2]),\
                                            (x3,c3[0]),(y3,c3[1]),(z3,c3[2])  ])
                b31=(diff(bonda,z1)).subs([ (x1,c1[0]),(y1,c1[1]),(z1,c1[2]),\
                                            (x2,c2[0]),(y2,c2[1]),(z2,c2[2]),\
                                            (x3,c3[0]),(y3,c3[1]),(z3,c3[2])  ])

                b12=(diff(bonda,x2)).subs([ (x1,c1[0]),(y1,c1[1]),(z1,c1[2]),\
                                            (x2,c2[0]),(y2,c2[1]),(z2,c2[2]),\
                                            (x3,c3[0]),(y3,c3[1]),(z3,c3[2])  ])
                b22=(diff(bonda,y2)).subs([ (x1,c1[0]),(y1,c1[1]),(z1,c1[2]),\
                                            (x2,c2[0]),(y2,c2[1]),(z2,c2[2]),\
                                            (x3,c3[0]),(y3,c3[1]),(z3,c3[2])  ])
                b32=(diff(bonda,z2)).subs([ (x1,c1[0]),(y1,c1[1]),(z1,c1[2]),\
                                            (x2,c2[0]),(y2,c2[1]),(z2,c2[2]),\
                                            (x3,c3[0]),(y3,c3[1]),(z3,c3[2])  ])

                b13=(diff(bonda,x3)).subs([ (x1,c1[0]),(y1,c1[1]),(z1,c1[2]),\
                                            (x2,c2[0]),(y2,c2[1]),(z2,c2[2]),\
                                            (x3,c3[0]),(y3,c3[1]),(z3,c3[2])  ])
                b23=(diff(bonda,y3)).subs([ (x1,c1[0]),(y1,c1[1]),(z1,c1[2]),\
                                            (x2,c2[0]),(y2,c2[1]),(z2,c2[2]),\
                                            (x3,c3[0]),(y3,c3[1]),(z3,c3[2])  ])
                b33=(diff(bonda,z3)).subs([ (x1,c1[0]),(y1,c1[1]),(z1,c1[2]),\
                                            (x2,c2[0]),(y2,c2[1]),(z2,c2[2]),\
                                            (x3,c3[0]),(y3,c3[1]),(z3,c3[2])  ])
                tmp1.append(b11)
                tmp1.append(b21)
                tmp1.append(b31)
                tmp1.append(b12)
                tmp1.append(b22)
                tmp1.append(b32)
                tmp1.append(b13)
                tmp1.append(b23)
                tmp1.append(b33)
                parB[i].append(tmp1)
                #print tmp1
                # For validation use
                #tmpBA=[]
                #tmpBA=GenBA(c1,c2,c3)
                #print tmpBA
                #print [b11,b21,b31,b12,b22,b32,b13,b23,b33]
                #print ""
            if parID[j][0]=="d": # dihedral:
                d1=parID[j][1] - 1
                d2=parID[j][2] - 1
                d3=parID[j][3] - 1
                d4=parID[j][4] - 1
                c1=np.array(Brlist[i].Vcoor[d1])
                c2=np.array(Brlist[i].Vcoor[d2])
                c3=np.array(Brlist[i].Vcoor[d3])
                c4=np.array(Brlist[i].Vcoor[d4])
                tmp1=[]
                tmp1.append("d")
                tmp1.append( d1 + 1 )
                tmp1.append( d2 + 1 )
                tmp1.append( d3 + 1 )
                tmp1.append( d4 + 1 )
                #tmpBD=[]

                tmpBD=GenBD(c1,c2,c3,c4) # calculate B matrix # old way

                #tmpBD=GenBD2(c1,c2,c3,c4) # Try new way

                #print tmpBD
                #tmpBD = nmGenBD(c1,c2,c3,c4)
                #print tmpBD
                for m in tmpBD:
                    tmp1.append(m)
                parB[i].append(tmp1)



    #sys.stdout.write("\n")
    #print ("B matrix calculation completed.")
    return parB


#######################################


## Calculate the u-vector of URVA
def calU(pID,bID,Brlist,parB):
    "calculate the u-vector \
     for a specific parameter(bID) at a specific \
     point(pID) on reaction path for componenet analysis"

    G = calG(pID,bID,Brlist,parB)
    Vmass = copy.deepcopy(Brlist[pID].Vmass)
    Natom = Brlist[pID].NAtom
    b = makB(Natom,pID,bID,parB)
    #print b
    tmp1=[]
    for i in range(len(b)):
        j = i / 3
        tmp1.append( b[i] * Vmass[j]**(-0.5) )   #G**(-0.5) )

    u = np.array(tmp1)
    unm = np.linalg.norm( u )
    u = u / unm
    return u


###################################




## Calculate Wilson G matrix element
def calG(pID,bID,Brlist,parB):
    "calculate the wilson G matrix element \
     for a specific parameter(bID) at a specfic \
     point(pID) on reaction path"
    #print pID
    Vmass = Brlist[pID].Vmass
    Natom = Brlist[pID].NAtom
    b = makB(Natom,pID,bID,parB)
    tmp1=[]
    for i in range(len(b)):
        j = i / 3
        tmp1.append( b[i] / Vmass[j] )
    tmp1 = np.array(tmp1)
    b = np.array(b)
    G = np.dot( tmp1, b)

    return G
    #print b[0:20]

###########################################


## Recover the full-length B vector
def makB(Natom,pID,bID,parB):
    "recover the full b vector"
    b=np.zeros(3*Natom)

    # many parameters waiting to be implemented here...
        # 1. pyramidalization angle
        # 2. out-of-plane angle
        # 3. ring coordinates (done)

    if parB[pID][bID][0] == "g": # ring coordinate
        atmspec = parB[pID][bID][1]
        bv =      parB[pID][bID][2]
        nring = len(atmspec) # number of atoms in this ring
        for j in range(len(atmspec)):
            r = atmspec[j] - 1
            for k in range(3):
                index1 = 3*r + k
                index2 = 3*j + k
                b[ index1 ] = bv[ index2 ]


    if parB[pID][bID][0] == "r": # bond
        r1 = parB[pID][bID][1] - 1
        r2 = parB[pID][bID][2] - 1
        b[ 3*r1+0 ] = parB[pID][bID][3]
        b[ 3*r1+1 ] = parB[pID][bID][4]
        b[ 3*r1+2 ] = parB[pID][bID][5]
        b[ 3*r2+0 ] = parB[pID][bID][6]
        b[ 3*r2+1 ] = parB[pID][bID][7]
        b[ 3*r2+2 ] = parB[pID][bID][8]

    if parB[pID][bID][0] == "a": # angle
        r1 = parB[pID][bID][1] - 1
        r2 = parB[pID][bID][2] - 1
        r3 = parB[pID][bID][3] - 1
        b[ 3*r1+0 ] = parB[pID][bID][4]
        b[ 3*r1+1 ] = parB[pID][bID][5]
        b[ 3*r1+2 ] = parB[pID][bID][6]
        b[ 3*r2+0 ] = parB[pID][bID][7]
        b[ 3*r2+1 ] = parB[pID][bID][8]
        b[ 3*r2+2 ] = parB[pID][bID][9]
        b[ 3*r3+0 ] = parB[pID][bID][10]
        b[ 3*r3+1 ] = parB[pID][bID][11]
        b[ 3*r3+2 ] = parB[pID][bID][12]

    if parB[pID][bID][0] == "d": # diheral
        r1 = parB[pID][bID][1] - 1
        r2 = parB[pID][bID][2] - 1
        r3 = parB[pID][bID][3] - 1
        r4 = parB[pID][bID][4] - 1
        b[ 3*r1+0 ] = parB[pID][bID][5]
        b[ 3*r1+1 ] = parB[pID][bID][6]
        b[ 3*r1+2 ] = parB[pID][bID][7]
        b[ 3*r2+0 ] = parB[pID][bID][8]
        b[ 3*r2+1 ] = parB[pID][bID][9]
        b[ 3*r2+2 ] = parB[pID][bID][10]
        b[ 3*r3+0 ] = parB[pID][bID][11]
        b[ 3*r3+1 ] = parB[pID][bID][12]
        b[ 3*r3+2 ] = parB[pID][bID][13]
        b[ 3*r4+0 ] = parB[pID][bID][14]
        b[ 3*r4+1 ] = parB[pID][bID][15]
        b[ 3*r4+2 ] = parB[pID][bID][16]




    return b



####################################





## calculate path direction and curvature for OLD browsing file
def dircur(Brlist):
    numP = len(Brlist)
    for i in range(numP):
        #print Brlist[i].Vgrad
        #force = np.linalg.norm(Brlist[i].Vgrad)
        #print (force)
        NAt3 = 3* Brlist[0].NAtom
        eta = []
        for j in range( NAt3 ):
            k = j / 3
            tmp1 = Brlist[i].Vgrad[j] *  (Brlist[i].Vmass[k])**(-0.5)
            eta.append (tmp1)
        eta = np.array(eta)
        mwforce = np.linalg.norm(eta)
        eta = eta / mwforce  # PATH DIRECTION
        #print (np.linalg.norm(eta))

        mwhess = LT2Sqr( NAt3 , Brlist[i].Vhess )
        for j in range(NAt3):
            for k in range(NAt3):
                jm = j / 3
                km = k / 3
                mwhess[j][k] = mwhess[j][k] /( (Brlist[i].Vmass[jm])**0.5 *  (Brlist[i].Vmass[km])**0.5 )

        mwhess = np.matrix( mwhess )
        left =   (mwhess.dot(eta)).tolist()[0]
        #print (left)
        left = np.array(left)
        efe = np.dot(eta,left)
        right = eta * efe
        #print (right)
        cur = (left - right) / mwforce   # CURVATURE
        #print Brlist[i].s
        #print np.linalg.norm(cur)
        Brlist[i].Veta  = copy.deepcopy ( eta )
        Brlist[i].Vkappa= copy.deepcopy ( cur )


#############################

## Solve Wilson Equation with projection
def solWilsProj(Vcoor,Vhess,Vmass,Proj):
    NAtom = len(Vmass)
    NAt3 = 3*NAtom
    mwhess = LT2Sqr( 3*NAtom, Vhess)
    for j in range(NAt3):
        for k in range(NAt3):
            jm = j / 3
            km = k / 3
            mwhess[j][k] = mwhess[j][k] /\
            ( (Vmass[jm])**0.5 *  (Vmass[km])**0.5 )

    M,VecOut = GenTRQProj(Vcoor,Vmass,Proj)
    DMat = VecOut[M:3*NAtom] # Take part of VecOut
    #print (len(DMat))
    df1 = (MatMp1(DMat,mwhess)).transpose()
    #print (df1)
    fint = MatMp1(df1,DMat)
    del df1
    gc.collect()
    #print (fint)
    Eig,Ev = np.linalg.eig(fint)
    Ev = Ev.transpose()

    #Correct Sign of eigenvectors
    Ev = CorSignVec(Ev)



    #for i in Eig:
    #   if i > 0:
    #      print ((i**0.5)*5140.48715246,"cm^-1")
    #print Ev
    NModes = np.array(MatMp1(np.array(DMat).transpose(),Ev))
    #print NModes

    for i in range(len(NModes)):
        v = copy.deepcopy( NModes[i] )
        vnorm = np.linalg.norm(v)
        NModes[i] = v / vnorm
        #print np.linalg.norm(NModes[i])

    CartDisp=[]
    for i in range(len(NModes)):
        v = copy.deepcopy( NModes[i] )
        for j in range(NAt3):
            jm = j / 3
            v[j] = v[j] / (Vmass[jm]**0.5)
        vnorm = np.linalg.norm(v)
        CartDisp.append( (v/vnorm).tolist() )
    #print np.linalg.norm(CartDisp[1])

    EigSort =  np.sort(Eig)
    NModesSort = NModes[Eig.argsort(),:]
    CartDispSort = np.array(CartDisp)[Eig.argsort(),:]
    #print EigSort
    #print NModesSort
    #print CartDispSort
    #print ""
    #print Eig
    #print NModes
    #print CartDisp
    #Eig = Eig[idx]
    #NModes = NModes[idx]
    #CartDisp = CartDisp[idx]

    return (EigSort,NModesSort,CartDispSort)
#


###########################
## Solve Wilson Equation
def solWils(Vcoor,Vhess,Vmass):
    NAtom = len(Vmass)
    NAt3 = 3*NAtom
    mwhess = LT2Sqr( 3*NAtom, Vhess)
    ek,vk = np.linalg.eig(mwhess)
    for j in range(NAt3):
        for k in range(NAt3):
            jm = j / 3
            km = k / 3
            mwhess[j][k] = mwhess[j][k] /\
            ( (Vmass[jm])**0.5 *  (Vmass[km])**0.5 )

    M,VecOut = GenTRQ(Vcoor,Vmass) # In most cases, M = 6
    DMat = VecOut[M:3*NAtom] # Take part of VecOut
    #print (len(DMat))
    df1 = (MatMp1(DMat,mwhess)).transpose()
    #print (df1)
    fint = MatMp1(df1,DMat)
    del df1
    gc.collect()
    #print (fint)
    Eig,Ev = np.linalg.eig(fint)
    Ev = Ev.transpose()

    #Correct Sign of eigenvectors
    Ev = CorSignVec(Ev)

    #for i in Eig:
    #   if i > 0:
    #      print ((i**0.5)*5140.48715246,"cm^-1")
    #print Ev
    NModes = np.array(MatMp1(np.array(DMat).transpose(),Ev))
    #print NModes

    for i in range(len(NModes)):
        v = copy.deepcopy( NModes[i] )
        vnorm = np.linalg.norm(v)
        NModes[i] = v / vnorm
        #print np.linalg.norm(NModes[i])

    CartDisp=[]
    for i in range(len(NModes)):
        v = copy.deepcopy( NModes[i] )
        for j in range(NAt3):
            jm = j / 3
            v[j] = v[j] / (Vmass[jm]**0.5)
        vnorm = np.linalg.norm(v)
        CartDisp.append( (v/vnorm).tolist() )
    #print np.linalg.norm(CartDisp[1])

    EigSort =  np.sort(Eig)
    NModesSort = NModes[Eig.argsort(),:]
    CartDispSort = np.array(CartDisp)[Eig.argsort(),:]
    #print EigSort
    #print NModesSort
    #print CartDispSort
    #print ""
    #print Eig
    #print NModes
    #print CartDisp
    #Eig = Eig[idx]
    #NModes = NModes[idx]
    #CartDisp = CartDisp[idx]

    return (EigSort,NModesSort,CartDispSort)
########################
def GenTRQProj(CC,AtMass,Proj):
    NAtom = len(AtMass)
    VecOut=[]
    for i in range(3*NAtom):
        VecOut.append([])
    for i in range(3*NAtom):
        for j in range(3*NAtom):
            VecOut[i].append(0)

    CMCoor,Rot = MOfI(CC,AtMass)
    M,VecOut = TRVect(AtMass,CC,CMCoor,Rot,VecOut)
    #print M
    #VecOut[M] = Proj
    #print VecOut
    for i in range(len(Proj)):
        VecOut[M][i] = Proj[i]
    M = M + 1
    #print VecOut
    VecOut = GSorth(0,3*NAtom,M,VecOut)  ###
    #print ("M=",M)
    #print VecOut
    return (M,VecOut)




########################
def GenTRQ(CC,AtMass):
    NAtom = len(AtMass)
    VecOut=[]
    for i in range(3*NAtom):
        VecOut.append([])
    for i in range(3*NAtom):
        for j in range(3*NAtom):
            VecOut[i].append(0)

    CMCoor,Rot = MOfI(CC,AtMass)
    M,VecOut = TRVect(AtMass,CC,CMCoor,Rot,VecOut)

    VecOut = GSorth(0,3*NAtom,M,VecOut)  ###
    #print ("M=",M)
    return (M,VecOut)


########################
def MOfI(C,AtMass):
    NAtom = len(AtMass)
    t1 = 0
    t2 = 0
    t3 = 0
    m = 0
    for i in range(NAtom):
        m = m + AtMass[i]
        t1 = t1 + C[i][0] * AtMass[i]
        t2 = t2 + C[i][1] * AtMass[i]
        t3 = t3 + C[i][2] * AtMass[i]
    COM=[t1/m,t2/m,t3/m]
    T=[]
    for i in range(6):
        T.append(0)

    for i in range(NAtom):
        Wt = AtMass[i]
        X = C[i][0] - COM[0]
        Y = C[i][1] - COM[1]
        Z = C[i][2] - COM[2]
        T[1-1] = T[1-1] + Wt * (Y*Y+Z*Z)
        T[3-1] = T[3-1] + Wt * (X*X+Z*Z)
        T[6-1] = T[6-1] + Wt * (X*X+Y*Y)
        T[2-1] = T[2-1] - Wt * X * Y
        T[4-1] = T[4-1] - Wt * X * Z
        T[5-1] = T[5-1] - Wt * Y * Z

    Tm=[]
    Tm.append([T[0],T[1],T[3]])
    Tm.append([T[1],T[2],T[4]])
    Tm.append([T[3],T[4],T[5]])

    #print Tm
    PMom,EigVec = np.linalg.eig(Tm)
    #print PMom,EigVec

    return (COM,EigVec)
##################################################


## TRVect
def TRVect(AtMass,C,CMCoor,RotEig,Vect):
    Cutoff = 1.0e-12
    NAtom = len(AtMass)
    NAt3 = 3*NAtom

    for i in range(NAtom):
        SqrtMI = (AtMass[i])**0.5
        CX = C[i][0] - CMCoor[0]
        CY = C[i][1] - CMCoor[1]
        CZ = C[i][2] - CMCoor[2]
        CXP = CX*RotEig[0][0] + CY*RotEig[0][1] + CZ*RotEig[0][2]
        CYP = CX*RotEig[1][0] + CY*RotEig[1][1] + CZ*RotEig[1][2]
        CZP = CX*RotEig[2][0] + CY*RotEig[2][1] + CZ*RotEig[2][2]
        Vect[0][3*i+0] = SqrtMI
        Vect[1][3*i+1] = SqrtMI
        Vect[2][3*i+2] = SqrtMI
        Vect[3][3*i+0] = (CYP*RotEig[2][0] - CZP*RotEig[1][0]) *SqrtMI
        Vect[3][3*i+1] = (CYP*RotEig[2][1] - CZP*RotEig[1][1]) *SqrtMI
        Vect[3][3*i+2] = (CYP*RotEig[2][2] - CZP*RotEig[1][2]) *SqrtMI
        Vect[4][3*i+0] = (CZP*RotEig[0][0] - CXP*RotEig[2][0]) *SqrtMI
        Vect[4][3*i+1] = (CZP*RotEig[0][1] - CXP*RotEig[2][1]) *SqrtMI
        Vect[4][3*i+2] = (CZP*RotEig[0][2] - CXP*RotEig[2][2]) *SqrtMI
        Vect[5][3*i+0] = (CXP*RotEig[1][0] - CYP*RotEig[0][0]) *SqrtMI
        Vect[5][3*i+1] = (CXP*RotEig[1][1] - CYP*RotEig[0][1]) *SqrtMI
        Vect[5][3*i+2] = (CXP*RotEig[1][2] - CYP*RotEig[0][2]) *SqrtMI

    NTrRo = 6
    I = 1
    while I <= NTrRo:
        X = np.dot(Vect[I-1],Vect[I-1])
        if X < Cutoff:
            NTrRo = NTrRo - 1
            Vect[I-1] = Vect[NTrTo]
        else:
            X = 1 / (X**0.5)
            for j in range(3*NAtom):
                Vect[I-1][j] = Vect[I-1][j] * X
            I = I + 1

    return (NTrRo,Vect)




## Correct sign of eigen vectors
def CorSignVec(a):
    la = len(a)
    for i in range(la):
        flag = 1.0
        v = (a[i]).tolist()
        vmin = min(v)
        vmax = max(v)
        if abs(vmin) > abs(vmax):
            flag = -1.0
            a[i] = flag * a[i]

    return a





# CalCdC for calculating Coriolis coupling coefficients
def CalCdC(natom,AtMass,CC,Eta):

    # FXMw should be a list FXMw[natom][3]
    # FXMw is actually Eta according to the comments in old version
    # So we need to convert Eta into FXMw with regard to the format

    FXMw = []
    for i in range(natom):
        FXMw.append( [] )
    for i in range(natom):
        for j in range(3):
            FXMw[i].append( 0 )


    for i in range(natom):
        for j in range(3):
            FXMw[i][j] = Eta[ 3*i + j ]


    Tresh = 1e-10
    nat3 = 3*natom

    # MOfI
    CMCoor,V = MOfI(CC,AtMass)

    # Initialization of C0 and dC0
    C0 = []
    dC0 = []
    for i in range(6):
        C0.append( [] )
        dC0.append( [] )
    for i in range(6):
        for j in range(natom):
            C0[i].append([])
            dC0[i].append([])
    for i in range(6):
        for j in range(natom):
            for k in range(3):
                C0[i][j].append( 0 )
                dC0[i][j].append( 0 )
    #
    for IAt in range(natom):
        SqrtM = math.sqrt( AtMass[IAt] )
        X = ( CC[IAt][0] - CMCoor[0] ) * SqrtM
        Y = ( CC[IAt][1] - CMCoor[1] ) * SqrtM
        Z = ( CC[IAt][2] - CMCoor[2] ) * SqrtM
        C0[0][IAt][0] = SqrtM
        C0[1][IAt][1] = SqrtM
        C0[2][IAt][2] = SqrtM

        # The following part might involve some bugs, but I wrote it with
        # the help of other parts
        C0[3][IAt][0] = -Z * V[1][0] + Y * V[2][0]
        C0[3][IAt][1] = -X * V[2][0] + Z * V[0][0]
        C0[3][IAt][2] = -Y * V[0][0] + X * V[1][0]

        C0[4][IAt][0] = -Z * V[1][1] + Y * V[2][1]
        C0[4][IAt][1] = -X * V[2][1] + Z * V[0][1]
        C0[4][IAt][2] = -Y * V[0][1] + X * V[1][1]

        C0[5][IAt][0] =  Z * V[1][2] - Y * V[2][2]
        C0[5][IAt][1] =  X * V[2][2] - Z * V[0][2] # fixed
        C0[5][IAt][2] =  Y * V[0][2] - X * V[1][2]

        #
        X = FXMw[IAt][0]
        Y = FXMw[IAt][1]
        Z = FXMw[IAt][2]

        #
        dC0[3][IAt][0] = -Z * V[1][0] + Y * V[2][0]
        dC0[3][IAt][1] = -X * V[2][0] + Z * V[0][0]
        dC0[3][IAt][2] = -Y * V[0][0] + X * V[1][0]

        dC0[4][IAt][0] = -Z * V[1][1] + Y * V[2][1]
        dC0[4][IAt][1] = -X * V[2][1] + Z * V[0][1]
        dC0[4][IAt][2] = -Y * V[0][1] + X * V[1][1]

        dC0[5][IAt][0] =  Z * V[1][2] - Y * V[2][2]
        dC0[5][IAt][1] =  X * V[2][2] - Z * V[0][2]
        dC0[5][IAt][2] =  Y * V[0][2] - X * V[1][2]


    #
    NTrRo = 0

    for I in range(6):
        X = math.sqrt( SProd2lay(C0[I],C0[I]) ) # Might lead to bug
        if X > Tresh:
            NTrRo = NTrRo + 1
            C0[NTrRo - 1] = copy.deepcopy(C0[I])  # Might lead to bug
            dC0[NTrRo - 1] = copy.deepcopy(dC0[I]) # Might lead to bug

    return (NTrRo,C0,dC0)



# CalQRS for calculating coriolis coupling coefficient
def CalQRS(natom, NTrRo,C0,dC0):
    nat3 = natom * 3
    Tresh = 1.0e-6

    # Initialize S,S1
    S=[]
    S1=[]
    for i in range(NTrRo):
        S.append( [] )
        S1.append( [] )
    for i in range(NTrRo):
        for j in range(NTrRo):
            S[i].append( 0 )
            S1[i].append( 0 )


    # Calculate S
    for i in range(NTrRo):
        for j in range(NTrRo):
            S[i][j] = SProd2lay(C0[i],C0[j])  # Might lead to bug


    V = copy.deepcopy( S )
    Eig,Ev = np.linalg.eig( V )
    Ev = Ev.transpose()
    for I in range(NTrRo):
        for J in range(NTrRo):
            S1IJ = 0.0
            for K in range(NTrRo):
                X = Eig[ K ]
                if X <= Tresh:
                    Y = 0.0
                else:
                    Y = 1.0 / X
                S1IJ = S1IJ + Ev[K][I] * Ev[K][J] * Y
            S1[J][I] = S1IJ

    # Initialize Tmp1
    Tmp1 = []
    for i in range(NTrRo):
        Tmp1.append( [] )
    for i in range(NTrRo):
        for j in range(nat3):
            Tmp1[i].append( 0 )



    # calculate Tmp1
    for loop1 in range(NTrRo):
        for loop2 in range(natom):
            for loop3 in range(3):
                #Tmp1[ loop1][ 3*loop2+loop3 ]
                # Manual scalar product
                tpp = 0
                for i in range(NTrRo):
                    pp = S1[loop1][i] * C0[i][loop2][loop3]    # Corrected
                    tpp = tpp + pp
                Tmp1[loop1][ 3*loop2+loop3 ] = tpp

    # Initialize Tmp2
    Tmp2 = []
    for i in range(nat3):
        Tmp2.append( [] )
    for i in range(nat3):
        for j in range(nat3):
            Tmp2[i].append( 0 )

    # calculate Tmp2
    for loop2 in range(nat3):
        for loop1 in range(natom):
            for loop0 in range(3):
                #Tmp2[loop1][3*loop2+loop3]
                # Manual scalar product
                tpp = 0
                for i in range(NTrRo):
                    pp = dC0[i][loop1][loop0] * Tmp1[i][loop2]
                    tpp = tpp + pp
                Tmp2[loop1*3+loop0][loop2] = tpp

    # Initialize R
    R = []
    for i in range(nat3):
        R.append( [] )
    for i in range(nat3):
        for j in range(nat3):
            R[i].append( 0 )

    # calculate R
    for loop2 in range(nat3):
        for loop1 in range(natom):
            for loop0 in range(3):
                #R[3*loop1+loop0][loop2]
                # Manual scalar product
                tpp = 0
                for i in range(NTrRo):
                    pp = C0[i][loop1][loop0] * Tmp1[i][loop2]
                    tpp = tpp + pp
                R[3*loop1+loop0][loop2] = tpp


    # Initialize Tmp3
    Tmp3 = []
    for i in range(nat3):
        Tmp3.append( [] )

    for i in range(nat3):
        for j in range(nat3):
            Tmp3[i].append( 0 )

    # calculate Tmp3
    for loop1 in range(nat3):
        for loop2 in range(nat3):
                #Tmp3[loop1][loop2]
                # manual scalar product
            tpp = 0
            for i in range(nat3):
                pp = R[loop1][i] * Tmp2[i][loop2]
                tpp = tpp + pp
            Tmp3[loop1][loop2] = tpp

    # Initialize Q
    Q = []
    for i in range(nat3):
        Q.append( [] )
    for i in range(nat3):
        for j in range(nat3):
            Q[i].append( 0 )

    for i in range(nat3):
        for j in range(nat3):
            Q[i][j] = Tmp2[i][j] - Tmp3[i][j]

    # Initialize dR
    dR = []
    for i in range(nat3):
        dR.append( [] )
    for i in range(nat3):
        for j in range(nat3):
            dR[i].append( 0 )

    for i in range(nat3):
        for j in range(nat3):
            dR[i][j] = Q[i][j] + Q[j][i]

    return (R,dR)



# CalPdP for calculate coriolis coupling coefficient
def CalPdP(natom, AtMass, R, dR, FXMw, RoMw, s):

    # FXMw is eta -> dimension correct
    # RoMw is curvature -> dimension correct




    nat3 = 3*natom

    # Initialize P
    P = []
    for i in range(nat3):
        P.append( [] )
    for i in range(nat3):
        for j in range(nat3):
            P[i].append( 0 )

    # calculate P
    for i in range(nat3):
        for j in range(nat3):
            P[j][i] = R[j][i] + FXMw[i] * FXMw[j] # might lead to bug

    # Initialize dP
    dP = []
    for i in range(nat3):
        dP.append( [] )
    for i in range(nat3):
        for j in range(nat3):
            dP[i].append( 0 )

    # calculate dP
    for i in range(nat3):
        for j in range(nat3):
            if s <= 0:
                dP[j][i] = dR[j][i] + FXMw[i] * RoMw[j] + RoMw[i] * FXMw[j] # the minus sign might lead to bugs
                                                                             # I have changed '-' signs into '+' signs
            else:
                dP[j][i] = dR[j][i] - FXMw[i] * RoMw[j] - RoMw[i] * FXMw[j]
#
    return (P,dP)




# CalKdK for coriolis coupling
def CalKdK(natom,AtMass,FFX,P,dP):
    nat3 = 3*natom


    # FFXMw should have the dimension of (nat3,nat3)
    # However, input FFX  might be the lower triangular form of Hessian (un-mass-weighted)

    FFXMw = LT2Sqr(nat3,FFX)

    # we need mass-weighting for FFXMw
    for j in range(nat3):
        for k in range(nat3):
            jm = j / 3
            km = k / 3
            FFXMw[j][k] = FFXMw[j][k] /\
                          ( (AtMass[jm])**0.5 * (AtMass[km])**0.5 )

    # Initialize P1
    P1 = []
    for i in range(nat3):
        P1.append( [] )
    for i in range(nat3):
        for j in range(nat3):
            P1[i].append( 0 )

    # calculate P1
    for i in range(nat3):
        for j in range(nat3):
            if i == j:
                P1[j][i] = 1.0 - P[j][i]
            else:
                P1[j][i] = -1 * P[j][i]


    # Initialize Tmp1, Tmp2 and Tmp3
    # Initialize FFXP,dFFXP
    Tmp1 = []
    Tmp2 = []
    Tmp3 = []
    FFXP = []
    dFFXP = []


    for i in range(nat3):
        Tmp1.append( [] )
        Tmp2.append( [] )
        Tmp3.append( [] )
        FFXP.append( [] )
        dFFXP.append( [] )
    for i in range(nat3):
        for j in range(nat3):
            Tmp1[i].append( 0 )
            Tmp2[i].append( 0 )
            Tmp3[i].append( 0 )
            FFXP[i].append( 0 )
            dFFXP[i].append( 0 )

    # calculate Tmp1
    for loop1 in range(nat3):  # corrected
        for loop2 in range(nat3):  # corrected
            Tmp1[loop1][loop2] = np.dot( P1[loop1],FFXMw[loop2] )

    # calculate FFXP
    for loop1 in range(nat3):
        for loop2 in range(nat3):
            tpp = 0
            for i in range(nat3):
                pp = Tmp1[loop1][i] * P1[i][loop2]
                tpp = tpp + pp
            FFXP[loop1][loop2] = tpp

    # calculate Tmp2
    for loop1 in range(nat3):
        for loop2 in range(nat3):
            Tmp2[loop1][loop2] = np.dot( dP[loop1], FFXMw[loop2] )

    # calculate dFFXP - step 1
    for loop1 in range(nat3):
        for loop2 in range(nat3):
            tpp = 0
            for i in range(nat3):
                pp = -1 * Tmp2[loop1][i] * P1[i][loop2]
                tpp = tpp + pp
            dFFXP[loop1][loop2] = tpp

    # calculate Tmp3
    for loop1 in range(nat3):
        for loop2 in range(nat3):
            Tmp3[loop1][loop2] = np.dot( P1[loop1], FFXMw[loop2] )

    # calculate dFFXP - step 2
    for loop1 in range(nat3):
        for loop2 in range(nat3):
            tpp = 0
            for i in range(nat3):
                pp = Tmp3[loop1][i] * dP[i][loop2]
                tpp = tpp + pp
            dFFXP[loop1][loop2] = dFFXP[loop1][loop2] - tpp

    #
    return (FFXP,dFFXP)



# CplNM1 for coriolis coupling
def CplNM1(natom, AtMass, DX, FFXP, dFFXP):
    nat3 = 3*natom
    Factor=26424608.1937383
    Tresh1=1.0e-3

    nvib = len(DX)

    # Initialize Tmp1,Tmp2,Tmp3,Tmp4,
    # Initialize B

    Tmp1 = []
    Tmp2 = []
    Tmp3 = []
    Tmp4 = []
    Frq2 = []
    GenFreqs = []
    B = []


    for i in range(nvib):
        Tmp1.append( [] )
        Tmp2.append( [] )
        Tmp3.append( [] )
        Tmp4.append( [] )
        Frq2.append( [] )
        GenFreqs.append( [] )
        B.append( [] )


    for i in range(nvib):
        for j in range(nat3):
            Tmp1[i].append( 0 )
            Tmp3[i].append( 0 )

        for j in range(nvib):
            Tmp2[i].append( 0 )
            Tmp4[i].append( 0 )
            B[i].append( 0 )

    # calculate Tmp1
    for loop1 in range(nvib):
        for loop2 in range(nat3):
            tpp = 0
            for i in range(nat3):
                pp = DX[loop1][i] * FFXP[loop2][i]
                tpp = tpp + pp
            Tmp1[loop1][loop2] = tpp

    # calculate Tmp2
    for loop1 in range(nvib):
        for loop2 in range(nvib):
            Tmp2[loop2][loop1] = np.dot( Tmp1[loop1], DX[loop2] )

    #
    DelMax = 0.0
    for Mu in range(nvib):
        X = Tmp2[Mu][Mu]
        Frq2[Mu] = X
        Y = Factor * X
        if Y > 0.0:
            Freq2 =  math.sqrt(Y)
        else:
            Freq2 = -math.sqrt(-Y)

        GenFreqs[Mu] = Freq2


    # calculate Tmp3
    for loop1 in range(nvib):
        for loop2 in range(nat3):
            tpp = 0
            for i in range(nat3):
                pp = DX[loop1][i] * dFFXP[i][loop2]
                tpp = tpp + pp
            Tmp3[loop1][loop2] = tpp


    # calculate Tmp4
    for loop1 in range(nvib):
        for loop2 in range(nvib):
            Tmp4[loop1][loop2] = np.dot( Tmp3[loop1], DX[loop2] )


    for Mu in range(nvib):
        for Nu in range(nvib):
            FrqMu = Factor * Frq2[Mu]
            FrqNu = Factor * Frq2[Nu]
            #print FrqMu,FrqNu         # Stuck here
            if FrqMu >= 0:
                FrqMu = math.sqrt(FrqMu)
            else:
                FrqMu = -math.sqrt(-FrqMu)
            if FrqNu >= 0:
                FrqNu = math.sqrt(FrqNu)
            else:
                FrqNu = -math.sqrt(-FrqNu)
            Delta = abs( FrqMu - FrqNu )

            if (Delta <= Tresh1 ) or (Mu == Nu):
                B[Nu][Mu] = 0.0
            else:
                B[Nu][Mu] = Tmp4[Nu][Mu]/(Frq2[Mu] - Frq2[Nu])
                if abs(B[Nu][Mu]) < 1e-8:
                    B[Nu][Mu] = 0.0
                B[Nu][Mu] = abs(B[Nu][Mu])



    #
    return B





# calcoriolis ( nm2,Brlist[j].Vmass,Brlist[j].Vkappa,Brlist[j].Veta,Brlist[j].Vhess,Brlist[j].Vcoor  )
def calcoriolis( DX, AtMass, kappa, eta, FFX, CC,s):

    natom = len(AtMass)
    nat3 = 3*natom


    # space for later correction
#    DX = []
#    for i in range(len(DXun)):
#       v = copy.deepcopy( DXun[i] )
#       for j in range(nat3):
#           jm = j / 3
#           v[j] = v[j] / (AtMass[jm]**0.5)
#       vnorm = np.linalg.norm(v)
#       DX.append( (v/vnorm).tolist() )



    #
    #natom = len(AtMass)
    NTrRo,C0,dC0 = CalCdC(natom,AtMass,CC,eta)
    R,dR = CalQRS(natom, NTrRo,C0,dC0)
    P,dP = CalPdP(natom, AtMass, R, dR, eta, kappa,s)
    FFXP,dFFXP = CalKdK(natom,AtMass,FFX,P,dP)
    B = CplNM1(natom, AtMass, DX, FFXP, dFFXP)

    sys.stdout.write( str(B[3][4]      )+"\n")



# calAvAM( oBlist, parB )
def calAvAM( Brlist, parB, parID, outfile):

    if os.path.isfile(outfile):
        stop("Error: this file already exists- "+outfile)
    fh = open(outfile,"w")

    names = ["s"]
    for i in parID:
        names.append( i[-1] )
    for i in range(len(names)-1):
        fh.write(names[i]+",")
    fh.write(names[-1])
    fh.write('\n')



    AvAMdata = []
    #[
    #[ s, qn_1, qn_2, ... ],
    #[                    ]
    #]


    nparm = len( parB[0] )
    numPn = len(Brlist)

    for i in range(numPn):
        datarow = []
        sval = Brlist[i].s
        natom = len(Brlist[i].Vmass)

        # calculate normal mode
        if sval == 0.0:
            e0,nm0,cdisp0 = solWilsProj( Brlist[i].Vcoor,Brlist[i].Vhess,Brlist[i].Vmass,Brlist[i].Vnm_0 )
        else:
            e0,nm0,cdisp0 = solWilsProj( Brlist[i].Vcoor,Brlist[i].Vhess,Brlist[i].Vmass,Brlist[i].Veta )
        L0 = copy.deepcopy( cdisp0 ) # maybe we need to use cdisp0 instead
        nvib = len(L0)

        # calculate normal mode force constant
        FCNM = []
        #for j in range(nvib):
        #    FCNM.append( 0 )

        FFX = LT2Sqr(3*natom, Brlist[i].Vhess  )

        for j in range(nvib):
            tmp1 = []
            for k in range(3*natom):
                tmp1.append( np.dot( L0[j], FFX[k] ) )
            FCNM.append( np.dot( tmp1, L0[j] ))
        #print FCNM


        # calculate SqMass
        SqMass = []
        for j in range(3*natom):
            SqMass.append( [] )
        for j in range(3*natom):
            for k in range(3*natom):
                SqMass[j].append( 0 )
        Vsqmass = []
        for j in range(natom):
            for k in range(3):
                Vsqmass.append( Brlist[i].Vmass[j] )

        for j in range(3*natom):
            SqMass[j][j] = Vsqmass[j]

        # calculate un-mass weighted kappa
        Vkappa = copy.deepcopy( Brlist[i].Vkappa )
        Lkappa = [Vkappa]
        for j in range(3*natom):
            Lkappa[0][j] = Lkappa[0][j] / math.sqrt( SqMass[j][j] )


#       # loop over all parameters
        for j in range(nparm):
            G = calG(i,j,Brlist,parB) # not sure..
            b = makB(natom,i,j,parB)
            #print b

            # calculate D matrix (vector)
            D = []
            for k in range(nvib):
                D.append( 0 )
            for k in range(nvib):
                D[k] = np.dot( b, L0[k] )

            # calculate Adiabatic Displacement (local mode in normal coordinates)
            AdiabDisp = []
            for k in range(nvib):
                AdiabDisp.append( 0 )

            tmp2 = 0
            for k in range(nvib):
                tmp2 = tmp2 + D[k]*D[k]/FCNM[k]
            for k in range(nvib):
                AdiabDisp[k] = D[k]/FCNM[k]/tmp2

            # calculate Adiabatic mode vector
            Ad = []
            for k in range(3*natom):
                Ad.append( 0 )

            for k in range(3*natom):
                tmp3 = 0
                for l in range(nvib):
                    tmp3 = tmp3 + AdiabDisp[l] * L0[l][k]
                Ad[k] = tmp3

            # put Ad in a list
            VAd = [ Ad ]
            nat3 = 3*natom
            amp1 = AnBInA( nat3, 1, 1, VAd, SqMass, Lkappa, True)
            #print amp1
            datarow.append( amp1[0][0] )

        #
#       print Brlist[i].s,datarow
        fh.write( str(Brlist[i].s)+","  )
        for k in range(len(datarow)-1):
            fh.write( str(datarow[k])+","   )
        fh.write(str(datarow[-1]) +"\n" )

    print(("Output: adiabatic mode coupling coefficient written to "+outfile))





## calKa2-P
def calKa2p( Brlist, parB):
    # This one is the projection version of calKa2
    pass
    au2dy=15.56893

    nparm = len(parB[0])
    numPn = len(Brlist)

    for i in range(numPn):
        sval = Brlist[i].s
        natom = len(Brlist[i].Vmass)

        nvib = 3*natom - 7 #

        e0,nm0,cdisp0 = solWils(Brlist[i].Vcoor, Brlist[i].Vhess, Brlist[i].Vmass)
        L1 = copy.deepcopy( nm0 )  # might need to be changed into cdisp0

        e2,nm2,cartdisp2=solWilsProj(Brlist[i].Vcoor, Brlist[i].Vhess, Brlist[i].Vmass,L1[0])
        L = copy.deepcopy( nm2 )

        f = LT2Sqr(3*natom, Brlist[i].Vhess )

        # Matrix products
        tmp1 = []
        for j in range(nvib):
            tmp1.append( [] )
        for j in range(nvib):
            for k in range(3*natom):
                tmp1[j].append( 0 )

        K = []
        for j in range(nvib):
            K.append( [] )
        for j in range(nvib):
            for m in range(nvib):
                K[j].append( 0 )


        for j in range(nvib):
            for k in range(3*natom):
                tmp1[j][k] = np.dot( L[j], f[k] )

        for j in range(nvib):
            for m in range(nvib):
                K[j][m] = np.dot( L[j], tmp1[m] )

        Kinv = np.linalg.inv(K)
        datarow = []
        for j in range(nparm):
            b = makB(natom,i,j,parB)
            d = []
            for k in range(nvib):
                d.append( np.dot(b, L[k]) )
            tmp1 = Kinv.dot(d)
            gamma = np.dot(d, tmp1)
            tmp3 = 1.0/gamma* au2dy
            datarow.append( tmp3 )

        print(sval, datarow)



## calKa2
def calKa2( Brlist, parB):
    # This one has no major difference from calKaNew
    pass
    au2dy=15.56893

    nparm = len(parB[0])
    numPn = len(Brlist)

    for i in range(numPn):
        sval = Brlist[i].s
        natom = len(Brlist[i].Vmass)
        nvib = 3*natom - 6 #
        e0,nm0,cdisp0 = solWils(Brlist[i].Vcoor, Brlist[i].Vhess, Brlist[i].Vmass)

        L = copy.deepcopy( nm0 )  # might need to be changed into cdisp0

        f = LT2Sqr(3*natom, Brlist[i].Vhess )

        # Matrix products
        tmp1 = []
        for j in range(nvib):
            tmp1.append( [] )
        for j in range(nvib):
            for k in range(3*natom):
                tmp1[j].append( 0 )

        K = []
        for j in range(nvib):
            K.append( [] )
        for j in range(nvib):
            for m in range(nvib):
                K[j].append( 0 )


        for j in range(nvib):
            for k in range(3*natom):
                tmp1[j][k] = np.dot( L[j], f[k] )

        for j in range(nvib):
            for m in range(nvib):
                K[j][m] = np.dot( L[j], tmp1[m] )

        Kinv = np.linalg.inv(K)
        datarow = []
        for j in range(nparm):
            b = makB(natom,i,j,parB)
            d = []
            for k in range(nvib):
                d.append( np.dot(b, L[k]) )
            tmp1 = Kinv.dot(d)
            gamma = np.dot(d, tmp1)
            tmp3 = 1.0/gamma* au2dy
            datarow.append( tmp3 )

        print(sval, datarow)

## calKaNew
def calKaNew( Brlist, parB, parID, outfile):

    if os.path.isfile(outfile):
        stop("Error: this file already exists- "+outfile)
    fh = open(outfile,"w")

    names = ["s"]
    for i in parID:
        names.append( i[-1] )


    for i in range(len(names)-1):
        fh.write( names[i]  )
        fh.write(",")
    fh.write( names[-1] )
    fh.write('\n')



    data = []

    au2dy=15.56893
    au2dy_a=4.359744

    nparm = len( parB[0] )
    numPn = len(Brlist)

    for i in range(numPn):
        sval = Brlist[i].s
        natom = len(Brlist[i].Vmass)
        f = LT2Sqr(3*natom, Brlist[i].Vhess )
        fpsinv =scipy.linalg.pinv2( f )
        #print fpsinv
        #stop ("")

        datarow = []

        for j in range(nparm):
            b = makB(natom,i,j,parB)

            tmp1 = fpsinv.dot(b)
            #print tmp1
            tmp2 = np.dot( tmp1,b)
            tmp3 = 1.0/tmp2 * au2dy
            datarow.append( tmp3 )

        #print sval, datarow
        fh.write( str(sval)+","   )
        for k in range(len(datarow) -1 ):
            fh.write( str( datarow[k]  )+","  )
        fh.write( str(datarow[-1] )+"\n"  )

    print(("Output: adiabatic force constant k^a written to "+outfile))
    print("Warning: above result should be used with caution.")




def TorsDispB(  xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd  ):
    pass

    # A' = [xap,yap,zap  ]

    xap = xapform( xa,ya,za,xb,yb,zb,xc,yc,zc )
    yap = yapform( xa,ya,za,xb,yb,zb,xc,yc,zc )
    zap = zapform( xa,ya,za,xb,yb,zb,xc,yc,zc )

    #print (xap,yap,zap)

    #A'A = [xa-xap, ya-yap, za-zap]
    ApA = [xa-xap, ya-yap, za-zap]


    # D' = [xdp, ydp, zdp]

    xdp = xdpform( xb,yb,zb,xc,yc,zc,xd,yd,zd )
    ydp = ydpform( xb,yb,zb,xc,yc,zc,xd,yd,zd )
    zdp = zdpform( xb,yb,zb,xc,yc,zc,xd,yd,zd )

    #print (xdp,ydp,zdp)

    #D'D = [xd-xdp, yd-ydp, zd-zdp]
    DpD = [xd-xdp, yd-ydp, zd-zdp]


    prod1 =  np.dot(ApA,DpD)


    # Direction on atom A = [ m1,n1,p1 ]

    m1 = m1form2( xa,ya,za,xb,yb,zb,xc,yc,zc )
    n1 = n1form2( xa,ya,za,xb,yb,zb,xc,yc,zc )
    p1 = p1form2( xa,ya,za,xb,yb,zb,xc,yc,zc )

    #print (m1,n1,p1)



    # Direction on atom D = [ m2,n2,p2 ]
    m2 = m2form2( xb,yb,zb,xc,yc,zc,xd,yd,zd )
    n2 = n2form2( xb,yb,zb,xc,yc,zc,xd,yd,zd )
    p2 = p2form2( xb,yb,zb,xc,yc,zc,xd,yd,zd )

    #print (m2,n2,p2)

    prod2 = m1*m2 + n1*n2 + p1*p2

    if prod1 < 0:
        if prod2 > 0:
            pass
        else:
            m2 = -1 *m2
            n2 = -1 *n2
            p2 = -1 *p2

    elif prod1 > 0:
        if prod2 < 0:
            pass
        else:
            m2 = -1 *m2
            n2 = -1 *n2
            p2 = -1 *p2



    bv =  [m1,n1,p1,m2,n2,p2]
    return bv








TorsDispB(0,-1,-2,0,0,-1,0,0,1,0,1,2)













def BendDispB1(xa,ya,za,xb,yb,zb,xc,yc,zc):
    pass
    # B matrix for Bending displacement - 1
    # no elements for atom B


    #           B
    #          / \
    #         /   \
    #        A     C
    #
    #

    # check whether 3 atoms are on the same line
    thresh = 0.001
    col = (ya - yb) * (xa - xc) - (ya - yc) * (xa - xb)
    if col < thresh:
        stop("Error: 3 atoms on the same line from BendDispB1...")


    # mu_A = [ m1, n1, p1 ]

    m1= m1formula( xa,ya,za,xb,yb,zb,xc,yc,zc )
    n1= n1formula( xa,ya,za,xb,yb,zb,xc,yc,zc )
    p1= p1formula( xa,ya,za,xb,yb,zb,xc,yc,zc )

    # mu_C = [ m2, n2, p2 ]

    m2 = m2formula( xa,ya,za,xb,yb,zb,xc,yc,zc )
    n2 = n2formula( xa,ya,za,xb,yb,zb,xc,yc,zc )
    p2 = p2formula( xa,ya,za,xb,yb,zb,xc,yc,zc )


    # check direction of mu_A and mu_C


    print("m1=",m1)
    print("n1=",n1)
    print("p1=",p1)
    print("m2=",m2)
    print("n2=",n2)
    print("p2=",p2)

    # calculate B matrix

    # D(m1, xa)
    #Dm1xa = dm1xaformula( xa,ya,za,xb,yb,zb,xc,yc,zc )
    #Dm1ya =

    #print "D(m1,xa)",Dm1xa


#BendDispB1(-1,-1,2,0,0,2,1,-1,2)




# Find atoms give large contribution to curvature
def findcontri(Brlist,thresh, maxatoms,outfile ):

    numPn = len(Brlist)
    for i in range(numPn):
        tmp = np.linalg.norm( Brlist[i].Vkappa )
        if tmp > thresh:
            sval = Brlist[i].s
            klist = []
            natom = Brlist[i].NAtom
            for j in range(natom):
                t1 = ( Brlist[i].Vkappa[ j*3 + 0 ] ) **2
                t2 = ( Brlist[i].Vkappa[ j*3 + 1 ] ) **2
                t3 = ( Brlist[i].Vkappa[ j*3 + 2 ] ) **2
                klist.append(   (t1+t2+t3)**0.5  )
            #
            klist = np.array( klist )
            ind = sorted(list(range(len(klist))), key=lambda x: klist[x])[-maxatoms:]
            ind2 = []
            for i in range(len(ind)):
                j = len(ind) - 1 - i
                ind2.append( ind[j] )
                ind2[i] = ind2[i] + 1


            print(sval, ind2)




            pass

    pass



# Get the points on IRC path
# GetPoints(Brlist)
def GetPoints(Brlist):
    pass
    data0 = []
    sval = []
    numPn = len(Brlist)
    for i in range(numPn):
        pts = []
        natom = Brlist[i].NAtom
        sval.append( Brlist[i].s )
        for j in range(natom):
            for k in range(3):
                pts.append( Brlist[i].Vcoor[j][k] * (Brlist[i].Vmass[j] )**0.5)
        #print pts
        data0.append( pts )

    #
    data1 = []
    for i in range(3*natom):
        data1.append( [] )
    for i in range(3*natom):
        for j in range(numPn):
            data1[i].append( data0[j][i] )

    #print data1[0],data1[1]
    #print sval

    for i in range(len(sval)):
        #print sval[i]
        #print data1[4][i]
        pass

    return (sval,data1)
