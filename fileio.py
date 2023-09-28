import sys
import struct
import os.path
from scipy.interpolate import CubicSpline
import numpy as np
import copy
import os
import time

def stop(words):
    print (words)
    sys.exit()


def logo():

    print("         _           ___           ___           ___           ___     ")
    print("        /\ \        /\__\         /\  \         /\__\         /\  \    ")
    print("       /  \ \      /:/  /        /::\  \       /:/  /        /::\  \   ")
    print("      / /\ \ \    /:/  /        /:/\:\  \     /:/  /        /:/\:\  \  ")
    print("     / / /\ \_\  /:/  /  ___   /::\~\:\  \   /:/__/  ___   /::\~\:\  \ ")
    print("    / / /_/ / / /:/__/  /\__\ /:/\:\ \:\__\  |:|  | /\__\ /:/\:\ \:\__\ ")
    print("   / / /__\/ /  \:\  \ /:/  / \/_|::\/:/  /  |:|  |/:/  / \/__\:\/:/  /")
    print("  / / /_____/    \:\  /:/  /     |:|::/  /   |:|__/:/  /       \::/  / ")
    print(" / / /            \:\/:/  /      |:|\/__/     \::::/__/        /:/  /  ")
    print("/ / /              \::/  /       |:|  |        ~~~~           /:/  /   ")
    print("\/_/                \/__/         \|__|                       \/__/    ")

    print ("")
    print("     The Standalone PROGRAM of UNIFIED REACTION VALLEY APPROACH ")
    print ("                       Ver. 0.1.61(2017/05/24)       ")
    print ("")
    print ("        Yunwen Tao, Wenli Zou, Elfi Kraka and Dieter Cremer")
    print ("      Computational and Theoretical Chemistry Group (CATCO), SMU  ")
    print ("                      Dallas, Texas 75275 USA                 ")
    print ("               Webpage: sites.smu.edu/dedman/catco/         ")
    print ("")


### Element symbol -> atomic number
def elem2IZ(sym):

    IZ = 0

    data =[
    [1,"H"],
    [2,"He"],
    [3,"Li"],
    [4,"Be"],
    [5,"B"],
    [6,"C"],
    [7,"N"],
    [8,"O"],
    [9,"F"],
    [10,"Ne"],
    [11,"Na"],
    [12,"Mg"],
    [13,"Al"],
    [14,"Si"],
    [15,"P"],
    [16,"Si"],
    [17,"Cl"],
    [18,"Ar"],
    [19,"K"],
    [20,"Ca"],
    [21,"Sc"],
    [22,"Ti"],
    [23,"V"],
    [24,"Cr"],
    [25,"Mn"],
    [26,"Fe"],
    [27,"Co"],
    [28,"Ni"],
    [29,"Cu"],
    [30,"Zn"],
    [31,"Ga"],
    [32,"Ge"],
    [33,"As"],
    [34,"Se"],
    [35,"Br"],
    [36,"Kr"],
    [37,"Rb"],
    [38,"Sr"],
    [39,"Y"],
    [40,"Zr"],
    [41,"Nb"],
    [42,"Mo"],
    [43,"Tc"],
    [44,"Ru"],
    [45,"Rh"],
    [46,"Pd"],
    [47,"Ag"],
    [48,"Cd"],
    [49,"In"],
    [50,"Sn"],
     [51,"Sb"],
     [52,"Te"],
     [53,"I"],
     [54,"Xe"],
     [55,"Cs"],
     [56,"Ba"],
     [57,"La"],
     [58,"Ce"],
     [59,"Pr"],
     [60,"Nd"],
     [61,"Pm"],
     [62,"Sm"],
     [63,"Eu"],
     [64,"Gd"],
     [65,"Tb"],
     [66,"Dy"],
     [67,"Ho"],
     [68,"Er"],
     [69,"Tm"],
     [70,"Yb"],
     [71,"Lu"],
     [72,"Hf"],
     [73,"Ta"],
     [74,"W"],
     [75,"Re"],
     [76,"Os"],
     [77,"Ir"],
     [78,"Pt"],
     [79,"Au"],
     [80,"Hg"],
     [81,"Tl"],
     [82,"Pb"],
     [83,"Bi"],
     [84,"Po"],
     [85,"At"],
     [86,"Rn"],
     [87,"Fr"],
     [88,"Ra"],
     [89,"Ac"],
     [90,"Th"],
     [91,"Pa"],
     [92,"U"],
     [93,"Np"],
     [94,"Pu"],
     [95,"Am"],
     [96,"Cm"],
     [97,"Bk"],
     [98,"Cf"],
     [99,"Es"],
     [100,"Fm"],
     [101,"Md"],
     [102,"No"],
     [103,"Lr"],
     [104,"Rf"],
     [105,"Db"],
     [106,"Sg"],
     [107,"Bh"],
     [108,"Hs"],
     [109,"Mt"]]

    tmp = sym.strip()

    for i in data:
        if i[1] == tmp:
            IZ = i[0]

#   print IZ
    return IZ



#elem2IZ("H")


####
# Atomic number -> mass converter
def massof(IZ):
    masslib=[ \
         1.007825037,   4.002603250,   7.016004500, \
         9.012182500,  11.009305300,  12.000000000, \
        14.003074008,  15.994914640,  18.998403250, \
        19.992439100,  22.989769700,  23.985045000, \
        26.981541300,  27.976928400,  30.973763400, \
        31.972071800,  34.968852729,  39.962383100, \
        38.963707900,  39.962590700,  44.955913600, \
        47.947946700,  50.943962500,  51.940509700, \
        54.938046300,  55.934939300,  58.933197800, \
        57.935347100,  62.929599200,  63.929145400, \
        68.925580900,  73.921178800,  74.921595500, \
        79.916520500,  78.918336100,  83.911506400, \
        84.911700000,  87.905600000,  88.905400000, \
        89.904300000,  92.906000000,  97.905500000, \
        98.906300000, 101.903700000, 102.904800000, \
       105.903200000, 106.905090000, 113.903600000, \
       114.904100000, 117.901800000, 120.903800000, \
       129.906700000, 126.900400000, 131.904200000, \
       132.905429000, 137.905000000, 138.906100000, \
       139.905300000, 140.907400000, 141.907500000, \
       144.912700000, 151.919500000, 152.920900000, \
       157.924100000, 158.925000000, 163.928800000, \
       164.930300000, 165.930400000, 168.934400000, \
       173.939000000, 174.940900000, 179.946800000, \
       180.948000000, 183.951000000, 186.956000000, \
       189.958600000, 192.963300000, 194.964800000, \
       196.966600000, 201.970600000, 204.974500000, \
       207.976600000, 208.980400000, 208.982500000, \
       210.987500000, 222.017500000, 223.019800000, \
       226.025400000, 227.027800000, 232.038200000, \
       231.035900000, 238.050800000, 237.048000000, \
       242.058700000, 243.061400000, 246.067400000, \
       247.070200000, 249.074800000, 252.082900000, \
       252.082700000, 255.090600000, 259.101000000, \
       262.109700000, 261.108700000, 262.114100000, \
       266.121900000, 264.124700000, 265.000000000, \
       268.138800000, 269.000000000, 272.000000000, \
       277.000000000, 280.000000000, 280.000000000, \
       280.000000000, 280.000000000, 280.000000000, \
       280.000000000, 280.000000000, 280.000000000 ]

    tmp = masslib[IZ-1]
    return tmp



def iscXYZ( filename, numberofpoints):
    "Initial Scan of the XYZ trajectory file"

    i = 0
    NP = 0
    with open(filename) as f:
        for line in f:
            i = i + 1
            if i == 1:
                try:
                    tmp1 = int(line)
                except ValueError:
                    print ("Error: Please check the number of atoms in the XYZ trajectory file! ")
                NP = NP + 1
            else:
                if (i-1)%(1+1+tmp1) == 0:
                    try:
                        tmp2 = int(line)
                    except:
                        print ("Error: Please check the number of atoms in the XYZ trajectory file! ")
                    if tmp1 == tmp2:
                        NP = NP + 1

    #print NP
    return NP



#xyzf="/scratch/users/yunwent/softs/cp2k-4.1/tutorials/chorismate-309/f-26/CHORISMATE-309-f26-pos-1.xyz"
#numP = 0
#iscXYZ( xyzf, numP )





##################
# initial scan of new browsing file

def iscNB( filename, numberofpoints):
    "Initial Scan of the new browsing file"

    numbegin=0
    numend=0
    with open(filename) as f:
        for line in f:
            if len(line)>0:
                #print line
                if line.split()[0]=="BEGIN":
                    numbegin = numbegin + 1
                if line.split()[0]=="END":
                    numend = numend + 1
    f.close()
    if numend < numbegin:
        print("Error: Incomplete browsing file.")
        sys.exit()
    else:
        numberofpoints = numend

    return numberofpoints

##################

# initial scan of old browsing file

def iscOB( filename, numberofpoints):
    "Initial Scan of the old browsing file"

    numbegin=0
    numend=0
    with open(filename) as f:
        for line in f:
            if len(line)>0:
                #print line
                if line.split()[0]=="BEGIN":
                    numbegin = numbegin + 1
                if line.split()[0]=="END":
                    numend = numend + 1
    f.close()
    if numend < numbegin:
        print("Error: Incomplete browsing file.")
        sys.exit()
    else:
        numberofpoints = numend

    return numberofpoints

##################




#### parse single XYZ data ####
def parseXYZ(strings):
    "Extract element and coordinates from XYZ data"

    elem = []
    Vcoor = []
    for i in range(len(strings)):
        if i > -1:
            if len(strings[i])>1:
                tmp1 = strings[i].split()[0]
                if tmp1.isdigit() == True:
                    print ("Error: Element symbols expected in XYZ file.")
                    sys.exit()
                el = strings[i].split()[0]
                c1 = float( strings[i].split()[1] ) * 1.889725989 # convert into a.u.
                c2 = float( strings[i].split()[2] ) * 1.889725989 # convert into a.u.
                c3 = float( strings[i].split()[3] ) * 1.889725989 # convert into a.u.
                elem.append( el )
                Vcoor.append( [ c1, c2, c3] )

    return (elem,Vcoor)




###### second read XYZ trajectory file ####
def rdXYZ(filename,Xlist):
    "Read in data from XYZ file"

    natom = 0
    i = 0
    with open( filename ) as f:
        for line in f:
            i = i + 1
            if i == 1:
                natom = int( line )
                break

    numPn = len(Xlist)
    bar = ProgressBar(total = numPn)
    #bar.move()
    #bar.log('We have arrived at: ' + str(i + 1))

    f.close()

    i = 0
    p = -1
    section = []
    # Read in data for each point
    with open( filename ) as f2:
        for line in f2:
            i = i + 1
            if i == 1 or i == 2:
                section = []
                #print line
            else:
                section.append( line )
            if i == (natom+2) :
                #print line
                i = 0
                p = p + 1
                #print section
                Vmass=[]
                elem,Vcoor = parseXYZ( section )
                #print "p",p
                #print len(elem),len(Voor)
                Xlist[p].Vcoor = Vcoor
                for j in range(len(elem)):
                    Vmass.append( massof( elem2IZ( elem[j] ) ) )
                Xlist[p].Vmass = Vmass

                bar.move()
                bar.log('We have arrived at: ' + str(p + 1))
#                sys.stdout.write('.')
#                sys.stdout.flush()

    #print natom
#    sys.stdout.write('\n')




#xyzf="/scratch/users/yunwent/softs/cp2k-4.1/tutorials/chorismate-309/f-26/CHORISMATE-309-f26-pos-1.xyz"
#rdXYZ( xyzf, 1)




####### second read new browsing file #####
def rdNB(filename,newBrowslist):
    "Read in data from new browsing file"

    fieldwidths = (22, 22, 22)
    fmtstring = ' '.join('{}{}'.format(abs(fw), 'x' if fw < 0 else 's')
                            for fw in fieldwidths)
    fieldstruct = struct.Struct(fmtstring)
    parse = fieldstruct.unpack_from

    numpts=len(newBrowslist)
    natom=0
    ct=-1
    irdA=0
    irdB=0
    pmass=0
    numrowmass=0
    tmpmass=[]
    irdC=0
    pcoor=0
    tmpcoor=[]
    irdD=0
    peta=0
    tmpeta=[]
    irdE=0
    pkap=0
    tmpkap=[]
    irdF=0

    with open(filename) as f:
        for i in f:
            #sys.stdout.write('.')
            if len(i) > 0:
                ## Number of atoms
                if i.split()[0]=="Natoms,NatomQ":
                    irdA = 1
                    ct = ct +1  # Counter
                    sys.stdout.write('.')
                    sys.stdout.flush()
                    continue
                if irdA == 1:
                    newBrowslist[ct].NAtom = int(i.split()[0])
                    natom = newBrowslist[ct].NAtom
                    irdA = 0
                ## Atomic mass
                if i.split()[0]=="Atomic":
                    if i.split()[1]=="masses":
                        irdB = 1
                        pmass=1
                        numrowmass = natom / 3
                        if natom % 3 >0:
                            numrowmass = numrowmass + 1   # Problems
                        continue
                if irdB == 1:
                    if pmass < numrowmass:
                        tmpmass.append(i.split())
                        pmass = pmass + 1
                    else:
                        j = natom % 3    # What if j == 0

                        tmp1=[]

                        if j!=0:
                            for k in range(j):
                                tmp1.append(i.split()[k])
                        else:
                            tmp1 = i.split()

                        tmpmass.append(tmp1)


                        irdB = 0
                        pmass = 0
                        tmpmass = [item for sublist in tmpmass for item in sublist]
                        tmpmass = [float(i) for i in tmpmass]
                        newBrowslist[ct].Vmass = tmpmass
                        tmpmass = []
                ## Cartesian coordinates
                if i.split()[0]=="CC":
                    irdC = 1
                    pcoor=1
                    continue
                if irdC == 1:
                    if pcoor < natom:
                        #print i
                        fields = parse(i.encode())
                        #print format(fields)
                        tmp2=[]
                        tmp2.append( float( fields[0].decode() ) )
                        tmp2.append( float( fields[1].decode() ) )
                        tmp2.append( float( fields[2].decode() ) )
                        #print tmp2
                        tmpcoor.append(tmp2)
                        pcoor = pcoor + 1
                    else:
                        #print i
                        fields = parse(i.encode())
                        #print format(fields)
                        #tmp2=list(format(fields))
                        #newBrowslist[ct].Vcoor = tmpcoor
                        tmp2=[]
                        tmp2.append( float( fields[0].decode() ) )
                        tmp2.append( float( fields[1].decode() ) )
                        tmp2.append( float( fields[2].decode() ) )
                        #print tmp2
                        tmpcoor.append(tmp2)
                        newBrowslist[ct].Vcoor = tmpcoor
                        irdC = 0
                        pcoor = 0
                        tmpcoor = []

                ## Path direction vector (Eta)[mass-weighted]
                if i.split()[0]=="Tangent":
                    irdD = 1
                    peta=1
                    continue
                if irdD == 1:
                    if peta < natom:
                        fields = parse(i.encode())
                        tmp3=[]
                        tmp3.append( float( fields[0].decode() ) )
                        tmp3.append( float( fields[1].decode() ) )
                        tmp3.append( float( fields[2].decode() ) )
                        tmpeta.append(tmp3)
                        peta = peta + 1
                    else:
                        fields = parse(i.encode())
                        tmp3=[]
                        tmp3.append( float( fields[0].decode() ) )
                        tmp3.append( float( fields[1].decode() ) )
                        tmp3.append( float( fields[2].decode() ) )
                        tmpeta.append(tmp3)
                        newBrowslist[ct].Veta = tmpeta
                        irdD = 0
                        peta = 0
                        tmpeta = []
                ## Curvature vector (Kappa)[mass-weighted]
                if i.split()[0]=="Curvature":
                    irdE = 1
                    pkap=1
                    continue
                if irdE == 1:
                    if pkap < natom:
                        fields = parse(i.encode())
                        tmp4=[]
                        tmp4.append( float( fields[0].decode() ) )
                        tmp4.append( float( fields[1].decode() ) )
                        tmp4.append( float( fields[2].decode() ) )
                        tmpkap.append(tmp4)
                        pkap = pkap + 1
                    else:
                        fields = parse(i.encode())
                        tmp4=[]
                        tmp4.append( float( fields[0].decode() ) )
                        tmp4.append( float( fields[1].decode() ) )
                        tmp4.append( float( fields[2].decode() ) )
                        tmpkap.append(tmp4)
                        newBrowslist[ct].Vkappa = tmpkap
                        irdE = 0
                        pkap = 0
                        tmpkap = []
                ## Energy, reaction parameter(s)
                if i.split()[0]=="IPoCou,Energy,XXIRC":
                    irdF = 1
                    continue
                if irdF == 1:
                    engy=float(i.split()[1])
                    rps =float(i.split()[2])  # Reaction Path S value
                    #print engy,rps
                    newBrowslist[ct].energy=engy
                    newBrowslist[ct].s     =rps
                    irdF=0
    f.close()
    sys.stdout.write('\n')

######################################################


####### second read old browsing file #####
def rdOB(filename,oldBrowslist):
    "Read in data from new browsing file"

    fieldwidths = (22, 22, 22)
    fmtstring = ' '.join('{}{}'.format(abs(fw), 'x' if fw < 0 else 's')
                            for fw in fieldwidths)
    fieldstruct = struct.Struct(fmtstring)
    parse = fieldstruct.unpack_from

    numpts=len(oldBrowslist)
    natom=0
    ct=-1
    irdA=0
    irdB=0
    pmass=0
    numrowmass=0
    tmpmass=[]
    irdC=0
    pcoor=0
    tmpcoor=[]
    irdD=0
    peta=0
    tmpeta=[]
    irdE=0
    pkap=0
    tmpkap=[]
    irdF=0

    with open(filename) as f:
        for i in f:
            #sys.stdout.write('.')
            if len(i) > 0:

                ## Number of atoms
                if i.split()[0]=="NZ,NSubs":
                    irdA = 1
                    ct = ct +1  # Counter
                    sys.stdout.write('.')
                    sys.stdout.flush()
                    continue
                if irdA == 1:
                    oldBrowslist[ct].NAtom = int(i.split()[0])
                    natom = oldBrowslist[ct].NAtom
                    irdA = 0

#                ## Atomic mass
                if i.split()[0]=="IAnZ,IZ1,IZ2,IZ3,IZ4,LBl,LAlpha,LBeta":
                    irdB = 1
                    pmass=1
                    numrowmass = natom
                    continue
                if irdB == 1:
                    if pmass < numrowmass:
                        IZ = int(i.split()[0])
                        tmpmass.append( massof(IZ) )
                        pmass = pmass + 1
                    else:
                        IZ = int(i.split()[0])
                        tmpmass.append( massof(IZ) )
                        irdB = 0
                        pmass = 0
                        #tmpmass = [item for sublist in tmpmass for item in sublist]
                        #tmpmass = [float(i) for i in tmpmass]
                        oldBrowslist[ct].Vmass = tmpmass
                        tmpmass = []


                ## Cartesian coordinates
                if i.split()[0]=="CC":
                    irdC = 1
                    pcoor=1
                    continue
                if irdC == 1:
                    if pcoor < natom:
                    #print i
                        fields = parse(i.encode())
                        #print(fields)
                        tmp2=[]
                        tmp2.append( float( fields[0].decode() ) )
                        tmp2.append( float( fields[1].decode() ) )
                        tmp2.append( float( fields[2].decode() ) )
                        #print tmp2
                        tmpcoor.append(tmp2)
                        pcoor = pcoor + 1
                    else:
                        #print i
                        fields = parse(i.encode())
                        #print format(fields)
                        #tmp2=list(format(fields))
                        #newBrowslist[ct].Vcoor = tmpcoor
                        tmp2=[]
                        tmp2.append( float( fields[0].decode() ) )
                        tmp2.append( float( fields[1].decode() ) )
                        tmp2.append( float( fields[2].decode() ) )
                        #print tmp2
                        tmpcoor.append(tmp2)
                        oldBrowslist[ct].Vcoor = tmpcoor
                        irdC = 0
                        pcoor = 0
                        tmpcoor = []

                ## Energy gradient [unmass-weighted]
                if i.split()[0]=="FX_ZMat_Orientation":
                    irdD = 1
                    peta=1
                    continue
                if irdD == 1:
                    if peta < natom:
                        fields = parse(i.encode())
                        tmp3=[]
                        tmp3.append( float( fields[0].decode() ) )
                        tmp3.append( float( fields[1].decode() ) )
                        tmp3.append( float( fields[2].decode() ) )
                        tmpeta.append(tmp3)
                        peta = peta + 1
                    else:
                        fields = parse(i.encode())
                        tmp3=[]
                        tmp3.append( float( fields[0].decode() ) )
                        tmp3.append( float( fields[1].decode() ) )
                        tmp3.append( float( fields[2].decode() ) )
                        tmpeta.append(tmp3)
                        tmpeta = [item for sublist in tmpeta for item in sublist]

                        oldBrowslist[ct].Vgrad = tmpeta
                        irdD = 0
                        peta = 0
                        tmpeta = []
                ## Hessian [unmass-weighted]
                if i.split()[0]=="FFX_ZMat_Orientation":
                    irdE = 1
                    nrowhess = (natom*(3*natom+1))/2
                    pkap=1
                    continue
                if irdE == 1:
                    if pkap < nrowhess: # corrected
                        fields = parse(i.encode())
                        tmp4=[]
                        tmp4.append( float( fields[0].decode() ) )
                        tmp4.append( float( fields[1].decode() ) )
                        tmp4.append( float( fields[2].decode() ) )
                        tmpkap.append(tmp4)
                        pkap = pkap + 1
                    else:
                        fields = parse(i.encode())
                        tmp4=[]
                        tmp4.append( float( fields[0].decode() ) )
                        tmp4.append( float( fields[1].decode() ) )
                        tmp4.append( float( fields[2].decode() ) )
                        tmpkap.append(tmp4)
                        tmpkap = [item for sublist in tmpkap for item in sublist]
                        oldBrowslist[ct].Vhess = tmpkap
                        irdE = 0
                        pkap = 0
                        tmpkap = []

                ## Energy, reaction parameter(s)
                if i.split()[0]=="IPoCou,Energy,XXIRC":
                    irdF = 1
                    continue
                if irdF == 1:
                    engy=float(i.split()[1])
                    rps =float(i.split()[2])
                    #print engy,rps
                    oldBrowslist[ct].energy=engy
                    oldBrowslist[ct].s     =rps
                    irdF=0
    f.close()
    sys.stdout.write('\n')

######################################################



## Read in user input file, take the atomic index for parameters
def rdPar(inputfile,parID):
    "Read in the user input and get the atomic index for parameters"
    iopar=0
    irdA=0
    parID=[]
    with open(inputfile) as f:
        for i in f:
            if len(i) > 0:
                if i.split("=")[0][0]=="@":
                    if i.split("=")[0][1:5]=="PARM":
                        if i.split("=")[1].lower().strip()=="on":
                            #print ("Parameter Analysis Enabled.")  # this part seems useless...
                            iopar = 1
                        else:
                            iopar = 0
                            #print ("Parameter Analysis Disabled.")
                #print i.split()
                if len(i.split())>0:
                    if (i.split()[0]=="PARAMETER"):
                        irdA = 1
                        continue
                if irdA == 1:
                    tmp=i.split()

                    if tmp[0]=="ring": # RING COORDINATES
                        tmp1=i.split(":")[0].split()[1]  # Number of atoms in a ring
                        try:
                            tmp1 = int(tmp1)
                        except ValueError:
                            print(("Error: Please check the number of atoms in a ring! \n"+i.strip() ))
                            sys.exit()
                        if tmp1 <= 0:
                            print(("Error: Please check the number of atoms in a ring! \n"+i.strip() ))
                            sys.exit()
                        tmp2 = i.split("(")[1].split(")")[0].split() # Ring specification
        #           for p2 in range(len(tmp2)):
        #                tmp2[p2] = int(tmp2[p2))
        #                print tmp2
                        for c1 in tmp2:
                            if c1.isdigit() is False:
                                print(("Error: Please check the ring specification! \n"+i.strip() ))
                                sys.exit()
                        if len(list(set(tmp2))) != tmp1:
                            print(("Error: Please check the ring specification! \n"+i.strip() ))
                            sys.exit()
                        for p2 in range(len(tmp2)):
                            tmp2[p2] = int(tmp2[p2])
        #           print tmp2

                        if ("[" in i) and ("]" in i):
                            tmp3 = i.split("[")[1].split("]")[0].split()
                        else:
                            print ("Error: Please specify ring parameters! \n")
                            sys.exit()
                            #print tmp3
                        for c1 in tmp3:
                            if c1.isdigit() is False:
                                print(("Error: Please check the ring coordinate parameter!\n"+i.strip()))
                                sys.exit()
                        for c2 in range(len(tmp3)):
                            tmp3[c2] = int( tmp3[c2] )

                        if len(tmp3) < 2 or len(tmp3) > 3:
                            print(("Error: Please check the ring coordinate parameter!\n"+i.strip()))
                            sys.exit()
                        if tmp3[0] < 0 or tmp3[0] > 4:
                            print(("Error: Please check type of ring parameter!\n "+str(tmp3[0]) ))
                            sys.exit()
                        if tmp3[0] == 1:
                            if tmp3[1] < 1 or tmp3[1] > (tmp1-2):
                                print(("Error: Please check parameter for planar deformation amplitude!\n "+i.strip()))
                                sys.exit()
                        if tmp3[0] == 2:
                            if tmp3[1] < 1 or tmp3[1] > (tmp1-2):
                                print(("Error: Please check parameter for planar deformation phase angle!\n "+i.strip()))
                                sys.exit()
                        if tmp3[0] == 3:
                            if tmp1 % 2 == 1: # corrected
                                if tmp3[1] < 2 or tmp3[1] > ((tmp1-1)/2):
                                    print(("Error: Please check parameter for puckering amplitude!\n"+i.strip()))
                                    sys.exit()
                            else:
                                if tmp3[1] < 2 or tmp3[1] > (tmp1/2):
                                    print(("Error: Please check parameter for puckering amplitude!\n"+i.strip()))
                                    sys.exit()
                        if tmp3[0] == 4:
                            if tmp1 % 2 == 1:
                                if tmp3[1] < 2 or tmp3[1] > ((tmp1-1)/2):
                                    print(("Error: Please check parameter for puckering phase angle!\n"+i.strip()))
                                    sys.exit()
                            else:
                                if tmp3[1] < 2 or tmp3[1] > ( (tmp1/2) -1): # A little bit different
                                    print(("Error: Please check parameter for puckering phase angle!\n"+i.strip()))
                                    sys.exit()
                        toput=[]
                        toput.append("g")
                        bname = i.split(":")[1].strip()
                        #print bname
                        toput.append(tmp1)
                        toput.append(tmp2)
                        toput.append(tmp3)
                        toput.append(bname)
                        #print toput
                        parID.append(toput)




                    if tmp[0]=="oop": # Out-of-plane angle
                        tmp1=i.split(":")[0].split()
                        tmp2=[]
                        tmp2.append("o")
                        tmp2.append( int(tmp1[1]))
                        tmp2.append( int(tmp1[2]))
                        tmp2.append( int(tmp1[3]))
                        tmp2.append( int(tmp1[4]))
                        bname=i.split(":")[1].strip()
                        tmp2.append(bname)
                        #print tmp2
                        parID.append(tmp2)

                    if tmp[0]=="pyr": # Pyramidalization angle
                        tmp1=i.split(":")[0].split()
                        tmp2=[]
                        tmp2.append("y")
                        tmp2.append( int(tmp1[1]))
                        tmp2.append( int(tmp1[2]))
                        tmp2.append( int(tmp1[3]))
                        tmp2.append( int(tmp1[4]))
                        bname=i.split(":")[1].strip()
                        tmp2.append(bname)
                        #print tmp2
                        parID.append(tmp2)



                    if tmp[0]=="std":
                        tmp1=i.split(":")[0].split()
                        #print tmp1
                        if len(tmp1) == 3:#(bond)
                            tmp2=[]
                            tmp2.append("r")
                            tmp2.append( int(tmp1[1]) )
                            tmp2.append( int(tmp1[2]) )
                            bname=i.split(":")[1].strip()
                            tmp2.append(bname)
                            #print tmp2
                            parID.append(tmp2)
                        if len(tmp1) == 4:#(bond angle)
                            tmp2=[]
                            tmp2.append("a")
                            tmp2.append( int(tmp1[1]) )
                            tmp2.append( int(tmp1[2]) )
                            tmp2.append( int(tmp1[3]) )
                            bname=i.split(":")[1].strip()
                            tmp2.append(bname)
                            #print tmp2
                            parID.append(tmp2)
                        if len(tmp1) == 5:#(dihedral)
                            tmp2=[]
                            tmp2.append("d")
                            tmp2.append( int(tmp1[1]) )
                            tmp2.append( int(tmp1[2]) )
                            tmp2.append( int(tmp1[3]) )
                            tmp2.append( int(tmp1[4]) )
                            bname=i.split(":")[1].strip()
                            tmp2.append(bname)
                            #print tmp2
                            parID.append(tmp2)



                if len(i.split())>0:
                    if i.split()[0]=="END":
                        irdA = 0

    f.close()
    return parID
#################################


## Print parameter value data to an external file
def printqn2file(parVal,Brlist,filename):

    if os.path.isfile(filename):
        stop("Error: this file already exists- "+filename)
    fh = open(filename,"w")

    numPn = len(parVal)
    nparm = len(parVal[0])

    names = ["s"]
    for i in parVal[0]:
        names.append( i[1] )

    fh.write("=========================\n")
    fh.write("| PARAMETER VALUE TABLE |\n")
    fh.write("=========================\n")

    for i in range(len(names) -1):
        fh.write( names[i] )
        fh.write( ',' )
    fh.write( names[-1] )
    fh.write( '\n' )

    for i in range(numPn):
        fh.write( str(Brlist[i].s) )
        fh.write( ',')
        for j in range(len(parVal[i]) - 1):
            fh.write(  str(parVal[i][j][-1]    ) )
            fh.write( ',' )
        fh.write( str(parVal[i][-1][-1]    ) )
        fh.write( "\n" )

    #
    print(("Output: q_n values written to "+filename))


### Print parameter value data as to the terminal ###

def printqn(parVal,Brlist):
    numPn = len(parVal)

    nparm = len(parVal[0])

    names = ["s"]
    for i in parVal[0]:
        names.append( i[1] )
    print("=========================")
    print("| PARAMETER VALUE TABLE |")
    print("=========================")

    for i in range(len(names) -1):
        sys.stdout.write(names[i])
        sys.stdout.write(',')
    sys.stdout.write(names[-1])
    sys.stdout.write('\n')


    for i in range(numPn):
        sys.stdout.write( str(Brlist[i].s) )
        sys.stdout.write(',')
        for j in range(len(parVal[i]) - 1):
            sys.stdout.write( str(parVal[i][j][-1]    ) )
            sys.stdout.write(',')
        sys.stdout.write( str(parVal[i][-1][-1]    )   )
        sys.stdout.write('\n')













#############################
import  time

class ProgressBar:
    def __init__(self, count = 0, total = 0, width = 50):
        self.count = count
        self.total = total
        self.width = width
    def move(self):
        self.count += 1
    def log(self, s):
        sys.stdout.write(' ' * (self.width + 9) + '\r')
        sys.stdout.flush()
        #print s
        progress = self.width * self.count / self.total
        sys.stdout.write('{0:3}/{1:3}: '.format(self.count, self.total))
        sys.stdout.write('>' * progress + '-' * (self.width - progress) + '\r')
        if progress == self.width:
            sys.stdout.write('\n')
        sys.stdout.flush()


#bar = ProgressBar(total = 10)
#for i in range(10):
#    bar.move()
#    bar.log('We have arrived at: ' + str(i + 1))
#    time.sleep(1)




def getbykey(kee,lst):
    value = []
    if isinstance(lst,list):
        if isinstance(lst[0],list):
            pass
        else:
            stop("Error: function \'getbykey\' expects a list of lists")

        for i in range(len(lst)):
            if lst[i][0] == kee:
                value.append( lst[i][1] )
    else:
        stop("Error: function \'getbykey\' expects a list...")

    if len(value) == 0:
        stop("Error: nothing is found by  \'getbykey\' ...")
    elif len(value) > 1:
        stop("Error: multiple hits found in \'getbykey\' ... ")

    return value





def smergel(smallf,largef,outf):
    # Merge small portion of data into a large data

    if os.path.isfile(outf):
        stop("Error: this file already exists- "+outf)
    if os.path.isfile(smallf):
        pass
    else:
        stop("Error: this file does not exist- "+smallf)

    if os.path.isfile(largef):
        pass
    else:
        stop("Error: this file does not exist- "+largef)

    smalla = []
    smallv = []
    largea = []
    largev = []
    labels = []

    outv = []

    with open(smallf) as f:
        for line in f:
            if len(line) >2:
                a = float( line.split()[0] )
                b = float( line.split()[1] )
                smallv.append([a,b])
                smalla.append(a)

    with open(largef) as f:
        for line in f:
            if len(line) >2:
                a = float( line.split()[0] )
                b = float( line.split()[1] )
                largev.append([a,b])
                largea.append(a)

    for i in largea:
        labels.append(0)

    #
    for i in range(len(smalla)):
        if smalla[i] in largea:
            pass
        else:
            stop("Error: Grid point not match- "+smallf+" - "+largef)

    for i in range(len(largea)):
        if largea[i] in smalla:
            labels[i] = 1


    for i in range(len(largea)):
        if labels[i] == 0:
            outv.append( largev[i] )
        else:
            a = largea[i]
            bl = getbykey(a, smallv)
            b = bl[0]
            outv.append( [a,b] )

    #return outv
    expList2File(outv,outf)


def expList2File(lst1,filename):
    number = len(lst1)
    if number == 0:
        stop("Error: empty list...")

    a = lst1[0]
    nestedlist = 0
    if isinstance( a,list ) == True:
        nestedlist = 1

    # check file exists or not
    if os.path.isfile(filename):
        stop("Error: this file already exists- "+filename)

    # start to write
    fh = open(filename,"w")
    if nestedlist == 0:
        for unit in lst1:
            fh.write( str(unit)+"\n" )
    else:
        for unit in lst1:
            for subunit in unit:
                fh.write( str(subunit) + " ")
            fh.write("\n")

    fh.close()




    #b = 0
    #for i in range(len(lst1)):
    #   a = isinstance( lst1[i], list)
    #    if a == True
    #isinstance( lst1, list )

    #pass


def catcover():
    logfile = '/users/yunwent/share/pURVA_log-01.txt'

    key = -1
    macname = os.uname()[1]
    if macname[0:3] == "mfc":
        key = 1
    if macname[0:7] == "mflogin":
        key = 1

    if key == -1:
        stop("Error: Copyright check failed.....")

    path = os.getcwd()
    userhome = os.environ['HOME']
    hostname = macname
    timestamp = "TIME "+time.strftime("%c")

    if not os.path.isfile(logfile):
        stop("Error: Unknown error.... ")
    else:
        if not os.access(logfile, os.W_OK):
            stop("Error: Unkown error.... ")


        fh = open(logfile, "a")
        fh.write("\n")
        fh.write("USER SOURCE: "+userhome+"\n")
        fh.write("MACHINE ID: "+hostname+"\n")
        fh.write("EXECUTION FOLDER: "+path+"\n")
        fh.write(timestamp+"\n\n")

    pass



def rdcmd(inpf,cmdgui): #,oldnew,nBrowsf,oBrowsf,xyzf, \


    # Verification functionality in CATCO group
    #key = catcover()


    #Ln, Rn, \
    #smthstepsize,smthleftn,smthrighn,smthcurthre,\
    #rmspkcuthigh,rmspkpercent,rmspkgratio ):

    # inpf -> user input file
    # cmdgui -> object including user inputs
    oldnew = -1
    nBrowsf = "InputDataFile"
    oBrowsf = "InputDataFile"
    xyzf    = "InputDataFile"


    # Return user input list
    # 0 -> Datfiletype
    # 1 -> Datfilename
    # 2 -> TFParm =(No, GeomOnly, All)
    # 3 -> TFvib

    # 4 -> TFcurvcorts
    # 5 -> NcurvcorLn
    # 6 -> NcurvcorRn

    # 7 -> TFautosmth
    # 8 -> NautosmthLn
    # 9 -> NautosmthRn
    #10 -> Nautosmthstpsize
    #11 -> Nautosmthd2y

    #12 -> TFremovespk
    #13 -> Nrmspkcut
    #14 -> Nrmspkperc
    #15 -> Nrmspkgratio
    #

    #16 -> TFdircurv
    #17 -> TFavam
    #18 -> TFcurvcpl
    #19 -> TFcoriolis
    #

    #20 -> TFenergy
    #21 -> TFadiabfc
    #

    #22 -> Basename



    #ADIABFC
    with open(inpf) as f:
        for line in f:
            if line[0] == "@":
                if "ADIABFC" in line:
                    if line.split()[1] == "=":
                        if line.split()[2] == "on":
                            tfadiabfc = 1
                            cmdgui.TFadiabfc = tfadiabfc
                        elif line.split()[2] == "off":
                            tfadiabfc = 0
                            cmdgui.TFadiabfc = tfadiabfc
                        else:
                            stop("Error: Please specify ADIABFC correctly in the input file.")

    f.close()
    #cmdgui.TFadiabfc = tfadiabfc


    # ENERGY
    with open(inpf) as f:
        for line in f:
            if len(line) > 4:
                if line[0] == "@":
                    if "ENERGY" in line:
                        if line.split()[1] == "=":
                            if line.split()[2] == "on":
                                tfenergy = 1

                                cmdgui.TFenergy = tfenergy #
                            elif line.split()[2] == "off":
                                tfenergy = 0
                                cmdgui.TFenergy = tfenergy
                            else:
                                stop("Error: Please specify ENERGY correctly in the input file.")

    f.close()
    #cmdgui.TFenergy = tfenergy


    # CORIOLIS
    with open(inpf) as f:
        for line in f:
            if len(line) > 4:
                if line[0] == "@":
                    if "CORIOLIS" in line:
                        if line.split()[1] == "=":
                            if line.split()[2] == "on":
                                tfcoriolis = 1
                            elif line.split()[2] == "off":
                                tfcoriolis = 0
                            else:
                                stop("Error: Please specify CORIOLIS correctly in the input file.")

    f.close()

    cmdgui.TFcoriolis = tfcoriolis
    #print "coriolis",tfcoriolis



    # CURVCPL
    with open(inpf) as f:
        for line in f:
            if len(line) > 4:
                if line[0] == "@":
                    if "CURVCPL" in line:
                        if line.split()[1] == "=":
                            if line.split()[2] == "on":
                                tfcurvcpl = 1
                            elif line.split()[2] == "off":
                                tfcurvcpl = 0
                            else:
                                stop("Error: Please specify CURVCPL correctly in the input file.")

    f.close()

    cmdgui.TFcurvcpl = tfcurvcpl
    #print "curvature coupling", tfcurvcpl





    # AVAM
    with open(inpf) as f:
        for line in f:
            if len(line) > 4:
                if line[0] == "@":
                    if "AVAM" in line:
                        if line.split()[1] == "=":
                            if line.split()[2] == "on":
                                tfavam = 1
                            elif line.split()[2] == "off":
                                tfavam = 0
                            else:
                                stop("Error: Please specify AVAM correctly in the input file.")

    f.close()

    cmdgui.TFavam = tfavam
    #print "avam", tfavam

    # DIRCURV
    with open(inpf) as f:
        for line in f:
            if len(line) > 4:
                if line[0] == "@":
                    if "DIRCURV" in line:
                        if line.split()[1] == "=":
                            if line.split()[2] == "on":
                                tfdircurv = 1
                            elif line.split()[2] == "off":
                                tfdircurv = 0
                            else:
                                stop("Error: Please specify DIRCURV correctly in the input file.")

    f.close()

    cmdgui.TFdircurv = tfdircurv
    #print "dircurv",tfdircurv




    # PARM
    with open(inpf) as f:
        for line in f:
            if len(line) > 4:
                if line[0] == "@":
                    if "PARM" in line:
                        if line.split()[1] == "=":
                            tmpans = line.split()[2]
                            if tmpans.upper() == "NO":
                                tfparm = 0
                            elif tmpans.upper() == "GEOMONLY":
                                tfparm = 1
                            elif tmpans.upper() == "ALL":
                                tfparm = 2
                            else:
                                stop("Error: Please specify PARM correctly in the input file.")
    f.close()

    cmdgui.TFparm = tfparm # <-
    #print "parm",tfparm


    # VIBRATION
    with open(inpf) as f:
        for line in f:
            if len(line) > 8:
                if line[0] == "@":
                    if "VIBRATION" in line:
                        if line.split()[1] == "=":
                            if line.split()[2] == "on":
                                tfvib = 1
                            elif line.split()[2] == "off":
                                tfvib = 0
                            else:
                                stop("Error: Please specify VIBRATION correctly in the input file.")

                            cmdgui.TFvib = tfvib # <-
    f.close()
    #print "vib",tfvib





    # DATAFILETYPE

    dftype = 'no'
    with open(inpf) as f:
        for line in f:
            if len(line) > 4:
                if line[0] == "@":
                    if "DATAFILETYPE" in line:
                        if line.split()[1] == "=":
                            dftype = line.split()[2]
    f.close()
    if dftype == "old":
        oldnew = 0
    elif dftype == "new":
        oldnew = 1
    elif dftype == "xyz":
        oldnew = 2
    else:
        stop("Error: Please specify DATAFILETYPE correctly in the input file.")


    #if oldnew == -1:
    #   stop("Error: Please specify DATAFILETYPE correctly in the input file.")
    cmdgui.Datfiletype = oldnew  # <-



    # BASEPATH

    with open(inpf) as f:
        for line in f:
            if len(line) > 4:
                if line[0] == "@":
                    if "BASEPATH" in line:
                        if line.split()[1] == "=":
                            tmpstr = line.split()[2]
                            strn = 0
                            if tmpstr[0] != "\"" and tmpstr[0] != "\'":
                                strn = 1
                            if tmpstr[-1] != "\"" and tmpstr[-1] != "\'":
                                strn = 1
                            if strn == 1:
                                stop("Error: Please specify BASEPATH correctly in the input file.")
                            bpath = line.split()[2]
                            if len(bpath) < 3:
                                stop("Error: BASEPATH specification is empty.")
                            cmdgui.Basename = bpath[1:-1]


    f.close()



    # DATAFILEPATH

    with open(inpf) as f:
        for line in f:
            if len(line) > 4:
                if line[0] == "@":
                    if "DATAFILEPATH" in line:
                        if line.split()[1] == "=":
                            tmpstr = line.split()[2]
                            strn = 0
                            #print tmpstr[0]
                            #print tmpstr[-1]
                            if tmpstr[0] != "\"" and tmpstr[0] != "\'":
                                strn = 1
                            if tmpstr[-1]!= "\"" and tmpstr[-1]!= "\'":
                                strn = 1
                            if strn == 1:
                                stop("Error: Please specify DATAFILEPATH correctly in the input file.")
                            dfilepath = line.split()[2]
                            if len(dfilepath) < 3:
                                stop ("Error: DATAFILEPATH specification is empty.")

                            nBrowsf = dfilepath
                            oBrowsf = dfilepath
                            xyzf    = dfilepath
                            cmdgui.Datfilename = dfilepath[1:-1] # <-

    f.close()


    # CURVCOR interface input
    cv = 0
    with open(inpf) as f:
        for line in f:
            if len(line) > 3:
                if line[0:7].upper() == "CURVCOR":
                    cv = 1
                    cmdgui.TFcurvcorts = 1 # <-
                    continue
                if cv == 1:
                    if "Ln" in line:
                        if line.split()[1] != "=":
                            stop("Error: Ln specification not correct in CURVCOR section.")
                        curvcorln = int(line.split()[2])
                        cmdgui.NcurvcorLn = curvcorln # <-
                        #print curvcorln
                        continue
                    if "Rn" in line:
                        if line.split()[1] != "=":
                            stop("Error: Rn specification not correct in CURVCOR section.")
                        curvcorrn = int(line.split()[2])
                        cmdgui.NcurvcorRn = curvcorrn # <-
                        #print curvcorrn
                        continue
                    if "END" in line:
                        cv = -1
                        continue

                if cv == -1:
                    break
    f.close()

    # AUTOSMOOTH interface input

    cv = 0
    istep = 0
    with open(inpf) as f:
        for line in f:
            if len(line) > 3:
                if line[0:8].upper() == "AUTOSMTH":
                    cv = 1
                    cmdgui.TFautosmth = 1 # <-
                    continue
                if cv == 1:
                    if "StepSize" in line:
                        if line.split()[1] != "=":
                            stop("Error: StepSize specification not correct in AUTOSMTH section.")
                        autosmthstepsize = float(line.split()[2])
                        if autosmthstepsize <= 0.0 or autosmthstepsize >=1.0:
                            stop("Error: StepSize value not valid in AUTOSMTH section.  ")
                        istep = 1
                        #print autosmthstepsize
                        cmdgui.Nautosmthstpsize = autosmthstepsize # <-
                        continue
                    if "Ln" in line:
                        if line.split()[1] != "=":
                            stop("Error: Ln specification not correct in AUTOSMTH section.")
                        autosmthln = int(line.split()[2])
                        if autosmthln <= 0:
                            stop("Error: Ln value not valid in AUTOSMTH section.")
                        if autosmthln > 10:
                            print ("Warning: Ln value large than normal setting in AUTOSMTH section.  ")
                        #print autosmthln
                        cmdgui.NautosmthLn = autosmthln # <-
                        continue
                    if "Rn" in line:
                        if line.split()[1] != "=":
                            stop("Error: Rn specification not correct in AUTOSMTH section.")
                        autosmthrn = int(line.split()[2])
                        if autosmthrn <= 0:
                            stop("Error: Rn value not valid in AUTOSMTH section.")
                        if autosmthrn > 10:
                            print ("Warning: Rn value large than normal setting in AUTOSMTH section.  ")
                        #print autosmthrn
                        cmdgui.NautosmthRn = autosmthrn # <-
                        continue
                    if "d2ythresh" in line:
                        if line.split()[1] != "=":
                            stop("Error: d2ythresh specification not correct in AUTOSMTH section.")
                        autosmthd2y = float(line.split()[2])
                        if autosmthd2y < 0:
                            stop("Error: d2ythresh value not valid in AUTOSMTH section.")
                        #print autosmthd2y
                        cmdgui.Nautosmthd2y = autosmthd2y # <-
                        continue

                    if "END" in line:
                        cv = -1
                        if istep == 0:
                            stop("Error: StepSize specification missing in AUTOSMTH section.")
                        continue

                if cv == -1:
                    break
    f.close()



    # RMSPK interface input

    cv = 0
    with open(inpf) as f:
        for line in f:
            if len(line) > 3:
                if line[0:5].upper() == "RMSPK":
                    cv = 1
                    cmdgui.TFremovespk = 1 # <-
                    continue
                if cv == 1:
                    if "CutHigh" in line:
                        if line.split()[1] != "=":
                            stop("Error: CutHigh specification not correct in RMSPK section.")
                        rmspkcuthigh = float( line.split()[2] )
                        #print rmspkcuthigh
                        cmdgui.Nrmspkcut = rmspkcuthigh # <-
                        continue
                    if "Percentage" in line:
                        if line.split()[1] != "=":
                            stop("Error: Percentage specification not correct in RMSPK section. ")
                        rmspkpercent = float( line.split()[2] )
                        if rmspkpercent < 0.0 or rmspkpercent > 1.0:
                            stop("Error: Percentage value not valid in RMSPK section.")
                        if rmspkpercent < 0.6:
                            print ("Warning: Percentage value smaller than normal setting in RMSPK section. ")
                        #print rmspkpercent
                        cmdgui.Nrmspkperc = rmspkpercent # <-
                        continue
                    if "GradRatio" in line:
                        if line.split()[1] != "=":
                            stop("Error: GradRatio specification not correct in RMSPK section.")
                        rmspkgradratio = float( line.split()[2] )
                        if rmspkgradratio < 1.0:
                            stop("Error: GradRatio value should be larger than 1.0 .")
                        #print rmspkgradratio
                        cmdgui.Nrmspkgratio = rmspkgradratio # <-
                        continue

                    if "END" in line:
                        cv = -1
                        continue

                if cv == -1:
                    break

    f.close()



    # Print out Title Section info.
    cv = 0
    titlelines = []
    #print "a"
    with open(inpf) as f:
        for line in f:
            if len(line) > 1:
                if cv == -1:
                    break

                if cv == 1:
                    titlelines.append( line )

                if line[0:5].upper() == "TITLE":
                    cv = 1
                    continue

                if len(line)>=9 and line[0:9].upper() == "END TITLE":
                    cv = -1
                    continue

    f.close()
    #print titlelines
    del titlelines[-1]
    #print titlelines
    maxlen = 0
    for mm in titlelines:
        if len(mm) > maxlen:
            maxlen = len(mm)

    #print maxlen

    for mm in range(maxlen):
        sys.stdout.write("-")
    sys.stdout.write("\n")

    for mm in range(len(titlelines)):
        sys.stdout.write(titlelines[mm])

    for mm in range(maxlen):
        sys.stdout.write("-")
    sys.stdout.write("\n")



    # Prepare the return list
    inplist = []


    # Return user input list
    inplist.append( cmdgui.Datfiletype ) # 0

    if cmdgui.Datfilename == "InputDataFile":
        stop("Error: DATAFILEPATH specification is missing...")
    inplist.append( cmdgui.Datfilename ) # 1


    inplist.append( cmdgui.TFparm ) # 2
    inplist.append( cmdgui.TFvib ) # 3

    inplist.append( cmdgui.TFcurvcorts ) # 4
    inplist.append( cmdgui.NcurvcorLn ) # 5
    inplist.append( cmdgui.NcurvcorRn ) # 6

    inplist.append( cmdgui.TFautosmth ) # 7
    inplist.append( cmdgui.NautosmthLn ) # 8
    inplist.append( cmdgui.NautosmthRn ) # 9
    inplist.append( cmdgui.Nautosmthstpsize ) # 10
    inplist.append( cmdgui.Nautosmthd2y ) # 11

    inplist.append( cmdgui.TFremovespk ) # 12
    inplist.append( cmdgui.Nrmspkcut ) # 13
    inplist.append( cmdgui.Nrmspkperc ) # 14
    inplist.append( cmdgui.Nrmspkgratio ) # 15

    inplist.append( cmdgui.TFdircurv ) # 16
    inplist.append( cmdgui.TFavam ) # 17
    inplist.append( cmdgui.TFcurvcpl ) # 18
    inplist.append( cmdgui.TFcoriolis ) # 19

    inplist.append( cmdgui.TFenergy ) # 20
    inplist.append( cmdgui.TFadiabfc ) # 21

    if cmdgui.Basename == "InputBaseName":
        stop("Error: BASEPATH specification is missing...")

    if cmdgui.Basename[-1] != "/":
        cmdgui.Basename = cmdgui.Basename + "/"

    # Test whether the folder exist
    if not os.path.isdir(cmdgui.Basename):
        stop("Error: folder not exist- "+cmdgui.Basename)
    if not os.path.exists(cmdgui.Basename+"main.py"):
        stop("Error: could not find pURVA in "+cmdgui.Basename)
    if not os.path.exists(cmdgui.Basename+"inter_removespikes_main.py"):
        stop("Error: RMSPK module missing in "+cmdgui.Basename)
    if not os.path.exists(cmdgui.Basename+"inter_autosmooth_main.py"):
        stop("Error: AUTOSMTH module missing in "+cmdgui.Basename)
    if not os.path.exists(cmdgui.Basename+"lm90.test.exe"):
        stop("Error: ADIABLM module missing in "+cmdgui.Basename)



    inplist.append( cmdgui.Basename ) # 22


    #


    return inplist

    #



def printqncomp(compdata, parID, filename,io):

    if os.path.isfile(filename):
        stop("Error: this file already exists- "+filename)
    fh = open(filename,"w")

    nparm = len(parID)

    names = ["s"]
    for i in parID:
        names.append( i[-1] )   # As the name of each parameter is specified
    #print names

    fh.write( "=============================================================\n")
    if io == 1: # eta
        fh.write("| DECOMPOSITION OF PATH DIRECTION INTO INTERNAL COORDINATES |\n")
    else:    # kappa
        fh.write("| DECOMPOSITION OF PATH CURVATURE INTO INTERNAL COORDINATES |\n")
    fh.write( "=============================================================\n")

    for i in range(len(names) -1):
        fh.write( names[i] )
        fh.write( ',' )
    fh.write( names[-1] )
    fh.write( '\n' )

    for i in range(len(compdata)):
        for j in range( len(compdata[i]) - 1 ):
            fh.write( str(compdata[i][j]) )
            fh.write( "," )
        fh.write( str(compdata[i][-1]) )
        fh.write( "\n" )

    if io == 1:
        print(("Output: Direction components written to "+filename))
    else:
        print(("Output: Curvature components written to "+filename))



# A general function to read in user input file with section
def rdkeywrds( inpfile, secname, keyname ):
    keyname = keyname.strip()
    secname = secname.strip()

    val = 'toberead'
    istart = 0
    iend   = 0
    texts = []

    with open(inpfile) as f:
        for line in f:
            if len(line) > 2:
                if secname in line:
                    if line.split()[0] == secname:
                        istart = 1
                        iend = -1
                        continue
                if istart == 1:
                    if line.split()[0] == "END" and line.split()[1] == secname:
                        iend = 1
                        break
                    else:
                        texts.append( line )
                        continue

    if istart == 0:
        print(("Warning: section "+secname+" not found in "+inpfile))
        val = 'toberead'

    if iend == -1:
        stop("Error: No ending line for the section of "+secname)


    label = 0
    for i in range(len(texts)):
        if len(texts[i]) > 1:
            if texts[i].split()[0] == keyname:
                if texts[i].split()[1] == "=":
                    sval = texts[i].split()[2]

                    if "E" in sval or "e" in sval or "." in sval:
                        val = float(sval)
                    else:
                        val = int(sval)
                    label = 1
    if label == 0:
        print(("Warning: "+keyname+" not found in the section of "+secname))

    #print val
    return val



#rdkeywrds("inp.inp","DMO","Slowest")





#PrtEnergy(nBlist, tag, "energy.csv", "energy_1_d.csv", "energy_2_d.csv")
def PrtEnergy(Brlist, energyfile, d1file, d2file):


    sval = []
    eval = []

    for i in range(len(Brlist)):
        sval.append( Brlist[i].s )
        eval.append( Brlist[i].energy )



    if os.path.isfile(energyfile):
        stop("Error: this file already exists- "+energyfile)

    f1h = open(energyfile, "w")
    f1h.write("   s  energy(a.u.)\n")


    for i in range(len(Brlist)):
        f1h.write("%8.3f %16.9f  %s" % (sval[i],eval[i],"\n")    )

    print(("Output: energy vs. s written to file- "+energyfile))



    # tag = 2 -> calculate derivatives


    cs = CubicSpline( sval, eval )
    x0 = sval[ 0]
    x1 = sval[-1]
    xs = np.arange(x0,x1,0.0005)

    ys = cs( xs )

#    print xs
#    print ys


#    for i in range(len(xs)):
#       f1h.write("%10.4f %16.9f  %s" % (xs[i],ys[i],"\n")    )


    if os.path.isfile(d1file):
        stop("Error: this file already exists- "+d1file)
    f2h = open( d1file,"w" )
    f2h.write("   s d_energy/d_s(a.u.)\n")


    d1x = []
    d1y = []

    for i in range(len(xs)):
        if i > 0 and i < (len(xs)-1):
            d1x.append( xs[i] )
            tmp = (ys[i+1]-ys[i-1])/(xs[i+1]-xs[i-1])
            d1y.append( tmp )
            #print xs[i], tmp
            f2h.write("%10.4f %16.9E %s" % (xs[i], tmp, "\n"   )   )

    # print out d1e to file
    print(("Output: d_energy/d_s vs. s written to file- "+d1file))


    if os.path.isfile(d2file):
        stop("Error: this file already exists- "+d2file)
    f3h = open(d2file,"w")


    d2x = []
    d2y = []

    for i in range(len(xs)):
        if i > 0 and i < (len(d1x)-1):
            d2x.append( d1x[i] )
            tmp = (d1y[i+1]-d1y[i-1])/(d1x[i+1]-d1x[i-1])
            d2y.append( tmp )
            #print d1x[i],tmp



    # smoothen -0.1 ~ 0.1 region with cubic spline fitting
    sL = -0.1
    sR =  0.1
    xtrain = []
    ytrain = []

    xpred = []
    ypred = []
    for i in range(len(d2x)):
        if d2x[i] < sL or d2x[i] > sR:
            xtrain.append( d2x[i] )
            ytrain.append( d2y[i] )
        else:
            xpred.append( d2x[i] )
            ypred.append( d2y[i] )

    skip = 0
    if len(xpred) == 0:
        skip = 1



    if skip == 0:
        cs2 = CubicSpline( xtrain, ytrain )
        xs = np.arange( xpred[0], xpred[-1], 0.0005 )
        #ys = cs2(xs)

        #d2x = []
        #d2y = []
        for i in range(len(d2x)):
            if d2x[i] < sL or d2x[i] > sR:
                pass
            else:
                d2y[i] = cs2( d2x[i] )



    f3h.write("   s d^2_energy/d_s^2(a.u.)\n")

    for i in range(len(d2x)):
        f3h.write( "%10.4f %16.9E %s" % (d2x[i],d2y[i],"\n") )
        #print d2x[i],d2y[i]

    print(("Output: d^2_energy/d_s^2 vs. s written to file- "+d2file))

    pass





# checkduply( Brlist ) - 05/17/2017
def checkduply( Brlist ):
    # Check duplicate points on reaction path

    deletelist = []
    newBrlist = []

    numPn = len( Brlist )

    for i in range(numPn-1):
        j = i + 1

        s1 = Brlist[i].s
        s2 = Brlist[j].s

        if abs(s1-s2)<1e-5:
            print(("Warning: duplicate s value found for "+str(s1)))

            c1 = copy.deepcopy( Brlist[i].Vcoor )
            c2 = copy.deepcopy( Brlist[j].Vcoor )
            c1 = [item for sublist in c1 for item in sublist]
            c2 = [item for sublist in c2 for item in sublist]

            c1 = np.array( c1 )
            c2 = np.array( c2 )
            #print c1
            diff = c1 - c2
            dif = np.linalg.norm( diff )
            if dif < 1e-5:
                print("         Info:  duplicate points identified")
                deletelist.append( j )
            else:
                stop("          Error: Cartesian coordinates not match...")
            #print("       norm = "+str(dif))






    if len(deletelist) ==  0:
        print("Info: no duplicate points found.")
        newBrlist = copy.deepcopy( Brlist)

    else:
        for i in range(numPn):
            if i not in deletelist:
                newBrlist.append( Brlist[i] )
        newNp = len(newBrlist)
        print(("Info: Total number of points in the data file: "+str(newNp)))
        pass





    #print len(newBrlist)

    pass
    return newBrlist
