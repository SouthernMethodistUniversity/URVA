import numpy as np
import sys



def directcut(y, label, thresh):
    num = len(label)
    for i in range(num):
        if y[i] > thresh:
            if label[i] == 0:
                label[i] = 1
    return label







def gradientfilter(x,y,label,thresh):
    # check both left and right part
    ################################  if both sides fail, then this is a spike point
    # if either side passes, change this possible spike point into a normal point.

    #print len(x), len(y), len(label)
    num = len(label)

    for i in range(num):
        if i <= 2 or i >= (num-1-2):
            continue
            pass
        else: # we do the gradient check
            if label[i] == 1:
                flag = 0
                gl0 = (y[i  ] - y[i-1])/(x[i  ] - x[i-1])
                gl1 = (y[i-1] - y[i-2])/(x[i-1] - x[i-2])
                gr0 = (y[i+1] - y[i  ])/(x[i+1] - x[i  ])
                gr1 = (y[i+2] - y[i+1])/(x[i+2] - x[i+1])
                if gl0 == 0.0 or gr0 ==0.0:
                    continue
                if (gl1/gl0) >= (1.0/thresh) and (gl1/gl0) <= thresh:
                    #print "Get one point",x[i]
                    flag = 1
                if (gr1/gr0) >= (1.0/thresh) and (gr1/gr0) <= thresh:
                    #print "Get one point",x[i]
                    flag = 1
                if flag == 1:
                    label[i] = 0
    return label



def checkpercent( x,y,diff,thresh ):

    num = len(y)
    label = list(range(num))
    for i in range(num):
        label[i] = 0

    if thresh < 1.0:
        thresh = 100 * thresh

    tmp = np.percentile(  np.array(diff) , thresh )
    #print tmp
    #print len(y), len(diff)

    for i in range(len(y)):
        #if diff[i] < tmp:
        #   print x[i],"     ",y[i]
        if diff[i] > tmp:
            label[i] = 1
    #
    #print sum(label)
    return label



def calcdiff(y):
    diff=[]
    num = len(y)
    #print "Number is num",num
    for i in range(1, (num-1)):
        x = y[i]
        xn =y[i-1]
        xp =y[i+1]
        tmp = ( abs(x-xn) + abs(x-xp) )/2
        diff.append( tmp )

    tmp1 = abs( y[0] - y[1] )
    diff.insert( 0, tmp1)
    tmp2 = abs( y[-1]- y[-2])
    diff.append( tmp2 )

    return diff






#a=[1,2,3,4]
#calcdiff(a)



def checkdupli(s):
    num = len(s)
    k = 0
    todelete=[]
    for i in range(num - 1):
        if s[i] == s[i+1]:
            print(("Warning: duplicate s found: "+str(s[i])))
            k = 1
            todelete.append( i+1 )
    if k == 1:
        pass
        #sys.exit()
    else:
        pass
        #print ("No duplicate points found. Proceeding...")

    return todelete


def readdat(filename):
    x=[]
    y=[]
    with open(filename) as f:
        for line in f:
            #print line
            if len(line) > 2:
                x.append( float( line.split()[0] ) )
                y.append( float( line.split()[1] ) )
    #print x
    #print y
    return (x,y)
