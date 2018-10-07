import sys
import math

nsnaps = 100
nsteps = 50
dt = 1.0
dx = 0.01
iconds = 1
k = 1.0

barrier =  2.0*math.pi/k

for k in xrange(iconds):

    c = open("_qprob"+str(k)+".txt","w")

    for i in xrange(nsnaps):
        sum_all = 0.0
        prob_all = 0.0
        sum_cross = 0.0
        prob_cross = 0.0
        a = open("res"+str(k)+"/wfc.state0.frame"+str(i),"r")

        for line in a:

            b = line.strip().split()
            if len(b) < 1:
                continue
            else:
                if   float(b[0]) > barrier:    
                    sum_cross += float(b[1])*dx
                elif float(b[0]) < -barrier:    
                    sum_cross += float(b[1])*dx


        c.write( str(i*dt*nsteps) + " " +  str(sum_cross) + "\n" )




