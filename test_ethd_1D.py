#*********************************************************************************
#* Copyright (C) 2017-2018  Brendan A. Smith, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/


import cmath
import math
import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *
import aux_functs

class tmp:
    pass    


def model1(q, params, full_id):
    """
    Hdia = 0.5*k*x^2   
    Sdia =  1.0
    Ddia  = 0.0
    """

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(1,1) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(1,1) )
  
    x = q.col(indx).get(0)
    x0,k = params["x0"], params["k"]

    Hdia.set(0,0, k*x*x*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) )

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j)


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj
    
def model2(q, params, full_id):
    """
    Morse Potential
    Hdia = D * ( 1 - exp ( alp(x-x0) ) )^2   
    Sdia = 1.0
    Ddia = 0.0
    """

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(1,1) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(1,1) )
 
    x = q.col(indx).get(0)
    k, D, x0 = params["k"], params["D"], params["x0"]

    # Calculate alpha
    alp = math.sqrt(0.5*k/D)    
    # Compute potential
    V = ( 1.0 - math.exp( -alp*(x-x0) ) )
    
    Hdia.set(0,0, D*V*V*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)
    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set( 0,0, ( 2.0*D*V*alp*math.exp(-alp*(x-x0)) )*(1.0+0.0j) )

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j)


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


def model3(q, params, full_id):
    """
    Sin Potential
    Hdia = A * sin( k*(q-B) ) + C    
    Sdia = 1.0
    Ddia = 0.0
    """

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(1,1) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(1,1) )

    x = q.col(indx).get(0)

    A = params["A"]
    B = params["B"]
    C = params["C"]
    k = params["k"]

    Hdia.set(0,0, (A*math.sin(k*(x-B))+C)*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, ( A*k*math.cos(k*(x-B)) )*(1.0+0.0j) )

        #  <dia| d/dR_0| dia>
        dc1_dia[i].set(0,0, 0.0+0.0j)


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


def model4(q, params, full_id):
    """
    1D Eckart Barrier Potential

    Define system specific parameters 
    Va = barrier height    
    Sdia = 1.0
    Ddia = 0.0
    """

    # Define potential specific constants
    Va = params["Va"]

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(1,1) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(1,1) )

    x = q.col(indx).get(0)
    x2 = x*x

    # z = sech(2x)
    # z2 = sech^2(2x)    
    z = (2.0 * math.cosh(2.0*x))/(1.0 + math.cosh(4.0*x))
    z2 = z*z

    Hdia.set(0,0, Va*z2*(1.0+0.0j))
    Sdia.set(0,0, 1.0+0.0j)
    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, ( ( -4.0*Va*math.tanh(2.0*x)*z2 ) )*(1.0+0.0j))

        #  <dia| d/dR_0| dia>
        dc1_dia[i].set(0,0, 0.0+0.0j)

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


def compute_model(q, params, full_id):
    """
    Generic calculator of the model Hamiltonians
    """

    model = params["model"]
    res = None

    if model==1:
        # HO
        res = model1(q, params, full_id)
    if model==2:
        # DW
        res = model2(q, params, full_id)
    if model==3:
        # 1D-sin 
        res = model3(q, params, full_id)
        #res = models_Libra.model7(q, params)
    if model==4:
        # 1D Eckart Barrier
        res = model4(q, params, full_id)

    return res


def run_QC_1D(params, case, ent_opt):

    """
    Runs tests for 1D model problems using ensembles of trajectories
    """
    
    ndia, nadi, nnucl, ntraj = 1, 1, 1, 5000

    # ======= Hierarchy of Hamiltonians =======
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.init_all(2)
    print "id=", ham.id, " level=", ham.level

    ham1 = []
    for tr in xrange(ntraj):
        ham1.append( nHamiltonian(ndia, nadi, nnucl) )
        #print ham1[tr].id, ham1[tr].level
        ham1[tr].init_all(2)
        ham.add_child(ham1[tr])
        #print Cpp2Py(ham1[tr].get_full_id())

    # Initialize Simulation 
    dt = params["dt"]
    nsnaps, nsteps = params["nsnaps"], params["nsteps"] 
    params.update({"ndof":nnucl, "ntraj":ntraj})
 
    print "ent_opt = ", ent_opt
    print "model = ", params["model"]
    print "case = ", case
   
    # Dynamical variables and system-specific properties
    mean_q = MATRIX(nnucl,1);   
    sigma_q = MATRIX(nnucl,1);
    mean_p = MATRIX(nnucl,1);   
    sigma_p = MATRIX(nnucl,1);  

    mean_q.set(0,0,params["q0"])  
    sigma_q.set(0,0,params["sq0"])
    mean_p.set(0,0,params["p0"])

    if ent_opt == 0:
        sigma_p.set(0,0,params["sp0"])
    elif ent_opt == 1:
        sigma_p.set(0,0,0.0)

    rnd = Random()
    q = MATRIX(nnucl,ntraj);  tsh.sample(q, mean_q, sigma_q, rnd)        
    p = MATRIX(nnucl,ntraj);  tsh.sample(p, mean_p, sigma_p, rnd)   

    # Set mass matricies  
    iM = MATRIX(nnucl,1);
    iM.set(0,0, 1.0/params["mass"])
   
    # Compute Hamiltonian properties
    ham.compute_diabatic(compute_model, q, params, 1)
    ham.compute_adiabatic(1, 1);
    if ent_opt == 1:
        ham.add_ethd_adi(q, iM, 1)

    os.system("mkdir _1D_dist_qc")
    out1 = open("_output_"+str(params["model"])+str(ent_opt)+str(case)+".txt", "w"); out1.close()   
    # Do the propagation
    for i in xrange(nsnaps):

        aux_functs.bin(q, -50.0, 50.0, 0.1, "_1D_dist_qc/_dist_"+str(params["model"])+str(ent_opt)+str(case)+"_"+str(i*int(dt)*nsteps)+".txt")
        expect_q, expect_p, expect_q2, expect_p2 = aux_functs.compute_properties(q,p)
        sx2  = expect_q2  - expect_q*expect_q
        spx2 = expect_p2  - expect_p*expect_p
        uncp = sx2 * spx2

        # Count the number of trajectories that cross the barrier 
        react_prob = aux_functs.traj_counter(q, params)

        #=========== Properties ==========
        Ekin, Epot, Etot = aux_functs.compute_etot(ham, p, iM)

        # Print the ensemble average - kinetic, potential, and total energies
        # Print the tunneling information. Here, we count each trajectory across the barrier.
        out1 = open("_output_"+str(params["model"])+str(ent_opt)+str(case)+".txt", "a")
        out1.write( "%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" % (i*dt*nsteps, Ekin, Epot, Etot, react_prob, expect_q, expect_p, expect_q2, expect_p2, sx2, spx2, uncp) )
        out1.close()
        
        for j in xrange(nsteps):
            Verlet1(dt, q, p, iM, ham, compute_model, params, ent_opt)



def make_fig(models,ent_opts,cases):

    # Need to specific what exactly the "cases" are. 
    data = MATRIXList()
    for i in xrange(len(models)):

        data.append( MATRIX(len(ent_opts),len(cases)) )

        for j in xrange(len(ent_opts)):
            for k in xrange(len(cases)):

                f = open("_output_"+str(models[i])+str(ent_opts[j])+str(cases[k])+".txt", "r")
                l = f.readlines()
                f.close()
                b = l[-1].strip().split()
                data[i].set(j,k,float(b[4]))
      
        out1 = open("_fig_model"+str(models[i])+".txt", "w"); out1.close()
        for j in xrange(len(cases)):
            out1 = open("_fig_model"+str(models[i])+".txt", "a")
            if len(ent_opts) == 1:
                out1.write( "%8.5f %8.5f\n" % ( cases[j] + 0.0, data[i].get(0,j) ) )
                out1.close()
            elif len(ent_opts) == 2:
                out1.write( "%8.5f %8.5f %8.5f\n" % ( cases[j] + 0.0, data[i].get(0,j), data[i].get(1,j) ) )
                out1.close()


