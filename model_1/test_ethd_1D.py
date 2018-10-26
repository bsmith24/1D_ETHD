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
import test_ho
import models

class tmp:
    pass    

def run_QC_1D(params, case, ent_opt):

    """
    Runs tests for 1D model problems using ensembles of trajectories
    """
    
    ndia, nadi, nnucl, ntraj = 1, 1, 1, 10000

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
    
    rnd = Random()
    q = MATRIX(nnucl,ntraj)
    p = MATRIX(nnucl,ntraj)

    if ent_opt == 0:

        for n in xrange(len(params["wfc"]["weights"])):
    
            params["wfc"]["weights"][n] = params["wfc"]["weights"][n].real

        sampling = test_ho.metropolis_gau(rnd, test_ho.HO_sup, q, params, ntraj, 2500, 0.05)
        print len(sampling)
        test_ho.bin(sampling, -1.5, 2.0, 0.01, 0, 0, "_distrib-"+str(case)+".txt") 

        for i in xrange(q.num_of_rows):
            for j in xrange(q.num_of_cols):
                q.set(i,j,sampling[j].get(i,0))
                p.set(i,j,params["p0"])    

        aux_functs.save_q_p_info(q,p,params,ent_opt,case) 

    else:
        # Now we are extracting the initial coordiates from the previous (classical) simulation 
    
        q, p = aux_functs.get_q_p_info(params, 0, case)

    # Set mass matricies  
    iM = MATRIX(nnucl,1);
    iM.set(0,0, 1.0/params["mass"])
   
    # Compute Hamiltonian properties
    ham.compute_diabatic(models.compute_model, q, params, 1)
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
            Verlet1(dt, q, p, iM, ham, models.compute_model, params, ent_opt)
