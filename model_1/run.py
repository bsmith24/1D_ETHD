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


import harmonic
import test_ethd_1D
import models


class tmp:
    pass


def potential(q, params):
    """
    Thin wrapper of the model Hamiltonians that can be used in 
    the fully-quantum calculations
    """

    return models.compute_model(q, params, Py2Cpp_int([0,0]))
   



def plot_pes_1D(params,case):
    """
    An auxiliary function to compute the PES profiles
    """

    x0 = -50.0
    q = MATRIX(1,1)

    f = open("_pes"+str(case)+".txt", "w")
    f.close()

    for i in xrange(10000):
        x = x0 + 0.01*i
        q.set(0,0,x)
        full_id = Py2Cpp_int([0,0])
        obj = models.compute_model(q, params, full_id)

        f = open("_pes"+str(case)+".txt", "a")
        f.write( "%8.5f   %8.5f   %8.5f \n" % (x, obj.ham_dia.get(0,0).real, obj.d1ham_dia[0].get(0,0).real))
        f.close()


def run_exact(params, case):
    """
    The main routine to run fully quantum calculations
    """

    dt = params["dt"]
    # Here we initialize the grid and wavefunction
    nsnaps, nsteps = params["nsnaps"], params["nsteps"]

    # 1-D Initilization
    wfc = Wfcgrid(-50.0, 50.0, 0.01, 1)
    wfc.init_wfc_1D(params["q0"], params["p0"], params["sq0"], 0)
    wfc.update_potential_1D(potential, params)
    wfc.update_propagator_1D(0.5*dt, params["mass"])  # this is important because we are using exp(-0.5*dt*H_loc)...
    wfc.update_propagator_K_1D(dt, params["mass"])    # ... together with exp(-dt*H_non-loc)    

    wfc.normalize_wfc_1D()

    f = open("_pops"+str(case)+".txt", "w")
    f.close()
  

    exp_pow = doubleList()
    exp_pow.append(2.0)

    flux = doubleList()
    flux.append(0.0)          # len(flux[i]) = nstates
    tot_flux = 0.0
    cum = 0.0
    os.system("mkdir _res"+str(case)+"")
    for i in xrange(nsnaps):  # time steps

        wfc.print_wfc_1D("_res"+str(case)+"/wfc", int(i*dt*nsteps), 0)   # Be sure to make a "res" directory
        #wfc.print_reci_wfc_1D("res/reci_wfc", i*dt*nsteps, 0)

        epot = wfc.e_pot_1D()
        ekin = wfc.e_kin_1D(params["mass"])
        etot = epot + ekin
        x, px   = wfc.get_x_1D(),  wfc.get_px_1D()
        x2, px2  = wfc.get_pow_x_1D(exp_pow[0]), wfc.get_pow_px_1D(2)

        # Compute the uncertainties and uncertainy products
        sx2  = x2  - x*x
        spx2 = px2 - px*px 
        uncp = sx2 * spx2

        f = open("_pops"+str(case)+".txt", "a")
        f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" % (i*nsteps*dt, ekin, epot, etot, cum, x, px, x2, px2, sx2, spx2, uncp, wfc.norm_1D())) 
        f.close()


        for j in xrange(nsteps):  # time steps 
            wfc.propagate_exact_1D(0)

            #res0 = Py2Cpp_double([0.0])
            #wfc.flux_1D( params["barrier"], res0, params["mass"])
            #cum += res0[0]*dt 

            #res1 = Py2Cpp_double([0.0])
            #wfc.flux_1D(-params["barrier"], res1, params["mass"])
            #cum -= res1[0]*dt


def init_params(model, case):
    """
    This function intiializes the params dictionary for a given case    
    """

    params = {"dt":1.0}
 

    if model == 0:
        # Constant Potential, V(q) = 0.0
        params.update( {"model":0} )
        params.update( {"barrier":50} )

        if  case == 0:
            params.update( {"q0":0.0, "p0":0.0, "sq0":0.25, "mass":2000.0} )


    elif model == 1:
        # Harmonic Potential
        params.update( {"mass":2000.0, "barrier":50.00, "k":0.2} )

        if case == 0:

            params.update( {"q0":0.0, "p0":0.0} )

            # params needed for numerically exact simulation
            wfc = {}
            wfc.update({"init_state":[0], "nu":[0], "weights":[1.0+0.0j], "x0":[params["q0"]], "px0":[params["p0"]]})
            alp = math.sqrt(params["k"] * params["mass"])
            wfc.update({"alpha":[ alp ] })
            params.update( {"model":1, "wfc": wfc} )

        elif case == 1:
            params.update( {"q0":0.0, "p0":0.0} )
            wfc = {}
            wfc.update({"init_state":[0,0], "nu":[0,1], "weights":[1.0+0.0j,1.0+0.0j], "x0":[params["q0"],params["q0"]], "px0":[params["p0"],params["p0"]]})
            alp = math.sqrt(params["k"] * params["mass"])
            wfc.update({"alpha":[ alp,alp ] })
            params.update( {"model":1, "wfc": wfc} )

        elif case == 2:
            params.update( {"q0":0.0, "p0":0.0} )
            wfc = {}
            wfc.update({"init_state":[0,0,0], "nu":[0,1,2], "weights":[1.0+0.0j,1.0+0.0j,1.0+0.0j], "x0":[params["q0"],params["q0"],params["q0"]], "px0":[params["p0"],params["p0"],params["p0"]]})
            alp = math.sqrt(params["k"] * params["mass"])
            wfc.update({"alpha":[ alp,alp,alp ] })
            params.update( {"model":1, "wfc": wfc} )

        elif case == 3:
            params.update( {"q0":0.0, "p0":0.0} )
            wfc = {}
            wfc.update({"init_state":[0,0,0,0], "nu":[0,1,2,3], "weights":[1.0+0.0j,1.0+0.0j,1.0+0.0j,1.0+0.0j], "x0":[params["q0"],params["q0"],params["q0"],params["q0"]], "px0":[params["p0"],params["p0"],params["p0"],params["p0"]]})
            alp = math.sqrt(params["k"] * params["mass"])
            wfc.update({"alpha":[ alp,alp,alp,alp ] })
            params.update( {"model":1, "wfc": wfc} )



        # Compute omega based on k
        omega = math.sqrt(params["k"]/params["mass"])
        params.update({"omega":omega})

    elif model==2:
        # 1D morse oscillator
        params.update( {"model":2} )
        params.update( {"barrier":50.0} )

        if   case == 0:
            params.update( {"q0":1.0, "p0":0.0, "sq0":0.25} )
            params.update( {"k":0.032, "D":0.1, "x0":0.0, "mass":2000.0} )
        if   case == 1:
            params.update( {"q0":1.5, "p0":0.0, "sq0":0.25} )
            params.update( {"k":0.032, "D":0.1, "x0":0.0, "mass":2000.0} )
        if   case == 2:
            params.update( {"q0":2.0, "p0":0.0, "sq0":0.25} )
            params.update( {"k":0.032, "D":0.1, "x0":0.0, "mass":2000.0} )


    elif model == 3:
        # 1D sin potential
        params.update( {"model":3,"mass":2000.0} )

        if   case == 0:
            params.update( {"A":0.0030, "B":0.5*math.pi, "C":0.0030, "k":1.0, "q0":0.0, "p0":0.0, "sq0":0.25, "sp0":0.0} )
        elif case == 1:
            params.update( {"A":0.0025, "B":0.5*math.pi, "C":0.0025, "k":1.0, "q0":0.0, "p0":0.0, "sq0":0.25, "sp0":0.0} )
        elif case == 2:
            params.update( {"A":0.0020, "B":0.5*math.pi, "C":0.0020, "k":1.0, "q0":0.0, "p0":0.0, "sq0":0.25, "sp0":0.0} )
        elif case == 3:
            params.update( {"A":0.0015, "B":0.5*math.pi, "C":0.0015, "k":1.0, "q0":0.0, "p0":0.0, "sq0":0.25, "sp0":0.0} )
        elif case == 4:
            params.update( {"A":0.0010, "B":0.5*math.pi, "C":0.0010, "k":1.0, "q0":0.0, "p0":0.0, "sq0":0.25, "sp0":0.0} )

        params.update( {"barrier":2.0*math.pi/params["k"]} )


    elif model == 4:
        # 1D eckart barrier
        params.update( {"model":4,"mass":2000.0} )
        params.update( {"barrier":0.0} )

        if   case == 0:
            params.update( {"Va":0.00625, "q0":-3.0, "p0":0.0, "sq0":0.25} )
        elif case == 1:
            params.update( {"Va":0.01250, "q0":-3.0, "p0":3.0, "sq0":0.25} )
        elif case == 2:
            params.update( {"Va":0.01875, "q0":-3.0, "p0":3.0, "sq0":0.25} )
        elif case == 3:
            params.update( {"Va":0.02500, "q0":-3.0, "p0":3.0, "sq0":0.25} )

    elif model == 5:
        # 1D eckart barrier
        params.update( {"model":5,"mass":2000.0} )
        params.update( {"barrier":0.67} )

        if   case == 0:
            params.update( {"q0":-0.2,  "p0":0.0, "sq0":0.1} )
            params.update( {"b":0.2981, "w": 0.01} ) 


    params.update( {"ETHD3_alpha":1.0} )
    return params


def run1D(nsnaps, nsteps):

    models   = [1]
    cases    = [1,2,3]
    ent_opts = [0,1]

    for model in models:
        for case in cases:

            # Initialize params dictionary
            params = init_params(model, case)
            params.update({"nsteps":nsteps, "nsnaps":nsnaps})

            print params
            #sys.exit(0)

            # run quanutm  
            if model == 1:
                # run analytical
                harmonic.run_analytical(params, case)

            else:
                # run numerically exact
                run_exact(params, case)

            # run classical or quantum-classical
            for ent_opt in ent_opts:
                test_ethd_1D.run_QC_1D(params, case, ent_opt)                 

            plot_pes_1D(params,case)


    # Post-Processing
    #test_ethd_1D.make_fig(models,ent_opts,cases)

run1D(100, 10)



