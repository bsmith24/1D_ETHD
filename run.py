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

import test_ho

class tmp:
    pass

def model1(q, params, full_id):
    """
    Harmonic potential
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
    #x0,k,D,V = params["q0"], params["k"], params["D"], params["V"]
    k = 1.0

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
    Symmetric Double Well Potential
    Hdia = 0.25*q^4 - 0.5*q^2   
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
    x2 = x*x

    Hdia.set(0,0, (0.25*x2*x2 - 0.5*x2)*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, x*(x2 - 1.0)*(1.0+0.0j) )

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



def potential(q, params):
    """
    Thin wrapper of the model Hamiltonians that can be used in 
    the fully-quantum calculations
    """

    return compute_model(q, params, Py2Cpp_int([0,0]))
    



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
        obj = compute_model(q, params, full_id)

        f = open("_pes"+str(case)+".txt", "a")
        f.write( "%8.5f   %8.5f   %8.5f \n" % (x, obj.ham_dia.get(0,0).real, obj.d1ham_dia[0].get(0,0).real))
        f.close()


def run_exact(nsnaps, nsteps, params, case):
    """
    The main routine to run fully quantum calculations
    """

    dt = params["dt"]
    # Here we initialize the grid and wavefunction

    # 1-D Initilization
    wfc = Wfcgrid(-50.0, 50.0, 0.01, 1)
    wfc.init_wfc_1D(params["q0"], params["p0"], params["sq0"], 0)
    wfc.update_potential_1D(potential, params)
    wfc.update_propagator_1D(0.5*dt, params["mass"])  # this is important because we are using exp(-0.5*dt*H_loc)...
    wfc.update_propagator_K_1D(dt, params["mass"])    # ... together with exp(-dt*H_non-loc)    
 
    wfc.normalize_1D()

    f = open("_pops"+str(case)+".txt", "w")
    f.close()
  
    flux = doubleList()
    flux.append(0.0)          # len(flux[i]) = nstates
    tot_flux = 0.0
    cum = 0.0
    os.system("mkdir _res"+str(case)+"")
    for i in xrange(nsnaps):  # time steps

        wfc.print_wfc_1D("_res"+str(case)+"/wfc", i, 0)   # Be sure to make a "res" directory
        #wfc.print_reci_wfc_1D("res/reci_wfc", i*dt*nsteps, 0)

        epot = wfc.e_pot_1D()
        ekin = wfc.e_kin_1D(params["mass"])
        etot = epot + ekin
        x, px   = wfc.get_x_1D(),  wfc.get_px_1D()
        x2, px2 = wfc.get_x2_1D(), wfc.get_px2_1D()

        f = open("_pops"+str(case)+".txt", "a")
        f.write("%8.5f   %8.5f   %8.5f   %8.5f  %8.5f  %8.5f    %8.5f    %8.5f  %8.5f  %8.5f\n" % (i*nsteps*dt, ekin, epot, etot, cum, x, px, x2, px2, wfc.norm_1D())) 
        f.close()

        for j in xrange(nsteps):  # time steps 
            wfc.propagate_exact_1D(0)

            res0 = Py2Cpp_double([0.0])
            res1 = Py2Cpp_double([0.0])
            wfc.flux_1D( params["barrier"], res0, params["mass"])
            wfc.flux_1D(-params["barrier"], res1, params["mass"])
            cum = cum + res0[0]*dt - res1[0]*dt


def init_params(model, case):
    """
    This function intiializes the params dictionary for a given case    
    """

    params = {"mass":2000.0, "dt":1.0}

    if   model == 2:
        # Double Well
        if   case == 0:
            params.update( {"model":2, "q0":-1.1, "p0":0.0, "sq0":0.04, "sp0":0.0} )
        elif case == 1:
            params.update( {"model":2, "q0":-1.1, "p0":1.0, "sq0":0.04, "sp0":0.0} )
        if   case == 2:
            params.update( {"model":2, "q0":-1.1, "p0":2.0, "sq0":0.04, "sp0":0.0} )
        elif case == 3:
            params.update( {"model":2, "q0":-1.1, "p0":3.0, "sq0":0.04, "sp0":0.0} )
        if   case == 4:
            params.update( {"model":2, "q0":-1.1, "p0":4.0, "sq0":0.04, "sp0":0.0} )
        elif case == 5:
            params.update( {"model":2, "q0":-1.1, "p0":5.0, "sq0":0.04, "sp0":0.0} )
        if   case == 6:
            params.update( {"model":2, "q0":-1.1, "p0":6.0, "sq0":0.04, "sp0":0.0} )
        elif case == 7:
            params.update( {"model":2, "q0":-1.1, "p0":7.0, "sq0":0.04, "sp0":0.0} )

    elif model == 3:        

        params.update( {"model":3} )

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
        # 1-D Eckart Barrier 
        params.update( {"model":4, "barrier":0.0} )

        if   case == 0:
            params.update( {"Va":0.00625, "q0":-1.0, "p0":1.0, "sq0":0.25, "sp0":0.0} )
        elif case == 1:
            params.update( {"Va":0.00625, "q0":-1.0, "p0":2.0, "sq0":0.25, "sp0":0.0} )
        elif case == 2:
            params.update( {"Va":0.00625, "q0":-1.0, "p0":3.0, "sq0":0.25, "sp0":0.0} )
        elif case == 3:
            params.update( {"Va":0.00625, "q0":-1.0, "p0":4.0, "sq0":0.25, "sp0":0.0} )
        elif case == 4:
            params.update( {"Va":0.00625, "q0":-1.0, "p0":5.0, "sq0":0.25, "sp0":0.0} )
        elif case == 5:
            params.update( {"Va":0.00625, "q0":-1.0, "p0":6.0, "sq0":0.25, "sp0":0.0} )
        elif case == 6:
            params.update( {"Va":0.00625, "q0":-1.0, "p0":7.0, "sq0":0.25, "sp0":0.0} )
        elif case == 7:
            params.update( {"Va":0.00625, "q0":-1.0, "p0":8.0, "sq0":0.25, "sp0":0.0} )


    return params


def run1D(nsnaps, nsteps):
    for model in [3]:
        for case in [0,1,2,3,4]:
            params = init_params(model, case)
            run_exact(nsnaps, nsteps, params, case)

            plot_pes_1D(params,case)

nsnaps = 100
nsteps = 50
run1D(nsnaps, nsteps)


