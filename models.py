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

class tmp:
    pass


def model0(q, params, full_id):
    """
    Constant potential
    Hdia = 0.0   
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

    Hdia.set(0,0, 0.0)
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0,0.0)

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j)

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


def model1(q, params, full_id):
    """
    Harmonic potential
    Hdia = 0.5*k*(x-x0)^2   
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
    x0, k = params["wfc"]["x0"][0], params["k"]

    Hdia.set(0,0, 0.5*k*(x-x0)*(x-x0)*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, k*(x-x0)*(1.0+0.0j) )

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


def model5(q, params, full_id):
    """
    #Params in:  q = list of trajectory positions
    Cubic potential with switching 
    Hdia = A*q^2 - B*q^3   
    Sdia = 1.0
    Ddia = 0.0
    """

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(1,1) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(1,1) )

    b, ww = 0.2981, 0.01
    m = params["mass"]
    mw2 = m * ww * ww
         
    x = q.col(indx).get(0)
    x2 = x*x

    fun = 0.5*mw2*x2 - (b/3.0)*x2*x
    dfun = mw2*x - b*x2

    sw = SWITCH(VECTOR(x,0.0,0.0),VECTOR(0.0, 0.0, 0.0), 0.9, 1.5);

    x = 1.5
    x2 = x * x
    fun_const = 0.5*mw2*x2 - (b/3.0)*x2*x
    dfun_const = mw2*x - b*x2

    en = fun*sw[0] + fun_const*(1.0 - sw[0])
    den = dfun *sw[0] + fun*sw[1].x - fun_const*sw[1].x
   

    Hdia.set(0,0, en*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, den*(1.0+0.0j) )

        #  <dia| d/dR_0| dia >
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

    if model==0:
        # Constant Potential
        res = model0(q, params, full_id)
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
    if model==5:
        # 1D Cubic Switching
        res = model5(q, params, full_id)

    return res

