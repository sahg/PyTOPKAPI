""" ODE_module.py

Routines used to solve the Ordinary Differential Equations in TOPKAPI

"""

__author__ = "Theo Vischel"
__version__ = "$Revision: 1.1 $"
__date__ = "$Date: 01/15/2007 $"

from copy import *
from math import * 
from scipy import *
from numpy import *


###=======================================###
###          ANALYTICAL SOLUTIONS         ###
###=======================================###
#### Solution if a=0
def input_zero_solution(b, alpha, V0, delta_t):
    V1=(V0**(1-alpha)+b*(alpha-1)*delta_t)**(1/(1-alpha))
    return V1
#### Solution if b=0
def coefb_zero_solution(a, V0, Dt):
    V1=a*Dt+V0
    return V1

###=======================================###
###          RUNGE KUTTA FEHLBERG         ###
###=======================================###
#```````````````````````````````````````````
def fonction(a,b,alpha):
    return lambda x:a-b*x**alpha

###```````````````````````````````````````````    
class RKF:
    """
    Adaptive solver
    """
    def __init__(self, min_step=10e-10, max_step=3600, min_tol=1e-7, max_tol=1e-3, init_time_step=3600):
        self.min_step = min_step
        self.max_step = max_step
        self.min_tol = min_tol    
        self.max_tol = max_tol
        self.resultbuffer = 0
        self.timebuffer = 0
        if init_time_step==0:
            self.delta_t = (max_step + min_step)/2.
        else:
            self.delta_t = init_time_step

    def compute_error(self, f, x, t, delta_t):
        # define all the constants for Runge-Kutta of order 4 and 5
        #From numerical recipes        
        a2 = 1/5.; a3 = 3/10.; a4 = 3/5.; a5 = 1.; a6 = 7/8.
        b21 = 1/5.
        b31 = 3/40.; b32 = 9/40.
        b41 = 3/10.; b42 = -9/10.; b43 = 6/5.
        b51 = -11/54.; b52 = 5/2.; b53 = -70/27.; b54 = 35/27.
        b61 = 1631/55296.; b62 = 175/512.; b63 = 575/13824.; b64 = 44275/110592.; b65 = 253/4096.
        c1 = 37/378.; c2 = 0.; c3 = 250/621.; c4 = 125/594.; c5 = 0.; c6 = 512/1771.
        d1 = 2825/27648.; d2 = 0.; d3 = 18575/48384.; d4 = 13525/55296.; d5 = 277/14336.; d6 = 1/4.

        k1 = delta_t*f(x)
        # compute the k2
        tmpx = x + b21*k1
        if tmpx<0:
            Err1=self.max_tol*1.1
            return Err1
        else:
            k2 = delta_t*f(tmpx)
        # compute the k3
        tmpx = x + b31*k1 + b32*k2
        if tmpx<0:
            Err1=self.max_tol*1.1
            return Err1
        else:
            k3 = delta_t*f(tmpx)
        # compute the k4
        tmpx = x + b41*k1 + b42*k2 + b43*k3
        if tmpx<0:
            Err1=self.max_tol*1.1
            return Err1
        else:
            k4 = delta_t*f(tmpx)
        # compute the k5
        tmpx = x + b51*k1 + b52*k2 + b53*k3 + b54*k4
        if tmpx<0:
            Err1=self.max_tol*1.1
            return Err1
        else:
            k5 = delta_t*f(tmpx)
        # compute the k6
        tmpx = x + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5
        if tmpx<0:
            Err1=self.max_tol*1.1
            return Err1
        else:
            k6 = delta_t*f(tmpx)
            
        # compute the new variables
        newx_rk5 = x + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6
        newx_rk4 = x + d1*k1 + d2*k2 + d3*k3 + d4*k4 + d5*k5 + d6*k6

        # Compute the error
        Err1=abs(newx_rk4-newx_rk5)

        self.resultbuffer = deepcopy(newx_rk5)
        self.timebuffer = t+delta_t

        return Err1

    def getnewdelta_t(self, f, x, t, max_delta_t=0):
##        print 'getnewdelta_t'
        if max_delta_t==0:
            max_delta_t = self.max_step
        if self.delta_t > max_delta_t:
            self.delta_t = max_delta_t

        done = 0
##        ntry=-1
        factor='initial'
        while done==0:
            error = self.compute_error(f, x, t, self.delta_t)
            if error > self.max_tol:
                if self.delta_t/2. < self.min_step:
                    self.delta_t = self.min_step
                    self.compute_error(f, x, t, self.delta_t)                    
                    done = 1
                else:
                    self.delta_t = self.delta_t / 2.
                    factor='already divided'
            elif error < self.min_tol:
                if factor=='already divided':
                    #Avoid to infintively loop by halving then multiplying by two then halving...
                    print 'The last time step was kept to avoid overlooping'
                    done = 1
                elif self.delta_t*2 > max_delta_t:
                    # do not go over the specified time step
                    if self.delta_t-max_delta_t < self.min_step:
                        done = 1
                    else:
                        self.delta_t = max_delta_t
                        self.compute_error(f, x, t, self.delta_t)
                        done = 1
                elif self.delta_t*2 > self.max_step:
                    self.delta_t = self.max_step
                    self.compute_error(f, x, t, self.delta_t)
                    done = 1
                else:
                    self.delta_t = self.delta_t*2
            else:
                done = 1
        return self.delta_t

    def step(self, f, x, t, delta_t=0):
        if delta_t==0:
	    # we should have called getnewdelta_t before, to obtain the time
	    # step used and hence the solution should already be computed and
	    # put in resultbuffer            
            if self.timebuffer!=t+self.delta_t:
#                print 'WARNING: This was not expected to ever occur!!!!'
                self.getnewdelta_t(f, x, t, delta_t)
            return self.resultbuffer
        else:
            # run the adaptive integration until we reach t+delta_t
            curtime = t
            curx = x
            while curtime < t+delta_t:
                self.getnewdelta_t(f, curx, curtime, (t+delta_t)-curtime)
                curx = self.step(f, curx, curtime)
                curtime = curtime + self.delta_t
            return curx

###=======================================###
###       QUASI ANALYTICAL SOLUTION       ###
###=======================================###
        
def qas(a, b, alpha, V0, delta_t,derivative=0):
    if alpha > 2:
        exposant=alpha/(alpha-1)
        y0=V0**(1-alpha)
        a_eq=b*(alpha-1)
        b_eq=a*(alpha-1)
    else:
        exposant=alpha
        y0=V0
        a_eq=a
        b_eq=b
        
    #   Definition of the variables alpha0 and beta0
    #   that approximate y**(exposant-1)=alpha0+beta0*y
    if derivative==1:
    ##    First solution with derivative function
        alpha0,beta0=adjust_derivative_line(exposant,y0)
    else:
    ##    Second solution by roughly estimating y1 and adjust the line
        alpha0,beta0=adjust_2points_line(a_eq,b_eq,exposant,y0,delta_t)

    if alpha > 2:
        A=-a*(alpha-1)*beta0
        B=alpha0/beta0
        C=-b/(a*beta0)
    else:
        A=-b*beta0
        B=alpha0/beta0
        C=-a/(b*beta0)
        
    
    y1=solution(A,B,C,y0,delta_t)
    if alpha > 2:
        V1=y1**(1/(1-alpha))
    else:
        V1=y1
        
    return V1

def adjust_2points_line(a_eq,b_eq,exposant,y0,delta_t):
    #Definition of the variables alpha0 and beta0
    #that approximate y**(exposant-1)=alpha0+beta0*y
    f=fonction(a_eq,b_eq,exposant)
    y1_estimat=y0+delta_t*a_eq #-->OK 1.
    beta0=(y1_estimat**(exposant-1)-y0**(exposant-1))/(y1_estimat-y0)
    alpha0=y0**(exposant-1)-beta0*y0
    return alpha0,beta0

def solution(A,B,C,y0,delta_t):
    p1=1/2.*(-B+sqrt(B**2-4*C))
    p2=1/2.*(-B-sqrt(B**2-4*C))
    y1=(p1-p2*((y0-p1)/(y0-p2)*exp(A*delta_t*(p1-p2))))\
        /(1-((y0-p1)/(y0-p2)*exp(A*delta_t*(p1-p2))))
    return y1                                       
