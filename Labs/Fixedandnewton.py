# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 12:42:12 2026

@author: Benja
"""

"""
 This script uses the bisection method to approximate the root of a 
 scalar function.
"""

############################################# 
"""
Copyright (C) 2025  Adrianna M. Gillman

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
############################################# 

# import libraries
import numpy as np
import matplotlib.pyplot as plt
def driver():

# use routines    
    f = lambda x: np.e**(x**2+7*x-30)
    fp = lambda x: np.e**(x**2+7*x-30)*(2*x+7)
    fp2 = lambda x: 2*np.e**(x**2+7*x-30)*2*x*np.e**(x**2+7*x-30)+7*np.e**(x**2+7*x-30)*(2*x+7)
    a = 2.5
    b = 4

#    f = lambda x: np.sin(x)
#    a = 0.1
#    b = np.pi+0.1

    tol = 1e-6

    [astar,ier] = bisection(f,fp,fp2,a,b,tol)
    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))
    
    p0=astar
    Nmax = 100
    
    (p,pstar,info,it) = newton(f,fp,p0,tol, Nmax)
    print('the approximate root is', '%16.16e' % pstar)
    print('the error message reads:', '%d' % info)
    print('Number of iterations:', '%d' % it)
    return [p,pstar,info,it]

# define routines
def bisection(f,fp,fp2,a,b,tol):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 
    def newt(x):
        return f(x)*fp2(x)/(fp2(x)**2)
    fa = f(a)
    fb = f(b);
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier]

#   verify end points are not a root 



    count = 0
    d = 0.5*(a+b)
    while (abs(newt(d))> tol):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
      
    astar = d
    ier = 0
    print('count = ', count)



def newton(f,fp,p0,tol,Nmax):
      """
      Newton iteration.
      
      Inputs:
        f,fp - function and derivative
        p0   - initial guess for root
        tol  - iteration stops when p_n,p_{n+1} are within tol
        Nmax - max number of iterations
      Returns:
        p     - an array of the iterates
        pstar - the last iterate
        info  - success message
              - 0 if we met tol
              - 1 if we hit Nmax iterations (fail)
         
      """
      p = np.zeros(Nmax+1);
      p[0] = p0
      for it in range(Nmax):
          p1 = p0-f(p0)/fp(p0)
          p[it+1] = p1
          if (abs(p1-p0) < tol):
              pstar = p1
              info = 0
              return [p,pstar,info,it]
          p0 = p1
          pstar = p1
          info = 1
driver()               


x=np.linspace(-5,10)
zeros = [0 for x in x]
def y(x):
    return x-4*np.sin(2*x)-3
plt.plot(x,y(x))
plt.plot(x,zeros)
plt.show()