
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import matplotlib.pyplot as plt

#Defining the function to integrate
def func(x):
    return (np.exp(-2*x))*(np.cos(10*x))

#Defining the integral
def func_integral(x):
    return (np.exp(-2*x))*(5*np.sin(10*x) - np.cos(10*x))/(52)


# TRAPEZOID METHOD

#This is the core program of the Trapezoid Method
def trapezoid_core(f,x,h):
    return 0.5*h*(f(x+h) + f(x))

#Wrapper function for the Trapezoid Method
def trapezoid_method(f,a,b,N):
    #f == function to integrate
    #a == lower limit of integration
    #b == upper limit of integration
    #N == number of function evaluations to use
    
    #define x values to perform trapezoid rule
    x = np.linspace(a,b,N)
    h = x[1]-x[0]
    
    #define the value of the integral
    Fint = 0.0
    
    #perform the integral using the trapezoid method
    for i in range(0, len(x)-1,1):
        Fint += trapezoid_core(f,x[i],h)
        
    #return the answer
    return Fint

#SIMPSON'S METHOD

#This is the core program of Simpson's Method
def simpson_core(f,x,h):
    return h*( f(x) + 4*f(x+h) + f(x+2*h))/3.

#Wrapper function for Simpson's Method
def simpsons_method(f,a,b,N):
    #f == function to integrate
    #a == lower limit of integration
    #b == upper limit of integration
    #N == number of function evaluations to use
    
    '''Note the number of chunks will be N-1
        so if N is odd, then we don't need to
        adjust the last segment'''
    
    #Define x values to perform simpson's rule
    x = np.linspace(a,b,N)
    h = x[1]-x[0]

    #define the value of the integral
    Fint = 0.0
    
    #perform the integral using simpson's method
    for i in range(0,len(x)-2,2):
        Fint += simpson_core(f,x[i],h)
        
    #apply simpson's rule over the last interval
    #if N is even
    if((N%2)==0):
        Fint += simpson_core(f,x[-2],0.5*h)
        
    return Fint

# ROMBERG INTEGRATION

#This is the Romberg Core
def romberg_core(f,a,b,i):
    
    #we need the difference b-a
    h = b-a
    
    #and the increment between new func evals
    dh = h/2.**(i)
    
    #we need the cofactor
    K = h/2.**(i+1)
    
    #and the function evaluations
    M = 0.0
    for j in range(2**i):
        M += f(a + 0.5*dh + j*dh)
        
    #return the answer
    return K*M

#Wrapper function for Romberg
def romberg_integration(f,a,b,tol):
    
    #define an iteration variable
    i = 0
    
    #define a maximum number of iterations
    imax = 1000
    
    #define an error estimate, set to a large value
    delta = 100.0*np.fabs(tol)
    
    #set an array of integral answers
    I = np.zeros(imax,dtype=float)
    
    #get the zeroth romberg iteration
    I[0] = 0.5*(b-a)*(f(a) + f(b))
    
    #iterate by 1
    i += 1
    
    while(delta>tol):
        
        #find this romberg iteration
        I[i] = 0.5*I[i-1] + romberg_core(f,a,b,i)
        
        #compute the new fractional error estimate
        delta = np.fabs( (I[i]-I[i-1])/I[i])
        
        print(i,I[i],I[i-1],delta)
        
        if(delta>tol):
            
            #iterate
            i+=1
            
            #if we've reached the maximum iterations
            if(i>imax):
                print("Max iterations reached.")
                raise StopIteration('Stopping iterations after ',i)
                
    #return the answer
    return I[i]

#Print Answers:

print("Answer:")
Answer = func_integral(np.pi)-func_integral(0)
print(Answer)
print("Trapezoid:")
print(trapezoid_method(func,0,1,10))
print("Simpson's Method:")
print(simpsons_method(func,0,1,10))
print("Romberg:")
tolerance = 1.0e-6
RI = romberg_integration(func,0,1,tolerance)      
print(RI, (RI-Answer)/Answer, tolerance)

