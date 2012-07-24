import scipy
import scipy.integrate as integrate

dy = .1 * .5 # /int_{\Delta y} y dy = \Delta y / 2
dz = .1 * .5 # /int_{\Delta z} z dz = \Delta z / 2
mua = 4.1
k =  0.57
deltat = 1.0
theta = 0.5

xlo= 0.0 
xhi= 1./20
print xlo,xhi
deltax = xhi - xlo

def uexact(x,t):
  return scipy.exp( -mua *x )/(mua *k) *( scipy.exp( k*mua*mua * t ) - 1. )

def phi1(x):
  return (xhi - x )/deltax 
def phi2(x):
  return (x - xlo )/deltax  

def dphi1dx(x):
  return -1./deltax 
def dphi2dx(x):
  return  1./deltax 

def uh(x,t):
  return phi1(x) * uexact(xlo,t) + phi2(x) * uexact(xhi,t)  

def duhdx(x,t):
  return ( uexact(xhi,t) - uexact(xlo,t) )/deltax 

def q(x):
  return mua * scipy.exp( -mua *x )

def pennesphi(x,phi,dphidx):
  deltau = uh(x,1.0) - uh(x,0.0)
  print x,"deltau =", deltau , "deltax =", theta*(duhdx(x,0.0) + duhdx(x,1.0)) ,"q(x) =", q(x)
  print "phi =", phi(x) , "dphidx =", dphidx(x)
  print  ( deltau / deltat - q(x) ) * phi(x) + k*theta*(duhdx(x,0.0) + duhdx(x,1.0)) *dphidx(x)
  return  ( deltau / deltat - q(x) ) * phi(x) + k*theta*(duhdx(x,0.0) + duhdx(x,1.0)) *dphidx(x)


#print integrate.quadrature(pennesphi1,xlo,xhi)
print  dy * dz * integrate.fixed_quad(pennesphi,xlo,xhi,(phi1,dphi1dx),2)[0] 
print  dy * dz * integrate.fixed_quad(pennesphi,xlo,xhi,(phi2,dphi2dx),2)[0] 
#print integrate.quadrature(pennesphi2,xlo,xhi)
#print integrate.fixed_quad(pennesphi,xlo,xhi,(2),2)

