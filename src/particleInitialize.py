from pylab import *


def particleInitialize(p, n, L):
  """
  addParticle(x, y, z, vx, vy, vz, r,
              thetax, thetay, thetaz, 
              omegax, omegay, omegaz): 
  """
  dx = 2.0*L / n
  d  = linspace(dx/2.0 - L, L - dx/2.0, n)

  for i in d:
    for j in d:
      for k in d:
        p.addParticle(i,j,k,0,0,0,1.0/4,0,0,0,0,0,0)
