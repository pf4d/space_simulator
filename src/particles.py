from pylab import *

# 9/23/08 JVJ and Tim Bocek

class GranularMaterialForce(object):

  def __init__(self, k=1.5, gamma=0.3):
    # parameters in force model
    self.k     = k           # elastic 'bounce'
    self.gamma = gamma       # energy dissipation/loss
    self.f     = 0.5         # angular exchange damping coefficient
    self.rho   = 1.0         # density
    self.G     = 6.67384E-11 # gravitational constant

  def __call__(self, p):
    # find position differences :
    d, dx, dy, dz = p.distanceMatrix(p.x, p.y, p.z)

    # compute overlap :
    dr = d - p.sumOfRadii
    fill_diagonal(dr, 0)

    # no forces arising in no overlap cases :
    dr[dr > 0] = 0

    # compute elastic particle/particle forces :
    mag_r = self.k * dr

    # velocity differences :
    dv, dvx, dvy, dvz = p.distanceMatrix(p.vx, p.vy, p.vz)
    da, dax, day, daz = p.distanceMatrix(p.ax, p.ay, p.az)
    
    # damping terms :
    vijDotrij        = dvx*dx + dvy*dy + dvz*dz
    vijDotrij[dr==0] = 0 

    # damping is subtracted from force :
    mag_r += self.gamma * vijDotrij / d

    # gravitational pull between particles :
    a_g = self.G * p.mi / d**2
    a_g[d < 0.0] = 0

    # Project onto components, sum all forces on each particle
    p.ax = sum(mag_r * dx/d * p.ratioOfRadii + a_g*dx/d, axis=1)
    p.ay = sum(mag_r * dy/d * p.ratioOfRadii + a_g*dy/d, axis=1)
    p.az = sum(mag_r * dz/d * p.ratioOfRadii + a_g*dz/d, axis=1)
    
    # differences in angular velocities :
    do, dox, doy, doz = p.distanceMatrix(p.omegax, p.omegay, p.omegaz)
    dox[dr==0] = 0
    doy[dr==0] = 0
    doz[dr==0] = 0
    do[dr==0]  = 0
    
    # projection of a onto the tangent plane to r (tangential acceleration) :
    atx = dax - (dax * dx) / d**2 * dx
    aty = day - (day * dy) / d**2 * dy
    atz = daz - (daz * dz) / d**2 * dz
    
    # projection of v onto the tangent plane to r (tangential velocity) : 
    vtx = dvx - (dvx * dx) / d**2 * dx
    vty = dvy - (dvy * dy) / d**2 * dy
    vtz = dvz - (dvz * dz) / d**2 * dz
    
    # calculate the radius vector to the point of torque :
    radx = dx/d * p.r
    rady = dy/d * p.r
    radz = dz/d * p.r
   
    # no angular forces where particles do not touch :
    atx[dr==0] = 0
    aty[dr==0] = 0
    atz[dr==0] = 0
    vtx[dr==0] = 0
    vty[dr==0] = 0
    vtz[dr==0] = 0

    # calculate torque (r x F) :
    taux = rady*atz - radz*aty
    tauy = radz*atx - radx*atz
    tauz = radx*aty - rady*atx

    # calculate tangential velocity parallel to torque (r x vt) :
    epix = rady*vtz - radz*vty
    epiy = radz*vtx - radx*vtz
    epiz = radx*vty - rady*vtx

    # angular acceleration coefficient:
    k = 0.5

    # linear velocity frictional contribution to torque coefficient :
    f = 0.0

    # angular velocity damping coefficient from medium (air) :
    g = 0.0

    # moment of inertia for a sphere :
    I = 0.4*p.r**2
    
    # angular momentum exchange coefficient (damping term) :
    #oijDotdij        = dox*dx + doy*dy + doz*dz
    #oijDotdij[dr==0] = 0
    kappa            = self.f * p.r
   
    # project onto components, sum all angular forces on each particle
    p.alphax = sum((k*taux + f*epix) / I + kappa*dox, axis=1) - g*p.omegax
    p.alphay = sum((k*tauy + f*epiy) / I + kappa*doy, axis=1) - g*p.omegay
    p.alphaz = sum((k*tauz + f*epiz) / I + kappa*doz, axis=1) - g*p.omegaz
   

class NebulaGranularMaterialForce(object):

  def __init__(self, k=1.5, gamma=0.3):
    # parameters in force model
    self.k     = k           # elastic 'bounce'
    self.gamma = gamma       # energy dissipation/loss

  def __call__(self, p):
    # find position differences :
    d, dx, dy, dz = p.distanceMatrix(p.x, p.y, p.z)

    # compute overlap :
    dr = d - p.sumOfRadii
    dr[0] = 0

    # no forces arising in no overlap cases :
    dr[dr > 0] = 0

    # compute elastic particle/particle forces :
    mag_r = self.k * dr

    # velocity differences :
    dv, dvx, dvy, dvz = p.distanceMatrix(p.vx, p.vy, p.vz)
    da, dax, day, daz = p.distanceMatrix(p.ax, p.ay, p.az)
    
    # damping terms :
    vijDotrij         = dvx*dx + dvy*dy + dvz*dz
    vijDotrij[dr==0]  = 0 

    # damping is subtracted from force :
    mag_r += self.gamma * vijDotrij / d

    # Project onto components, sum all forces on each particle
    p.ax = mag_r * dx/d * p.ratioOfRadii
    p.ay = mag_r * dy/d * p.ratioOfRadii
    p.az = mag_r * dz/d * p.ratioOfRadii
    

class VerletIntegrator(object):

  def __init__(self, dt=0.01):
    # time step
    self.dt = dt

  def __call__(self, force, p):
    dt = self.dt

    # position update
    p.x = p.x + p.vx*dt + 0.5*p.ax*dt**2
    p.y = p.y + p.vy*dt + 0.5*p.ay*dt**2
    p.z = p.z + p.vz*dt + 0.5*p.az*dt**2

    # angular roation update :
    p.thetax = p.omegax*dt + 0.5*p.alphax*dt**2
    p.thetay = p.omegay*dt + 0.5*p.alphay*dt**2
    p.thetaz = p.omegaz*dt + 0.5*p.alphaz*dt**2
    p.update_theta()

    # update periodic BC
    p.pbcUpdate()

    # store accelerations for averaging that is done 
    ax = p.ax
    ay = p.ay
    az = p.az
    
    # store angular accelerations for crank-nicolson scheme :
    alphax = p.alphax
    alphay = p.alphay
    alphaz = p.alphaz

    force(p) # force update with new positions

    # velocity updates
    p.vx = p.vx + 0.5*(ax + p.ax)*dt
    p.vy = p.vy + 0.5*(ay + p.ay)*dt
    p.vz = p.vz + 0.5*(az + p.az)*dt

    # angular velocity updates :
    p.omegax = p.omegax + 0.5*(alphax + p.alphax)*dt
    p.omegay = p.omegay + 0.5*(alphay + p.alphay)*dt
    p.omegaz = p.omegaz + 0.5*(alphaz + p.alphaz)*dt


class NebulaVerletIntegrator(object):

  def __init__(self, dt=0.01):
    # time step
    self.dt = dt

  def __call__(self, force, p):
    dt = self.dt

    # position update
    p.x = p.x + p.vx*dt + 0.5*p.ax*dt**2
    p.y = p.y + p.vy*dt + 0.5*p.ay*dt**2
    p.z = p.z + p.vz*dt + 0.5*p.az*dt**2

    # update periodic BC
    p.pbcUpdate()

    # store accelerations for averaging that is done 
    ax = p.ax
    ay = p.ay
    az = p.az
    
    force(p) # force update with new positions

    # velocity updates
    p.vx = p.vx + 0.5*(ax + p.ax)*dt
    p.vy = p.vy + 0.5*(ay + p.ay)*dt
    p.vz = p.vz + 0.5*(az + p.az)*dt


class Particles(object):

  def __init__(self, L, force, periodicX=1, periodicY=1, periodicZ=1):
    """
    """
    # container size
    self.L = L
    # total Number of particles
    self.N = 0
    # type
    self.type = 'float32'
    # positions
    self.x  = array([],dtype=self.type)
    self.y  = array([],dtype=self.type)
    self.z  = array([],dtype=self.type)
    # velocities
    self.vx = array([],dtype=self.type)
    self.vy = array([],dtype=self.type)
    self.vz = array([],dtype=self.type)
    # forces
    self.ax = array([],dtype=self.type)
    self.ay = array([],dtype=self.type)
    self.az = array([],dtype=self.type)
    # radii
    self.r  = array([],dtype=self.type)
    # angular rotation :
    self.theta  = []
    self.thetax = array([],dtype=self.type)
    self.thetay = array([],dtype=self.type)
    self.thetaz = array([],dtype=self.type)
    # angular velocity :
    self.omegax = array([],dtype=self.type)
    self.omegay = array([],dtype=self.type)
    self.omegaz = array([],dtype=self.type)
    # angular acceleration :
    self.alphax = array([],dtype=self.type)
    self.alphay = array([],dtype=self.type)
    self.alphaz = array([],dtype=self.type)
    # periodic on?
    self.periodicX = periodicX 
    self.periodicY = periodicY 
    self.periodicZ = periodicZ 
    # force function :
    self.f = force
    # mass / gravity :
    self.m     = array([],dtype=self.type)
    self.V     = array([],dtype=self.type)
    self.mi    = array([],dtype=self.type)
    self.rho   = array([],dtype=self.type)
     
  def addParticle(self, x, y, z, vx, vy, vz, r, rho,
                  thetax, thetay, thetaz, 
                  omegax, omegay, omegaz): 
    """
    """
    self.x   = hstack((self.x,x))
    self.y   = hstack((self.y,y))
    self.z   = hstack((self.z,z))
    self.vx  = hstack((self.vx,vx))
    self.vy  = hstack((self.vy,vy))
    self.vz  = hstack((self.vz,vz))
    self.ax  = hstack((self.ax,0))
    self.ay  = hstack((self.ay,0))
    self.az  = hstack((self.az,0))
    self.r   = hstack((self.r,r))
    self.rho = hstack((self.rho,rho))
    I        = identity(3)
    R        = self.rotate(I, array([thetax, thetay, thetaz]))
    self.theta.append(R)
    self.thetax = hstack((self.thetax,thetax))
    self.thetay = hstack((self.thetay,thetay))
    self.thetaz = hstack((self.thetaz,thetaz))
    self.omegax = hstack((self.omegax,omegax))
    self.omegay = hstack((self.omegay,omegay))
    self.omegaz = hstack((self.omegaz,omegaz))
    self.alphax = hstack((self.alphax,0))
    self.alphay = hstack((self.alphay,0))
    self.alphaz = hstack((self.alphaz,0))
    self.N  = self.N+1
    temp    = tile(self.r,(self.N,1))
    self.sumOfRadii   = temp + temp.T
    self.ratioOfRadii = temp / temp.T
    # gravitational pull :
    V          = 4.0/3.0 * pi * r**3
    self.V     = hstack((self.V, V))
    self.m     = self.rho * self.V
    self.mi    = tile(self.m, (self.N,1))
    fill_diagonal(self.mi, zeros(self.N))
    self.f(self)

  def update_theta(self):
    """
    """
    theta_n = []
    for M,x,y,z in zip(self.theta,self.thetax,self.thetay,self.thetaz):
      v = array([x,y,z])
      M = self.rotate(M, v)
      theta_n.append(M)
    self.theta = theta_n

  def rotate(self, M, v):
    """
    rotate the particle's orientation matrix <M> about the x, y, and z axes by 
    angles provided in <v> array.
    """
    rx = v[0]
    ry = v[1]
    rz = v[2]
    c  = cos(rx)
    s  = sin(rx)
    Rx = array([[1, 0,  0],
                [0, c, -s],
                [0, s,  c]])
    c  = cos(ry)
    s  = sin(ry)
    Ry = array([[ c, 0, s],
                [ 0, 1, 0],
                [-s, 0, c]])
    c  = cos(rz)
    s  = sin(rz)
    Rz = array([[c, -s, 0],
                [s,  c, 0],
                [0,  0, 1]])
    R  = dot(Rx, dot(Ry, Rz))
    return dot(M, R)

  def pbcUpdate(self):
    """
    Moves paricles across periodic boundary
    """
    if self.periodicX:
      self.x[self.x >  self.L/2] = self.x[self.x >  self.L/2] - self.L
      self.x[self.x < -self.L/2] = self.x[self.x < -self.L/2] + self.L
    if self.periodicY:
      self.y[self.y >  self.L/2] = self.y[self.y >  self.L/2] - self.L
      self.y[self.y < -self.L/2] = self.y[self.y < -self.L/2] + self.L
    if self.periodicZ:
      self.z[self.z >  self.L/2] = self.z[self.z >  self.L/2] - self.L
      self.z[self.z < -self.L/2] = self.z[self.z < -self.L/2] + self.L

  def distanceMatrix(self, x, y, z):
    """
    Computes distances between all particles and places the result in a 
    matrix such that the ij th matrix entry corresponds to the distance 
    between particle i and j
    """ 
    xtemp = tile(x, (self.N,1))
    ytemp = tile(y, (self.N,1))
    ztemp = tile(z, (self.N,1))
    dx    = xtemp - xtemp.T
    dy    = ytemp - ytemp.T
    dz    = ztemp - ztemp.T
  
    # Particles 'feel' each other across the periodic boundaries
    if self.periodicX:
      dx[dx >  self.L/2] = dx[dx >  self.L/2] - self.L
      dx[dx < -self.L/2] = dx[dx < -self.L/2] + self.L
    if self.periodicY:
      dy[dy >  self.L/2] = dy[dy >  self.L/2] - self.L
      dy[dy < -self.L/2] = dy[dy < -self.L/2] + self.L
    if self.periodicZ:
      dz[dz >  self.L/2] = dz[dz >  self.L/2] - self.L
      dz[dz < -self.L/2] = dz[dz < -self.L/2] + self.L

    # Total Distances
    d = sqrt(dx**2 + dy**2 + dz**2)

    # Mark zero entries with negative 1 to avoid divergences
    d[d==0] = -1

    return d, dx, dy, dz


class Nebula(Particles):

  def __init__(self, L, force, periodicX=1, periodicY=1, periodicZ=1):
    super(Nebula, self).__init__(L, force, periodicX, periodicY, periodicZ)
     
  def addParticle(self, x, y, z, vx, vy, vz, r): 
    self.x   = hstack((self.x,x))
    self.y   = hstack((self.y,y))
    self.z   = hstack((self.z,z))
    self.vx  = hstack((self.vx,vx))
    self.vy  = hstack((self.vy,vy))
    self.vz  = hstack((self.vz,vz))
    self.ax  = hstack((self.ax,0))
    self.ay  = hstack((self.ay,0))
    self.az  = hstack((self.az,0))
    self.thetax = hstack((self.thetax,0))
    self.thetay = hstack((self.thetay,0))
    self.thetaz = hstack((self.thetaz,0))
    self.omegax = hstack((self.omegax,0))
    self.omegay = hstack((self.omegay,0))
    self.omegaz = hstack((self.omegaz,0))
    self.alphax = hstack((self.alphax,0))
    self.alphay = hstack((self.alphay,0))
    self.alphaz = hstack((self.alphaz,0))
    self.theta.append(identity(3))
    self.r   = hstack((self.r,r))
    self.N   = self.N+1
    self.sumOfRadii   = self.r[0] + self.r
    self.ratioOfRadii = self.r[0] / self.r
    self.f(self)

  def update_theta(self):
    """
    """
    super(Nebula, self).update_theta()

  def rotate(self, M, v):
    """
    rotate the particle's orientation matrix <M> about the x, y, and z axes by 
    angles provided in <v> array.
    """
    return super(Nebula, self).rotate(M, v)

  def pbcUpdate(self):
    """
    Moves paricles across periodic boundary
    """
    super(Nebula, self).pbcUpdate()

  def distanceMatrix(self, x, y, z):
    """
    Computes distances between all particles and places the result in a 
    matrix such that the ij th matrix entry corresponds to the distance 
    between particle i and j
    """
    dx = x[0] - x
    dy = y[0] - y
    dz = z[0] - z
  
    # Particles 'feel' each other across the periodic boundaries
    if self.periodicX:
      dx[dx >  self.L/2] = dx[dx >  self.L/2] - self.L
      dx[dx < -self.L/2] = dx[dx < -self.L/2] + self.L
    if self.periodicY:
      dy[dy >  self.L/2] = dy[dy >  self.L/2] - self.L
      dy[dy < -self.L/2] = dy[dy < -self.L/2] + self.L
    if self.periodicZ:
      dz[dz >  self.L/2] = dz[dz >  self.L/2] - self.L
      dz[dz < -self.L/2] = dz[dz < -self.L/2] + self.L

    # Total Distances
    d = sqrt(dx**2 + dy**2 + dz**2)

    # Mark zero entry with negative 1 to avoid divergences
    d[0] = -1

    return d, dx, dy, dz


class Specter(object):

  def __init__(self, n, L):
    """
    """
    self.n = int(n)       # number of points
    self.L = L            # extent of field
    self.x = array([])    # x-coords
    self.y = array([])    # y-coords
    self.z = array([])    # z-coords
    for i in range(self.n):
      #x,y,z = L*rand(3) - L/2.0
      x,y,z = L*randn(3)
      self.addSpecter(x,y,z)
   
  def addSpecter(self, x, y, z): 
    """
    """
    self.x  = hstack((self.x,x))
    self.y  = hstack((self.y,y))
    self.z  = hstack((self.z,z))


def initialize_grid(p, n, r, rho, L):
  """
  addParticle(x, y, z, vx, vy, vz, r, rho,
              thetax, thetay, thetaz, 
              omegax, omegay, omegaz): 
  """
  dx = 2.0*L / n
  d  = linspace(dx/2.0 - L, L - dx/2.0, n)

  p.addParticle(0,0,2*L,0,0,0,3,10,0,0,0,0,0,0)
  for i in d:
    for j in d:
      for k in d:
        p.addParticle(i,j,k,0,0,0,r,rho,0,0,0,0,0,0)


def initialize_random(p, n, r, rho, L):
  """
  addParticle(x, y, z, vx, vy, vz, r, rho,
              thetax, thetay, thetaz, 
              omegax, omegay, omegaz): 
  """
  p.addParticle(0,0,2*L,0,0,0,3,rho,0,0,0,0,0,0)
  for i in range(n):
    r     = 0.1*randn() + r
    x,y,z = L*randn(3)
    p.addParticle(x,y,z,0,0,0,r,rho,0,0,0,0,0,0)

def initialize_system(p):
   """
   """
   G         = 6.67e-11
   au        = 1.49597871e11
   
   r_sun     = 6.955e8
   r_earth   = 6.3781e6
   r_moon    = 1.7374e6
   r_jupiter = 7.1492e7
   r_saturn  = 6.0268e7
   
   d_earth   = au
   d_moon    = 4.0e8
   d_jupiter = 5.2*au
   d_saturn  = 9.5*au

   m_sun     = 1.98844e30
   m_earth   = 5.97219e24
   m_moon    = 7.2377e22
   m_jupiter = 1.8986e27
   m_saturn  = 5.6846e26

   V_sun     = 4/3.0 * pi * r_sun**3
   V_earth   = 4/3.0 * pi * r_earth**3
   V_moon    = 4/3.0 * pi * r_moon**3
   V_jupiter = 4/3.0 * pi * r_jupiter**3
   V_saturn  = 4/3.0 * pi * r_saturn**3

   rho_sun     = m_sun / V_sun
   rho_earth   = m_earth / V_earth
   rho_moon    = m_moon / V_moon
   rho_jupiter = m_jupiter / V_jupiter
   rho_saturn  = m_saturn / V_saturn

   F_earth     = G * m_earth * m_sun / d_earth**2
   F_moon      = G * m_earth * m_moon / d_moon**2
   F_jupiter   = G * m_jupiter * m_sun / d_jupiter**2
   F_saturn    = G * m_saturn * m_sun / d_saturn**2

   vt_earth    = sqrt(G * m_sun / d_earth)
   vt_moon     = sqrt(G * m_earth / d_moon)
   vt_jupiter  = sqrt(G * m_sun / d_jupiter)
   vt_saturn   = sqrt(G * m_sun / d_saturn)

   #addParticle(x, y, z, vx, vy, vz, r, rho,
   #            thetax, thetay, thetaz, 
   #            omegax, omegay, omegaz): 
   p.addParticle(au + 6.3781e6 + 1e6,0,0,0,0,0,3,0,0,0,0,0,0,0)
   p.addParticle(0, 0, 0, 0, 0, 0, r_sun, rho_sun, 0, 0, 0, 0, 0, 0)
   p.addParticle(d_earth, 0, 0, 0, 0, vt_earth, r_earth, rho_earth,
                 0, 0, 0, 0, 0, 0)
   p.addParticle(d_earth+d_moon, 0, 0, 0, 0, vt_moon, r_moon, rho_moon,
                 0, 0, 0, 0, 0, 0)
   p.addParticle(d_jupiter, 0, 0, 0, 0, vt_jupiter, r_jupiter, rho_jupiter,
                 0, 0, 0, 0, 0, 0)
   p.addParticle(d_saturn, 0, 0, 0, 0, vt_saturn, r_saturn, rho_saturn,
                 0, 0, 0, 0, 0, 0)


def initialize_earth(p):
   """
   """
   r_earth   = 6.3781e6
   m_earth   = 5.97219e24
   V_earth   = 4/3.0 * pi * r_earth**3
   rho_earth = m_earth / V_earth

   #addParticle(x, y, z, vx, vy, vz, r, rho,
   #            thetax, thetay, thetaz, 
   #            omegax, omegay, omegaz): 
   p.addParticle(0,0,r_earth + 1e2,0,0,0,3,0,0,0,0,0,0,0)
   p.addParticle(0, 0, 0, 0, 0, 0, r_earth, rho_earth, 0, 0, 0, 0, 0, 0)


def initialize_planet(p):
   """
   """
   G         = p.f.G
   r_earth   = 10000
   rho_earth = 10000000
   rho_moon  = 10000000
   r_moon    = 750
   d_moon    = r_earth + r_moon + 10000
   
   V_earth   = 4.0/3.0 * pi * r_earth**3
   V_moon    = 4.0/3.0 * pi * r_moon**3
   m_earth   = rho_earth * V_earth
   m_moon    = rho_moon * V_moon
   vt_moon   = sqrt(G * m_earth / d_moon)
   vt_earth  = sqrt(G * m_moon / d_moon)
   print vt_moon

   #addParticle(x, y, z, vx, vy, vz, r, rho,
   #            thetax, thetay, thetaz,
   #            omegax, omegay, omegaz):
   p.addParticle(0,0, 30000, 0,0,0,3,1e-16,0,0,0,0,0,0)
   p.addParticle(0, 0, 0, 0, 0, 0, r_earth, rho_earth, 
                 0, 0, 0, 0, 0, 0)
   p.addParticle(0, 0, d_moon, vt_moon, 0, 0, r_moon, rho_moon, 
                 0, 0, 0, 0, 0, 0)

def initialize_nebula(p, n, r, L):
  """
  """
  dx = 2.0*L / n
  d  = linspace(dx/2.0 - L, L - dx/2.0, n)

  p.addParticle(0,0,-2*L,0,0,0,3)
  for i in d:
    for j in d:
      for k in d:
        p.addParticle(i,j,k,0,0,0,r)



