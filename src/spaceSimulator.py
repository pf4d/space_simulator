from particles      import *
from OpenGL.GL      import *
from OpenGL.GLUT    import *
from OpenGL.GLE     import *
from OpenGL.GLU     import *
from FTGL           import *
from pylab          import *
from objLoader      import *
from time           import time
import sys

class Camera(object):

  def __init__(self):
    """
    """
    self.initialize()

  def initialize(self):
    """
    reset the rotation to the identity matrix.
    """
    self.M = identity(3)
    self.M[2,2] = -1
    #self.update([.2, 0, 0])  # provide an initial rotation up

  def update(self, v):
    """
    rotate the camera orientation matrix about the x, y, and z axes by angles 
    provided in <v> array.
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
    
    self.M = dot(self.M, R)

fwd       = False          # craft moving forward
back      = False          #   "     "    back
left      = False          #   "     "    left
right     = False          #   "     "    right
up        = False          #   "     "    up
down      = False          #   "     "    down
rollLeft  = False          #   "   rolling left 
rollRight = False          #   "     "     right
yawLeft   = False          #   "   rotating left
yawRight  = False          #   "   rotating right
ascend    = False          #   "   ascending
descend   = False          #   "   descending
paused    = False          # pauses the game
                           
                           
taccel    = 20.0           # translational acceleration to apply 
raccel    = 1.0            # rotational acceleration to apply
                           
rotx      = 0              # camera x rotation
roty      = 0              # camera y rotation
rotz      = 0              # camera z rotation
camDist   = 4.0            # camera distance coef.

# create viewpoint camera :
camera    = Camera()

# colors used :
shipColor  = 0.5*ones(3)             # color of the ship
axesColor  = array([0.2, 0.2, 0.2])  # color of the readout axes
aColor     = array([0.0, 0.0, 1.0])  # acceleration vector color
vColor     = array([0.0, 1.0, 0.0])  # velocity vector color
alphaColor = array([1.0, 0.0, 0.0])  # angular acceleration vector color
omegaColor = array([1.0, 1.0, 0.0])  # angular velocity vector color

# variable used :
frames     = 0                       # for spf calculation
lastTime   = time()                  # current time
fps        = 1.0                     # current frames per second
w          = 700                     # screen width
h          = 700                     # screen height
                                     
r_earth    = 6.3781e6
au         = 1.49597871e11
dt         = 0.1                     # time step
L          = 120.0                   # size of the box
t          = 0                       # initial time
                                    
## "pool balls" :                    
#k          = 30.0                    # elastic 'bounce'
#gamma      = 0.1                     # energy dissipation/loss
## "squishy balls" :                  
#k          = 1.5                     # elastic 'bounce'
#gamma      = 0.1                     # energy dissipation/loss
# "space balls" :                   
k          = 30.0                    # elastic 'bounce'
gamma      = 0.4                     # energy dissipation/loss
                                    
rho        = 1e9                     # particle denisty
g          = 0.0                     # downward acceleration

# particle update data:
COUNT         = 1          # number of time steps computed
UPDATE_FRAMES = 1          # how often to redraw screen

# how resolved are the spheres?
STACKS = 10
SLICES = 15

# create specter field :
specter = Specter(1e3, 4*L)

# create star field :
#star    = Specter(1e3, 10*r_earth)
star    = Specter(1e3, 1000*L)

# instantiate the forces function between particles
f = GranularMaterialForce(k=k, gamma=gamma)
#f = NebulaGranularMaterialForce(k=k, gamma=gamma)

# create some particles and a box
p = Particles(L, f, periodicY=0, periodicZ=0, periodicX=0)
#p = Nebula(L, f, periodicY=0, periodicZ=0, periodicX=0)

#  addParticle(x, y, z, vx, vy, vz, r, rho,
#              thetax, thetay, thetaz, 
#              omegax, omegay, omegaz): 
#initialize_nebula(p, 5, 8.0, L/2)
k          = 30.0                    # elastic 'bounce'
gamma      = 0.1                     # energy dissipation/loss
initialize_grid(p, 4, 3.0, rho, L)
#initialize_random(p, 100, 3.0, rho, L/2)
#initialize_system(p)
#initialize_earth(p)
#initialize_planet(p)

# instantiate Integrator
integrate = VerletIntegrator(dt)
#integrate = NebulaVerletIntegrator(dt)

def init():
  """
  """
  # general properties :
  glEnable(GL_COLOR_MATERIAL)
  glEnable(GL_BLEND)
  glShadeModel(GL_SMOOTH)
  #glShadeModel(GL_FLAT)
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE)

  glEnable(GL_POLYGON_OFFSET_FILL) # Prevents hidden line problems when drawing
  glPolygonOffset(1.0, 1.0)        # a wireframe on top of filled polygons.

  # 3d parameters :
  glEnable(GL_CULL_FACE)
  glEnable(GL_DEPTH_TEST)
 
  # lights : 
  glEnable(GL_LIGHTING)
  glEnable(GL_LIGHT0)
  glEnable(GL_LIGHT1)

  # point config :
  glEnable(GL_POINT_SPRITE)
  glEnable(GL_POINT_SMOOTH)
  glEnable(GL_BLEND)
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
  glPointSize(10.0)


def draw_ship_vectors(dx, dy):
  """
  """
  global shipColor, axesColor, aColor, vColor, alphaColor, omegaColor

  w = float(glutGet(GLUT_WINDOW_WIDTH))
  h = float(glutGet(GLUT_WINDOW_HEIGHT))
  ar = w / h

  # save projection matrix
  glMatrixMode(GL_PROJECTION)
  glPushMatrix()

  # switch to orthographic projection :
  glLoadIdentity()
  glOrtho(-L*ar, L*ar, -L, L, -4*L, 4*L)

  # back to the modelview matrix mode, so that we can translate/scale text :
  glMatrixMode(GL_MODELVIEW)

  # draw the ship statistics :
  glMaterial(GL_FRONT, GL_EMISSION,  [0.0, 0.0, 0.0, 0.0])
  glMaterial(GL_FRONT, GL_SPECULAR,  [0.5, 0.5, 0.5, 0.0])
  glMaterial(GL_FRONT, GL_SHININESS, 100.0)
  glDisable(GL_LIGHTING)   # disable lighting
  
  # parameters for positioning :
  xt =  ar*L-1.5*dx
  yt = -L+2*dy
  zt =  L
  xr =  ar*L-3.5*dx
  yr =  yt
  zr =  zt
  
  # font parameters :
  font = BitmapFont('ProggyTinySZ.ttf')
  font.FaceSize(15)
  
  # vectors of orientation :
  R      = rotate_vector(array([-pi/6, pi/6, 0]))
  M      = dot(p.theta[0].T, R)
  av     = array([p.ax[0], p.ay[0], p.az[0]])
  vv     = array([p.vx[0], p.vy[0], p.vz[0]]) + 1e-16
  alphav = -array([p.alphax[0], p.alphay[0], p.alphaz[0]])
  omegav = -array([p.omegax[0], p.omegay[0], p.omegaz[0]]) + 1e-16
  x      = array([1,0,0])
  y      = array([0,1,0])
  z      = array([0,0,1])
  xyz1   = zeros(3)
  ca     = 15.0
  xyz2_x = xyz1 + ca*x
  xyz2_y = xyz1 + ca*y
  xyz2_z = xyz1 + ca*z
  
  #=============================================================================
  # draw the translational info :
  glPushMatrix()
  glLoadIdentity()
  glTranslate(xt, yt, zt)
  mvm = glGetFloatv(GL_MODELVIEW_MATRIX)   # current rotation / translation
  mvm[:3,:3] = R                           # rotate the view by R
  vec = dot(p.theta[0], vv)                # find components of velocity
  glLoadMatrixf(mvm)                       # reload the modelView matrix
  glColor(shipColor)
  glCallList(obj.gl_list)

  glColor(vColor)
  glRasterPos3f(sign(vec[0])*xyz2_x[0], xyz2_x[1], xyz2_x[2])
  font.Render("%.1f" % vec[0])
  glRasterPos3f(xyz2_y[0], sign(vec[1])*xyz2_y[1], xyz2_y[2])
  font.Render("%.1f" % vec[1])
  glRasterPos3f(xyz2_z[0], xyz2_z[1], sign(vec[2])*xyz2_z[2])
  font.Render("%.1f" % -vec[2])
  
  glLineWidth(1.0)
  glColor(axesColor)
  glBegin(GL_LINES)
  glVertex3fv(xyz1 - xyz2_x)
  glVertex3fv(xyz2_x)
  glVertex3fv(xyz1 - xyz2_y)
  glVertex3fv(xyz2_y)
  glVertex3fv(xyz1 - xyz2_z)
  glVertex3fv(xyz2_z)
  glEnd()
  glPopMatrix()
 
  # draw the rotational info : 
  glPushMatrix()
  glLoadIdentity()
  glTranslate(xr, yr, zr)
  mvm = glGetFloatv(GL_MODELVIEW_MATRIX)   # current rotation / translation
  mvm[:3,:3] = R                           # rotate the view by R
  glLoadMatrixf(mvm)                       # reload the modelView matrix
  glColor(shipColor)
  glCallList(obj.gl_list)
  
  glColor(omegaColor)
  glRasterPos3f(sign(omegav[0])*xyz2_x[0], xyz2_x[1], xyz2_x[2])
  font.Render("%.1f" % omegav[0])
  glRasterPos3f(xyz2_y[0], sign(omegav[1])*xyz2_y[1], xyz2_y[2])
  font.Render("%.1f" % omegav[1])
  glRasterPos3f(xyz2_z[0], xyz2_z[1], sign(omegav[2])*xyz2_z[2])
  font.Render("%.1f" % -omegav[2])
  
  glLineWidth(1.0)
  glColor(axesColor)
  glBegin(GL_LINES)
  glVertex3fv(xyz1 - xyz2_x)
  glVertex3fv(xyz2_x)
  glVertex3fv(xyz1 - xyz2_y)
  glVertex3fv(xyz2_y)
  glVertex3fv(xyz1 - xyz2_z)
  glVertex3fv(xyz2_z)
  glEnd()
  glPopMatrix()

  #=============================================================================
  # translational statistics :
  glPushMatrix()
  glLoadIdentity()
  glTranslate(xt, yt, zt)
  mvm = glGetFloatv(GL_MODELVIEW_MATRIX)   # current rotation / translation
  mvm[:3,:3]  = M                          # rotate the view by M
  glLoadMatrixf(mvm)                       # reload the modelView matrix
  glLineWidth(2.0)
  glBegin(GL_LINES)
 
  # acceleration vector :
  glColor(aColor)
  c    = ca/taccel
  xyz2 = xyz1 + c*av
  glVertex3fv(xyz1)
  glVertex3fv(xyz2)
  
  # velocity vector :
  glColor(vColor)
  c    = 1.0
  xyz2 = xyz1 + c*vv
  glVertex3fv(xyz1)
  glVertex3fv(xyz2)
  
  glEnd()
  glPopMatrix()

  #=============================================================================
  # rotational statistics : 
  glPushMatrix()
  glLoadIdentity()
  glTranslate(xr, yr, zr)
  mvm = glGetFloatv(GL_MODELVIEW_MATRIX)   # current rotation / translation
  mvm[:3,:3]  = M                          # rotate the view by M
  glLoadMatrixf(mvm)                       # reload the modelView matrix
  glBegin(GL_LINES)
 
  # angular acceleration vector :
  glColor(alphaColor)
  c    = ca/raccel
  xyz2 = xyz1 + c*alphav 
  glVertex3fv(xyz1)
  glVertex3fv(xyz2)
  
  # angular velocity vector :
  glColor(omegaColor)
  c    = 20.0
  xyz2 = xyz1 + c*omegav
  glVertex3fv(xyz1)
  glVertex3fv(xyz2)
  
  # done :
  glEnd()
  glPopMatrix()
  
  # re-enable lighting :
  glEnable(GL_LIGHTING)

  # get back to old perspective matrix :
  glMatrixMode(GL_PROJECTION)
  glPopMatrix()
  glMatrixMode(GL_MODELVIEW)


def print_ship_stats(dx, dy):
  """
  """
  # save projection matrix
  glMatrixMode(GL_PROJECTION)
  glPushMatrix()

  # switch to orthographic projection :
  glLoadIdentity()
  glOrtho(-L, L, -L, L, -4*L, 4*L)

  # back to the modelview matrix mode, so that we can translate/scale text :
  glMatrixMode(GL_MODELVIEW)

  glDisable(GL_LIGHTING)   # disable lighting
  
  # print statistics :
  glPushMatrix()
  glLoadIdentity()
  
  glColor(1.0,1.0,1.0,1.0) 
  font = BitmapFont('ProggySquareSZ.ttf')
  font.FaceSize(16)
  #font = TextureFont('ProggySquareSZ.ttf')
  #font.FaceSize(13)
  #glScale(0.05, 0.05, 0.05)
  t1 = 'ship statistics (x,y,z) :'
  t2 = 'theta        : %.2E, %.2E, %.2E' % (p.thetax[0],p.thetay[0],p.thetaz[0])
  t3 = 'omega        : %.2E, %.2E, %.2E' % (p.omegax[0],p.omegay[0],p.omegaz[0])
  t4 = 'alpha        : %.2E, %.2E, %.2E' % (p.alphax[0],p.alphay[0],p.alphaz[0])
  t5 = 'position     : %.2E, %.2E, %.2E' % (p.x[0],  p.y[0],  p.z[0])
  t6 = 'velocity     : %.2E, %.2E, %.2E' % (p.vx[0], p.vy[0], p.vz[0])
  t7 = 'acceleration : %.2E, %.2E, %.2E' % (p.ax[0], p.ay[0], p.az[0])
  glRasterPos2f(-L+dy,-(L-5)+dy*3)
  font.Render(t1)
  glRasterPos2f(-L+dy,-(L-5)+dy*2.5)
  font.Render(t2)
  glRasterPos2f(-L+dy,-(L-5)+dy*2)
  font.Render(t3)
  glRasterPos2f(-L+dy,-(L-5)+dy*1.5)
  font.Render(t4)
  glRasterPos2f(-L+dy,-(L-5)+dy*1.0)
  font.Render(t5)
  glRasterPos2f(-L+dy,-(L-5)+dy*0.5)
  font.Render(t6)
  glRasterPos2f(-L+dy,-(L-5)+dy*0.0)
  font.Render(t7)
  glPopMatrix()
  
  # re-enable lighting :
  glEnable(GL_LIGHTING)

  # get back to old perspective matrix :
  glMatrixMode(GL_PROJECTION)
  glPopMatrix()
  glMatrixMode(GL_MODELVIEW)

def print_paused(dx,dy):
  """
  """
  # save projection matrix
  glMatrixMode(GL_PROJECTION)
  glPushMatrix()

  # switch to orthographic projection :
  glLoadIdentity()
  glOrtho(-L, L, -L, L, -4*L, 4*L)

  # back to the modelview matrix mode, so that we can translate/scale text :
  glMatrixMode(GL_MODELVIEW)

  glDisable(GL_LIGHTING)   # disable lighting
  
  # print statistics :
  glPushMatrix()
  glLoadIdentity()
  
  glColor(1.0,1.0,1.0,1.0) 
  #glRasterPos2f(L-dx-10, L-dy)
  font = BitmapFont('ProggySquareSZ.ttf')
  font.FaceSize(30)
  #font = TextureFont('ProggySquareSZ.ttf')
  #font.FaceSize(13)
  #glScale(0.05, 0.05, 0.05)
  #font.Render("n = %i" % p.N)
  
  glRasterPos2f(-L+dy + 20, L-dy-20)
  font.Render("Forward:        w")
  glRasterPos2f(-L+dy+20, L-dy-30)
  font.Render("Back:           s")
  glRasterPos2f(-L+dy+20, L-dy-40)
  font.Render("Left:           d")
  glRasterPos2f(-L+dy+20, L-dy-50)
  font.Render("Right:          a")
  glRasterPos2f(-L+dy+20, L-dy-60)
  font.Render("Yaw Left:       q")
  glRasterPos2f(-L+dy+20, L-dy-70)
  font.Render("Yaw Right:      e")
  glRasterPos2f(-L+dy+20, L-dy-80)
  font.Render("Pitch Up:       Up Arrow")
  glRasterPos2f(-L+dy+20, L-dy-90)
  font.Render("Pitch Down:     Down Arrow")
  glRasterPos2f(-L+dy+20, L-dy-100)
  font.Render("Roll Left:      Left Arrow")
  glRasterPos2f(-L+dy+20, L-dy-110)
  font.Render("Roll Right:     Right Arrow")
  glRasterPos2f(-L+dy+20, L-dy-120)
  font.Render("Descend:        Page Down")
  glRasterPos2f(-L+dy+20, L-dy-130)
  font.Render("Ascend:         Page Up")
  glRasterPos2f(-L+dy+20, L-dy-140)
  font.Render("Pause/Resume:   p")
  glRasterPos2f(-L+dy+20, L-dy-150)
  font.Render("Rotate Camera:  Left-click")
  glRasterPos2f(-L+dy+20, L-dy-160)
  font.Render("Zoom In/Out:    Scroll Wheel")
  glRasterPos2f(-L+dy+20, L-dy-170)
  font.Render("Exit:           esc")
  glPopMatrix()
  
  # re-enable lighting :
  glEnable(GL_LIGHTING)

  # get back to old perspective matrix :
  glMatrixMode(GL_PROJECTION)
  glPopMatrix()
  glMatrixMode(GL_MODELVIEW)

def print_stats(dx, dy):
  """
  """
  # save projection matrix
  glMatrixMode(GL_PROJECTION)
  glPushMatrix()

  # switch to orthographic projection :
  glLoadIdentity()
  glOrtho(-L, L, -L, L, -4*L, 4*L)

  # back to the modelview matrix mode, so that we can translate/scale text :
  glMatrixMode(GL_MODELVIEW)

  glDisable(GL_LIGHTING)   # disable lighting
  
  # print statistics :
  glPushMatrix()
  glLoadIdentity()
  
  glColor(1.0,1.0,1.0,1.0) 
  glRasterPos2f(L-dx, L-dy)
  font = BitmapFont('ProggySquareSZ.ttf')
  font.FaceSize(16)
  #font = TextureFont('ProggySquareSZ.ttf')
  #font.FaceSize(13)
  #glScale(0.05, 0.05, 0.05)
  font.Render("n = %i" % p.N)
  glRasterPos2f(-L+dy, L-dy)
  font.Render("%i FPS" % fps)
  glRasterPos2f(-L+dy, -L+dy)
  font.Render("Press p to pause/view controls")
  glPopMatrix()
  
  # re-enable lighting :
  glEnable(GL_LIGHTING)

  # get back to old perspective matrix :
  glMatrixMode(GL_PROJECTION)
  glPopMatrix()
  glMatrixMode(GL_MODELVIEW)


def draw_velocity_vectors():
  """
  draw velocity vectors.
  """
  global vColor

  glLineWidth(1.0)
  glPopMatrix()  
  glColor(vColor)
  glDisable(GL_LIGHTING)
  glBegin(GL_LINES)
  for i in range(p.N):
    v_mag = sqrt(p.vx[i]**2 + p.vy[i]**2 + p.vz[i]**2) + 1e-16
    xyz1 = array([p.x[i],  p.y[i],  p.z[i]])
    vxyz = array([p.vx[i], p.vy[i], p.vz[i]])
    vxyz = vxyz / v_mag * 2.0*p.r[i]
    xyz2 = xyz1 + vxyz
    glVertex3fv(xyz1)
    glVertex3fv(xyz2)
  glEnd()
  glEnable(GL_LIGHTING)
  glPushMatrix()


def draw_acceleration_vectors():
  """
  draw acceleration vectors.
  """
  global aColor

  glLineWidth(1.0)
  glPopMatrix()  
  glColor(aColor)
  glDisable(GL_LIGHTING)
  glBegin(GL_LINES)
  for i in range(p.N):
    a_mag = sqrt(p.ax[i]**2 + p.ay[i]**2 + p.az[i]**2) + 1e-16
    xyz1 = array([p.x[i],  p.y[i],  p.z[i]])
    axyz = array([p.ax[i], p.ay[i], p.az[i]])
    axyz = axyz / a_mag * 2.0*p.r[i]
    xyz2 = xyz1 + axyz
    glVertex3fv(xyz1)
    glVertex3fv(xyz2)
  glEnd()
  glEnable(GL_LIGHTING)
  glPushMatrix()


def draw_rotation_vectors():
  """
  draw rotation vectors.
  """
  # FIXME: broken
  glLineWidth(1.0)
  glPopMatrix()  
  glColor4f(0.0,1.0,1.0,1.0)
  glDisable(GL_LIGHTING)
  glBegin(GL_LINES)
  for i in range(p.N):
    xyz1  = array([p.x[i],  p.y[i],  p.z[i]])
    fxyz  = array([0, 0, 1])            # initial forward vector
    uxyz  = array([0, 1, 0])            # initial up vector
    txyz  = array([p.thetax[i], p.thetay[i], p.thetaz[i]])
    fxyz  = rotate_vector(fxyz, txyz)
    uxyz  = rotate_vector(uxyz, txyz)
    fxyz  = fxyz / norm(fxyz) * 2.0*p.r[i]
    uxyz  = uxyz / norm(uxyz) * 2.0*p.r[i]
    fxyz2 = xyz1 + fxyz
    uxyz2 = xyz1 + uxyz
    glVertex3fv(xyz1)
    glVertex3fv(fxyz2)
    glVertex3fv(xyz1)
    glVertex3fv(uxyz2)
  glEnd()
  glEnable(GL_LIGHTING)
  glPushMatrix()


def draw_angular_velocity_vectors():  
  """
  draw angular velocity vectors.
  """
  global omegaColor
  
  glLineWidth(1.0)
  glPopMatrix()
  glColor(omegaColor)
  glDisable(GL_LIGHTING)
  glBegin(GL_LINES)
  for i in range(p.N):
    omega_mag = sqrt(p.omegax[i]**2 + p.omegay[i]**2 + p.omegaz[i]**2) + 1e-16
    xyz1 = array([p.x[i],      p.y[i],      p.z[i]])
    vxyz = -array([p.omegax[i], p.omegay[i], p.omegaz[i]]) 
    vxyz = vxyz / omega_mag * (p.r[i]+0.5)
    xyz2 = xyz1 + vxyz
    glVertex3fv(xyz1)
    glVertex3fv(xyz2)
  glEnd()
  glEnable(GL_LIGHTING)
  glPushMatrix()

  
def draw_angular_acceleration_vectors():
  """
  draw angular acceleration vectors.
  """ 
  global alphaColor

  glLineWidth(1.0)
  glPopMatrix()
  glColor(alphaColor)
  glDisable(GL_LIGHTING)
  glBegin(GL_LINES)
  for i in range(p.N):
    alpha_mag = sqrt(p.alphax[i]**2 + p.alphay[i]**2 + p.alphaz[i]**2) + 1e-16
    xyz1 = array([p.x[i],      p.y[i],      p.z[i]])
    vxyz = -array([p.alphax[i], p.alphay[i], p.alphaz[i]])
    vxyz = vxyz / alpha_mag * (p.r[i]+0.5)
    xyz2 = xyz1 + vxyz
    glVertex3fv(xyz1)
    glVertex3fv(xyz2)
  glEnd()
  glEnable(GL_LIGHTING)
  glPushMatrix()


def rotate_vector(r):
  """
  rotate vector <v> about the x, y, and z axes by angles provided in <r> array.
  """
  rx = r[0]
  ry = r[1]
  rz = r[2]
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
  return R


def draw_specter_field(S, r):
  """
  """
  glDisable(GL_LIGHTING)   # disable lighting
  glPointSize(r)
  for i in range(S.n):
    mag = sqrt(S.x[i]**2 + S.y[i]**2 + S.z[i]**2)
    tra = (S.L - mag)/S.L
    glColor(1.0,1.0,1.0,tra)
    glBegin(GL_POINTS)
    glVertex3f(S.x[i], S.y[i], S.z[i])
    glEnd()
  glEnable(GL_LIGHTING)    # enable lighting

def draw_solid_particles():
  """
  draw the spheres.
  """
  glPushMatrix() 
  glMaterial(GL_FRONT, GL_EMISSION,  [0.0, 0.0, 0.0, 0.0])
  glMaterial(GL_FRONT, GL_SPECULAR,  [0.5, 0.5, 0.5, 0.0])
  glMaterial(GL_FRONT, GL_SHININESS, 100.0)
  for i in range(p.N):
    a = array([p.ax[i], p.ay[i], p.az[i]])
    if norm(a) > 20000.0 and i != 0:
      glColor(1/2.0, 1/2.0, 1/2.0, 1.0)
    elif i != 0:
      glColor(1, 1/2.0, 0.0, 1.0)
    
    glPushMatrix()
    glTranslate(p.x[i], p.y[i], p.z[i])
    
    # rotation :
    mvm = glGetFloatv(GL_MODELVIEW_MATRIX)    # current rotation / translation
    mvm[:3,:3] = dot(p.theta[i], mvm[:3,:3])  # rotate the view by theta[i]
    glLoadMatrixf(mvm)                        # reload the modelView matrix
    
    # draw particles as spheres or the ship if index == 0 :
    if i == 0 and camDist != 0:
      glColor(shipColor)
      glCallList(obj.gl_list)
    else:  
      glutSolidSphere(p.r[i], SLICES, STACKS)
    
    glPopMatrix()

def draw_wireframe_particles():
  """
  draw the spheres.
  """
  glPushMatrix() 
  glMaterial(GL_FRONT, GL_SPECULAR,  [0.0, 0.0, 0.0, 0.0])
  glMaterial(GL_FRONT, GL_SHININESS, 0.0)
  glMaterial(GL_FRONT, GL_EMISSION,  [0.0, 0.0, 0.0, 0.0])
  glLineWidth(1.0)
  for i in range(p.N):
    glPushMatrix()
    glTranslate(p.x[i], p.y[i], p.z[i])
    
    # rotation :
    mvm = glGetFloatv(GL_MODELVIEW_MATRIX)    # current rotation / translation
    mvm[:3,:3] = dot(p.theta[i], mvm[:3,:3])  # rotate the view by theta[i]
    glLoadMatrixf(mvm)                        # reload the modelView matrix
    
    # draw particles as spheres or the ship if index == 0 :
    if i == 0 and camDist != 0:
      glColor(shipColor)
      glCallList(obj.gl_list)
    elif i == 0:
      pass
      #glDisable(GL_LIGHTING)
      #glColor(1.0,0.0,0.0,1.0)
      #glBegin(GL_POINTS)
      #glVertex3f(p.x[i], p.y[i], -p.z[i])
      #glEnd()
      #glEnable(GL_LIGHTING)
    else:
      glColor(0.0,1.0,0.0,1.0)
      glutWireSphere(p.r[i]*1.01, 10, 5)
    
    glPopMatrix()

def draw_point_particles():
  """
  draw the spheres.
  """
  glPushMatrix() 
  glMaterial(GL_FRONT, GL_SPECULAR,  [0.0, 0.0, 0.0, 0.0])
  glMaterial(GL_FRONT, GL_SHININESS, 0.0)
  glMaterial(GL_FRONT, GL_EMISSION,  [0.0, 0.0, 0.0, 0.0])
  for i in range(p.N):
    glPushMatrix()
    glTranslate(p.x[i], p.y[i], p.z[i])
    
    # draw particles as spheres or the ship if index == 0 :
    if i == 0 and camDist != 0:
      glColor(shipColor)
      glCallList(obj.gl_list)
    else:  
      glColor(0.0,0.0,1.0,1.0)
      glBegin(GL_POINTS)
      glVertex3f(p.x[i], p.y[i], p.z[i])
      glEnd()

    glPopMatrix()

def display():
  """
  """
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

  dx = 0.2 * L
  dy = 0.1 * L
  
  if paused:
    print_paused(dx,dy)
  # camera viewpoint :
  glLoadIdentity()
  M  = dot(p.theta[0].T, camera.M)
  pt = array([p.x[0], p.y[0], p.z[0]]) + 1e-15
  rt = M[:,0]
  #up = cross(M[:,2], rt)
  up = M[:,1]
  fr = M[:,2] + pt
  pt = pt - M[:,2] / norm(M[:,2]) * camDist * p.r[0]  # moves the camera back
  gluLookAt(pt[0], pt[1], pt[2],              # Camera Position
            fr[0], fr[1], fr[2],              # Point the Camera looks at
            up[0], up[1], up[2])              # the Up-Vector
  
  #mvm = glGetFloatv(GL_MODELVIEW_MATRIX)
  #mvm[:3,:3] = dot(p.theta[0], mvm[:3,:3])
  #glLoadMatrixf(mvm)
  
  ## use fixed camera at center of domain :
  #gluLookAt(0,0,L*camDist,                   # Camera Position
  #          0,0,0,                           # Point the Camera looks at
  #          0,1,0)                           # the Up-Vector
  #glRotate(rotx,1,0,0)
  #glRotate(roty,0,1,0)
  #glRotate(rotz,0,0,1)

  # draw specter field :
  #draw_specter_field(specter, 1.0)
 
  # draw star field :
  draw_specter_field(star, 4.0)
 
  draw_solid_particles()
  #draw_wireframe_particles()
  #draw_point_particles()

  # draw vectors on particles :
  #draw_velocity_vectors()
  #draw_acceleration_vectors()
  #draw_rotation_vectors()
  #draw_angular_velocity_vectors()
  #draw_angular_acceleration_vectors()   

  # print statistics :
  print_stats(dx,dy)

  # print ship stats :
  #print_ship_stats(dx,dy)

  # draw ship vectors :
  draw_ship_vectors(dx,dy)
 
  # draw the lights : 
  lx1 = 0.0
  ly1 = 30000
  lz1 = 0.0
  
  lx2 = L
  ly2 = 30000
  lz2 = L
 
  glMaterial(GL_FRONT, GL_EMISSION,  [1.0, 1.0, 1.0, 0.0])
  glTranslate(lx1, ly1, lz1)
  glutSolidSphere(0.5, SLICES, STACKS)
  glTranslate(lx2, ly2, lz2)
  glutSolidSphere(0.5, SLICES, STACKS)
  glPopMatrix()

  # spot light :
  glLight(GL_LIGHT0, GL_AMBIENT,               [0.2,  0.2, 0.2, 1.0])
  glLight(GL_LIGHT0, GL_SPECULAR,              [1.0,  1.0, 1.0, 1.0])
  glLight(GL_LIGHT0, GL_DIFFUSE,               [1.0,  1.0, 1.0, 1.0])
  glLight(GL_LIGHT0, GL_POSITION,              [lx1,  ly1, lz1, 1.0])
  glLight(GL_LIGHT0, GL_SPOT_DIRECTION,        [0.0, -1.0, 0.0])
  glLight(GL_LIGHT0, GL_SPOT_CUTOFF,           10.0)
  glLight(GL_LIGHT0, GL_SPOT_EXPONENT,          3.0)
  glLight(GL_LIGHT0, GL_CONSTANT_ATTENUATION,   1.2)
  glLight(GL_LIGHT0, GL_LINEAR_ATTENUATION,     0.0)
  glLight(GL_LIGHT0, GL_QUADRATIC_ATTENUATION,  0.002)
  
  # regular light
  glLight(GL_LIGHT1, GL_AMBIENT,               [0.2,  0.2, 0.2, 1.0])
  glLight(GL_LIGHT1, GL_SPECULAR,              [1.0,  1.0, 1.0, 1.0])
  glLight(GL_LIGHT1, GL_DIFFUSE,               [1.0,  1.0, 1.0, 1.0])
  glLight(GL_LIGHT1, GL_POSITION,              [lx2,  ly2, lz2, 1.0])
  glLight(GL_LIGHT1, GL_CONSTANT_ATTENUATION,   2.0)
  glLight(GL_LIGHT1, GL_LINEAR_ATTENUATION,     0.0)
  glLight(GL_LIGHT1, GL_QUADRATIC_ATTENUATION,  0.0)
  
  glutSwapBuffers()
  glFlush()


def reshape(width, height):
  """
  """
  glViewport(0, 0, width, height)
  glMatrixMode(GL_PROJECTION)
  glLoadIdentity()
  ar = width/float(height)
  gluPerspective(90.0, ar, 1, 10*au)
  #glOrtho(-L*ar, L*ar, -L, L, -4*L, 4*L)
  glMatrixMode(GL_MODELVIEW)
  glLoadIdentity()


def idle():
  """
  """
  global COUNT, frames, lastTime, fps, yawLeft, yawRight, ascend, descend
  global fwd, back, left, right, up, down, taccel, raccel, rollLeft, rollRight
  global paused
  
  if not paused:
    # integrate the system forward in time :
    for i in range(UPDATE_FRAMES):
      integrate(f,p)
      COUNT = COUNT + 1 
  
    # move forward :
    if fwd == True:
      vf = p.theta[0][2,:]
      pf = vf / norm(vf) * taccel
      p.ax[0] -= pf[0]
      p.ay[0] -= pf[1]
      p.az[0] -= pf[2]
  
    # move backward :
    elif back == True:
      vf = p.theta[0][2,:]
      pf = vf / norm(vf) * taccel
      p.ax[0] += pf[0]
      p.ay[0] += pf[1]
      p.az[0] += pf[2]
  
    # move left :
    if left == True:
      vf = p.theta[0][0,:]
      pf = vf / norm(vf) * taccel
      p.ax[0] -= pf[0]
      p.ay[0] -= pf[1]
      p.az[0] -= pf[2]
  
    # move right :
    elif right == True:
      vf = p.theta[0][0,:]
      pf = vf / norm(vf) * taccel
      p.ax[0] += pf[0]
      p.ay[0] += pf[1]
      p.az[0] += pf[2]

    # pitch up :
    if up == True:
      vf = p.theta[0][0,:]
      pf = vf / norm(vf) * raccel
      p.alphax[0] += pf[0]
      p.alphay[0] += pf[1]
      p.alphaz[0] += pf[2]
  
    # pitch down :
    elif down == True:
      vf = p.theta[0][0,:]
      pf = vf / norm(vf) * raccel
      p.alphax[0] -= pf[0]
      p.alphay[0] -= pf[1]
      p.alphaz[0] -= pf[2]

    # roll left :
    if rollLeft == True:
      vf = p.theta[0][2,:]
      pf = vf / norm(vf) * raccel
      p.alphax[0] -= pf[0]
      p.alphay[0] -= pf[1]
      p.alphaz[0] -= pf[2]
  
    # roll right :
    elif rollRight == True:
      vf = p.theta[0][2,:]
      pf = vf / norm(vf) * raccel
      p.alphax[0] += pf[0]
      p.alphay[0] += pf[1]
      p.alphaz[0] += pf[2]

    # yaw left :
    if yawLeft == True:
      vf = p.theta[0][1,:]
      pf = vf / norm(vf) * raccel
      p.alphax[0] -= pf[0]
      p.alphay[0] -= pf[1]
      p.alphaz[0] -= pf[2]
  
    # yaw right :
    elif yawRight == True:
      vf = p.theta[0][1,:]
      pf = vf / norm(vf) * raccel
      p.alphax[0] += pf[0]
      p.alphay[0] += pf[1]
      p.alphaz[0] += pf[2]

    # ascend :
    if ascend == True:
      vf = p.theta[0][1,:]
      pf = vf / norm(vf) * taccel
      p.ax[0] += pf[0]
      p.ay[0] += pf[1]
      p.az[0] += pf[2]
  
    # descend :
    elif descend == True:
      vf = p.theta[0][1,:]
      pf = vf / norm(vf) * taccel
      p.ax[0] -= pf[0]
      p.ay[0] -= pf[1]
      p.az[0] -= pf[2]

  # redraw the screen :
  glutPostRedisplay()

  # calculate fps :
  currentTime = time()
  frames     += 1
  if (currentTime - lastTime) >= 1.0:
    fps       = frames
    frames    = 0
    lastTime += 1.0
    #print fps


def key(k, x, y):
  """
  """
  global camera, fwd, back, left, right, yawLeft, yawRight, paused
  
  # reset the camera :
  if k == 'c':
    print "'c' was pressed, reseting camera"
    camera.initialize()
  
  # quit the game :
  if k == '\033':
    print "'esc' was pressed"
    exit(0)

  if k == 'p':
    print "'p' was pressed"
    if paused:
      paused = False
    else:
      paused = True
  
  # print the number of particles :
  if k == 'n':
    print "'n' was pressed: n =", p.N
  
  # move forward :
  if k == 'w':
    fwd = True
  
  # move backward :
  if k == 's':
    back = True

  # move left :
  if k == 'a':
    left = True

  # move right :
  if k == 'd':
    right = True

  # yaw right :
  if k == 'e':
    yawRight = True

  # yaw left :
  if k == 'q':
    yawLeft = True


def keyUp(k,x,y):
  """
  """
  global fwd, back, left, right, yawLeft, yawRight
  
  # stop moving forward :
  if k == 'w':
    fwd = False

  # stop moving backward :
  if k == 's':
    back = False

  # stop moving left :
  if k == 'a':
    left = False

  # stop moving right :
  if k == 'd':
    right = False

  # stop yaw left :
  if k == 'q':
    yawLeft = False

  # stop yaw right :
  if k == 'e':
    yawRight = False
    

def special(k, x, y):
  """
  """
  global up, down, rollLeft, rollRight, ascend, descend
  
  # pitch up :
  if k == GLUT_KEY_UP:
    up = True
  
  # pitch down :
  if k == GLUT_KEY_DOWN:
    down = True
  
  # roll right :
  if k == GLUT_KEY_RIGHT:
    rollRight = True
  
  # roll left :
  if k == GLUT_KEY_LEFT:
    rollLeft = True

  # descend :
  if k == GLUT_KEY_PAGE_DOWN:
    descend = True
  
  # ascend :
  if k == GLUT_KEY_PAGE_UP:
    ascend = True
  # quit the game :
  if k == 'x':
    print "'x' was pressed"
    exit(0)



def specialUp(k,x,y):
  """
  """
  global up, down, rollLeft, rollRight, ascend, descend

  # stop pitching up :
  if k == GLUT_KEY_UP:
    up = False
  
  # stop pitching down :
  if k == GLUT_KEY_DOWN:
    down = False
  
  # stop rolling right :
  if k == GLUT_KEY_RIGHT:
    rollRight = False
  
  # stop rolling left :
  if k == GLUT_KEY_LEFT:
    rollLeft = False

  # stop descending :
  if k == GLUT_KEY_PAGE_DOWN:
    descend = False
  
  # stop ascending :
  if k == GLUT_KEY_PAGE_UP:
    ascend = False


def mouse(button,state,x,y):
  """
  """
  global beginx,beginy,rotate,camDist

  if button == GLUT_LEFT_BUTTON and state == GLUT_DOWN:
    #print "Mouseclick <x,y> : <%i,%i>" % (x,y)
    rotate = 1
    beginx = x
    beginy = y
  if button == GLUT_LEFT_BUTTON and state == GLUT_UP:
    rotate = 0
  # move camera out :
  if button == 4:
    camDist += 0.2
  # move camera in :
  if button == 3:
    if camDist > 0.2:
      camDist -= 0.2
    else:
      camDist = 0.0


def motion(x,y):
  """
  """
  global rotx,roty,beginx,beginy,rotate,camera

  if rotate:
    diffx  = (x - beginx) / 100.0
    diffy  = (y - beginy) / 100.0
    rotx   = rotx + diffy
    roty   = roty + diffx
    beginx = x
    beginy = y
    camera.update(array([-diffy, -diffx, 0]))
    glutPostRedisplay()
    #print "Mouse movement <x,y> : <%i,%i>" % (x,y)
  

# main method :
if __name__ == '__main__':
  
  # initial window position :
  sx = 600
  sy = 300

  # open a window
  glutInit(sys.argv)
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH)
  glutInitWindowPosition(sx, sy)
  glutInitWindowSize(w, h)
  glutCreateWindow("spaceSimulator")
  glutDisplayFunc(display)
  glutMouseFunc(mouse)
  glutMotionFunc(motion)
  glutReshapeFunc(reshape)
  glutIdleFunc(idle)
  glutKeyboardFunc(key)
  glutKeyboardUpFunc(keyUp)
  glutSpecialFunc(special)
  glutSpecialUpFunc(specialUp)

  obj = OBJ('SpaceShip.obj', revz=True)
  
  # initialize
  init()

  # hand off control to event loop
  glutMainLoop()


