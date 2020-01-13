### IMPORT STATEMENTS
import math
import time
import threading

###GLOBAL CONSTANTS
PI, E, LIGHT = math.pi, math.e, 299458792.0
DEG = PI/180
INF=math.inf
pi, e, deg, inf = PI, E, DEG, INF
JIM = complex(0,1)

MIN_HELM = 0.0000001 * DEG # curves with less helm than this
# can be considered straight lines

LAND_PREC = 0.1 # Recommended value for Land.delta_space

INCH = 0.0254
FOOT = 12 * INCH
MILE = 5280 * FOOT

MINUTE = 60.0
HOUR = 60 * MINUTE
DAY = 24 * HOUR

MPH = float(MILE) / HOUR

GRAV = 9.81 # acceleration due to gravity at Earth's surface in SI units

### BASIC MATH FUNCTIONS

sin,cos,tan = math.sin, math.cos, math.tan
arcsin, arccos, arctan = math.asin, math.acos, math.atan
ln, exp = math.log, math.exp

def avg(v):
  '''<List of Numbers>
  Returns their Average (mean)'''
  nu=float(sum(v))
  return nu/len(v)

def sgn(x):
  if x<0:
    return -1
  else:
    return 1


def prod(v):
  '''<List> of Numbers
  Multiplies together every element in <v>
  and returns the result'''
  p = 1
  for i in range(len(v)):
    p *= v[i]
  return p

def xmult(lhop, x1, x2):
  '''《lHospital limit》, 《scalar》, 《other scalar》
  lhop is winner of l'Hospital's
  rule
  Multiplies the two POSSIBLY INFINITE scalars 《x1》 and 《x2》 in accordance with whatever lHospital's rule
  says 0*inf would be in this situation
  '''
  if 0 in [x1, x2]:
    if inf in [x1, x2]:
      return lhop
    elif -inf in [x1, x2]:
      return -1*lhop
    else:
      return x1*x2
  else:
    return x1*x2

def xdiv(lhop, nu, den):
  '''《lHospital limit》, 《numerator》, 《denominator》
  Similar to xmult, except it does Division rather than multiplication
  Can divide by zero without a math error
  '''
  if abs(nu)==abs(den):
    if abs(nu)==inf:
      return sgn(nu)*sgn(den)*lhop
    elif abs(nu)==0:
      return lhop
  elif den==0 and nu!=0:
    return nu*inf
  else:
    return nu/den

def eul(theta):
  '''<Angle> in radians
  Returns phasor (complex number)
  with magnitude of one and phase angle <theta>'''
  return complex(cos(theta), sin(theta))


def conj(p):
  '''<Complex Number>
  Returns its Complex Conjugate'''
  return complex(p.real , -1 * p.imag)

def angle(x):
  '''<Complex number>
  Returns the phase angle, in radians'''
  k = xdiv(0, x.imag, x.real)
  return arctan(k)

### COMPLEX CONSTANTS

#In this file, directions of helm
#will sometimes be stored as unit phasors
#for example 0+j will mean
#due right or due east.
#Sometimes these directions will be absolute to the land,
#other times they will be relative to the direction of the car 

# ABSOLUTE DIRECTIONS

NORTH = eul(0)
EAST = eul(PI / 2)
SOUTH = eul(PI)
WEST = eul(-1 * PI / 2)

# RELATIVE DIRECTIONS
#These can be converted to absolute directions
#if you multiply them by an absolute direction phasor
#to "rotate" them.
#For example,   NORTH * eul(30*DEG)  
#means 30 degrees east of north.

O_CLOCK = eul(30 * DEG)
#People sometimes say "o'clock" as a unit of helm
#for example, "1 o'clock" means 30 degrees rightward
#O_CLOCK is that unit as a direction phasor
#To find the direction phasor for 2 o'clock,
#you  would go  
#O_CLOCK ** 2

ONE_DEGREE = eul( DEG ) #one degree rightwards as a direction phasor
ONE_RADIAN = eul(1) #one radian, as a direction phasor


### VECTOR FUNCTIONS

def funxvec(funx, v): #ver 
  '''<Python function> with one input, <List of Inputs>
  Returns list of outputs of <funx> given the inputs in <v>
  in corresponding order'''
  vy = []
  for i in range(len(v)):
    a=v[i]
    vy.append(funx(a))
  return vy

def funxvec_2(funx, v1, v2): #ver 
  n=min(len(v1), len(v2))
  vf=[]
  for i in range(n):
    a1, a2 = v1[i], v2[i]
    f=funx(a1,a2)
    vf.append(f)
  return vf

def funxvec_3(funx, v1, v2, v3):
  n = min(len(v1), len(v2))
  n = min(n, len(v3))
  vf = []
  for i in range(n):
    a1, a2, a3 = v1[i], v2[i], v3[i]
    f = funx(a1, a2, a3)
    vf.append( f )
  return vf

def kvec(k, vec):
  v=[]
  for i in range(len(vec)):
    v.append(k * vec[i])
  return v

def podvec(el, size):
  '''<Element Value>, <List Length>
  Returns list with <size> elements, all with value of <el>'''
  v=[]
  for i in range(size):
    v.append(el)
  return v

def riemvec(ti, tf, n):
  '''<Initial Value> on independent axis, <Final Value> on independent axis, <Desired number of subdivisions>
  Returns Python list of midpoints on x-axis for the subdivisions of a Riemann sum
  running from x=ti to x=tf'''
  bar=float(tf-ti)/n
  v=[]
  for i in range(n):
    n_x=0.5+i
    x = ti + (bar * n_x)
    v.append(x)
  return v

def copy(vec ):
  '''<Python list>
  Returns copy of this list'''
  v=[]
  for i in range(len(vec)):
    v.append(vec[i])
  return v
  
def closest(x, vec):
  '''<Number>, <List>
  Returns value in <vec> closest to <x>'''
  vm=[]
  for i in range(len(vec)):
    vm.append(abs(x-vec[i]))
  ndx = vm.index(min(vm))
  return vec[ndx]

def num_in_list(x, vec):
  '''<Value>, <List>
  Returns number of times <x> appears in <vec>'''
  s=0
  for a in vec:
    s+=int(a==x)   
  return s

def val_index_vec(x, vec):
  '''<Value>, <List>
  Returns list with indices where <x> appears in <vec>'''
  ret=[]
  for i in range(len(vec)):
    if vec[i]==x:
      ret.append(i)   
  return ret
  

### MATH FUNCTIONS FOR THIS PROGRAM

#WARNING:  curve_ray   doesn't work on perfectly straight lines.
#It only works for curved paths.
def curve_ray(theta, kurv_r):
  '''<Angle in radians> subtended by a curve, <Curvature Radius> of that curve
  Returns phasor (complex number) 
  whose Magnitude is the distance from startpoint to endpoint of the curve as the crow flies,
  and whose Angle is the angle of that ray WRT driver's direction at startpoint
  '''
  mag_ang = 0.5 * (PI - theta)
  mag = 2 * kurv_r * cos(mag_ang)
  return mag * eul(theta / 2)


def to_ray(d_path, helm):
  '''#@param   d_path   is distance travelled on a path 
  #@param   helm   is the helm of a car travelling that path
  #Returns phasor representing net displacement of that path,
  #be it straight or curved'''
  ang = d_path * helm
  if ang < MIN_HELM:
    return d_path * eul(0)
  else:
    return curve_ray( ang , 1.0/helm)

def unity(x):
  '''<Number>
  Returns number equal to <x>'''
  return 1*x

def always_zero(x):
  '''<Number>
  Returns 0 no matter the input'''
  return 0

def always_one(x):
  '''<Number>
  Returns 1 no matter what'''
  return 1

def kinetic(mass, speed):
  '''<Mass>, <Speed> of moving particle in SI units   
  Returns the Kinetic Energy of that particle, in watt-seconds
  '''
  return 0.5 * mass * (speed ** 2)

def velocity(speed, helm):
  '''<Speed>, <Helm> of moving object in SI units 
  Returns phasor representing the velocity,
  as a phasor representing the arc subtended in 1 second'''
  if helm < MIN_HELM:
    return speed * eul(0)
  else:
    return curve_ray(speed, speed * helm)



### FUNXION CLASS

class Funxion:

  def __init__(self, funx, initial_val):
    '''<Python function> with one arg, <Initial value> for Python function
    The function <funx> represents a mathematical function. Instantiates
    an object to hold that function, initialized to have input value <init_val>
    '''
    self.funx = funx
    self.in_val = 0.0 + initial_val
    self.out_val = self.funx(self.in_val)
    self.name = "funxion"
    self.BAR_NUM = 10000 #number of subdivisions for Riemann sum

  def sweep(self):
    self.out_val = self.funx(self.in_val)

  def feed(self, new_val):
    self.in_val = new_val
    self.sweep()
  
  def riemann(self, ti, tf, n):
    '''<Initial Value> on independedt axis, <Final Value> on independent axis, <Number of Subdivisions> for Riemann sum
    Returns Riemann sum, or approximate integral of funx(t) from t=ti to t=tf'''
    vx=riemvec(ti,tf,n)
    vy=funxvec(self.funx, vx)
    yav=avg(vy)
    return yav*(tf-ti)

  def riem(self, ti, tf):
    '''<Initial Value>, <Final Value>
    Returns the approximate integral of self.funx(t)
    from t = <ti> through t = <tf>'''
    return self.riemann(ti, tf, self.BAR_NUM) 


  def deriv(self, delta_x):
    '''<Differential on x-axis> object will use to find nearby values
    Approximates the derivative of self.funx's mathematical equivalent
    at t = self.in_val using nearby x-values'''
    #delta_y = self.funx(self.in_val+delta_x) - self.funx(self.in_val-delta_x)
    #print("\nself.in_val = "+str(self.in_val))
    #print("\nself.out_val = "+str(self.out_val))
    x_over, x_under = self.in_val + delta_x, self.in_val - delta_x
    #print("\nx_over = "+str(x_over)+"\nx_under = "+str(x_under))
    y_over, y_under = self.funx(x_over), self.funx(x_under)
    #print("\ny_over = "+str(y_over)+"\ny_under = "+str(y_under))
    ddx = float(y_over - y_under) / (x_over - x_under)
    #print("\nddx = "+str(ddx)+"\n\n")
    return ddx

  def refunx(self, funx_new):
    '''<Python function> to replace self.funx
    Changes self.funx and readjusts this Funxion object accordingly'''
    self.funx = funx_new
    self.sweep()

  def rename(self, new_name):
    '''<String> representing new name 
    Renames this Funxion object'''
    self.name = new_name
    self.sweep() 

  def copy(self):
    '''Creates new Funxion object with same instance variables as this one'''
    return Funxion(self.funx, self.in_val)

  def copy_name(self, new_name):
    '''<String> representing name of new Funxion object
    Copies this Funxion and gives that object its own name'''
    fu = self.copy() 
    fu.rename(new_name)

  def litter(self, n):
    '''<Number of desired copies> of this Funxion object
    Returns list of Funxion objects containing <n> copies of this Funxion object'''
    v=[]
    for i in range(n):
      v.append(self.copy())
    return v   

  def scale(self, k):
    '''<Scaling factor>
    Returns Funxion object just like this one,
    except the encapsulated function is scaled by <k>'''
    def dummy(x):
      return k * self.funx(x)
    fu = Funxion(dummy, self.in_val)
    return fu

  def scale_x(self, k):
    '''<Scaling factor>
    Returns Funxion object like this one 
    except shrunk in the x-axis by a factor of <k>'''
    def dummy(x):
      return self.funx(k * x )
    return Funxion(dummy, self.in_val)

  def shift(self, k):
    '''<Number>
    Returns Funxion object like this one,
    except with <k> units aded to the output'''
    def dummy(x):
      return k + self.funx(x )
    return Funxion(dummy, self.in_val)

  def shift_x(self, k):
    '''<NUmber>
    Returns FUnxion object like this one,
    but representing y(x + <k>) instead of y(x)'''
    def dummy(x):
      return self.funx(x + k )   
    return Funxion(dummy, self.in_val)





  def tostring(self):
    strega="Funxion   "+self.name+"  has these instance variables:\n"
    strega+="self.in_val = "+str(self.in_val)+"\nself.out_val = "+str(self.out_val)
    strega+="\nEND Funxion   "+self.name+"   instance variables\n\n"
    return strega

### FUNCTIONS THAT RETURN FUNXION OBJECTS 

def Konstant(k):
  '''<Constant>
  Returns a Funxion object
  whose output is always <k>
  '''
  def dummy(x):
    return k 
  fu = Funxion(dummy, 0)
  return fu

def Ramp(k):
  '''<Constant>
  Returns Funxion object
  representing the math function 
  y(x) = <k> * x'''
  def dummy(x):
    return k * unity( x )   
  fu = Funxion(dummy, 0)
  return fu

def SumFunxion(funxionvec, in_val):
  '''<List of Funxion objects>, <Initial Input Value>
  Returns a Funxion object whose output is the sum
  of several Funxion outpuots with the same input '''
  def dummy(x):
    yvec = []
    for funxion in funxionvec:
      #funxion.feed( x )
      #yvec.append(funxion.out_val)
      #EXPERIMENTAL LOOP BODY
      yvec.append(funxion.funx(x))
      #END EXPERIMENTAL LOOP BODY
    return sum(yvec)
  fu = Funxion(dummy, in_val)
  return fu

def Linear(slope, y_i):
  '''<slope>, <y-intercept>
  Returns Funxion object encapsulating linear function
  y = <y_i> + (<slope> * x)'''
  fv = [Ramp(slope), Konstant(y_i)]
  return SumFunxion(fv, 0)
  
def ProdFunxion(funxionvec, in_val):
  '''<List of Funxion objects>, <Initial Input Value>
  Returns a Funxion object whose output is the Product
  of several Funxion outputs with the same input '''
  def dummy(x):
    yvec = []
    for funxion in funxionvec:
      #funxion.feed( x )
      #yvec.append(funxion.out_val)
      #EXPERIMENTAL LOOP BODY
      yvec.append(funxion.funx(x))
      #END EXPERIMENTAL LOOP BODY
    return prod(yvec)
  fu = Funxion(dummy, 0)
  return fu


def Heav(x_rise):
  '''<Rise time> on independent axis   
  Returns Funxion object representing unit step function
  u(x - <x_rise>)'''
  def dummy(x):
    if x < x_rise:
      return 0
    else:
      return 1
  return Funxion(dummy, 0)

def Pulse(x_rise, x_fall):
  '''<Rise TIme>, <Fall Time> on independent axis   
  Returns Funxion object representing rectangular pulse function'''
  #return SumFunxion([Heav(x_rise), Heav(x_fall).scale(-1)], 0)
  #EXPERIMENTAL CODE BODY:
  def dummy(x):
    return int((x >= x_rise) and (x_fall >= x))
  return Funxion(dummy, 0)

def Exp():
  '''Returns Funxion object representing 
  equation y(x) = e^x'''
  return Funxion(exp, 0)

def Sinusoid(w, ph):
  '''<Angular Frequency>, <Phase Shift>
  Retuerns Funxion object representing sinusiodal function
  y(x) = sin(<w> x + <ph> )'''
  def dummy(x):
    return sin( (w *x) + ph)
  return Funxion(dummy, 0)

def Pow(k):
  '''<Power or exponent>
  Returns Funxion object representing equationy(x) = x ^ <k> '''
  def dummy(x):
    return x**k
  return Funxion(dummy, 0)

def Log():
  '''Returns Funxion object
  representing natural logarithm function.
  Uses absolute value of input to avoid crashing at negative numbers'''
  def dummy(x):
    return ln( abs(x) )
  return Funxion(dummy, 1) #initial input is 1 instead of 0

def Periodic(funxion, period):
  '''<Funxion object>, <Period>
  Returns Funxion object representing periodic function
  with period <period>, 
  defined as y(x) = <funxion>.funx(x) for 0 <= x <= <period>'''
  def dummy(x):
    return funxion.funx(x % period   )  
  return Funxion(dummy, 0)



### OBSTACLE CLASS

class Obstacle:
  
  def __init__(self, funxion, position):
    '''<Funxion object> representing spatial borders of obstacle, <phasor> representing initial position of obstacle
    Returns object representing
    the space occupied by an obstacle or vehicle'''
    self.name = "obstacle"
    self.funxion = funxion
    self.position = position
    self.orientation = NORTH 
    
  def ray_length(self, theta):
    '''<Angle> in radians
    Returns distance from centerpoint to border at that angle'''
    th = theta % (2 * PI)
    return self.funxion.funx(th)

  def move(self, delta_position):
    '''<Phasor> representing change in position
    Changes the absolute position of this object'''
    self.position += delta_position

  def turn(self, ang):
    '''<Angle> in radians rightward
    Rotates orientation of this Obstacle by <ang> radians'''
    self.orientation *= eul(ang)

  def within(self, p):
    '''<Position phasor>
    Retuns boolean indicating whether <p>
    is within borders of this Obstacle'''
    r = p - self.position
    r_adj = r / self.orientation
    ang = angle(r_adj)
    return (self.ray_length(ang) >= abs(r))

  def borderpt(self, ang):
    '''<Angle>
    Returns position phasor representing position phasor of point 
    on this Obstacle's border, <ang> radians rightward from self.orientatiin'''
    pt = 0
    pt += self.position
    pt += self.ray_length(ang) * eul(ang)
    return pt

  def bordervec(self):
    '''Returns list of points 
    on border of this Obstacle,
    as position phasors'''
    xv=riemvec(0, 2*PI, self.funxion.BAR_NUM)
    k = 0 + xv[0]
    for x in xv:
      x -= k
    return funxvec(self.funxion.funx, xv)


  def copy(self):
    '''Returns copy of this Obstacle objects'''
    obst = Obstacle(self.funxion, self.position)
    obst.name, obst.orientation = self.name, self.orientation
    return obst

  def rename(self, new_name):
    self.name = new_name

  def copy_name(self, new_name):
    obst = self.copy() 
    obst.rename(new_name)
    return obst

  def tostring(self):
    strega = "Obstacle   "+self.name+"   has these instance variables:\n"
    strega += "self.funxion = "+str(self.funxion)+"\n"
    strega += "self.position = "+str(self.position)+"\n"
    strega += "self.orientation = "+str(self.orientation)+"\n"
    strega += "END Obstacle   "+self.name+"   instance variables\n\n\n"
    return strega
    

### FUNCTIONS THAT RETURN OBSTACLE OBJECTS
def Circle(radius):
  '''<Radius> in meters open
  Returns Obstacle object representing a circle shape
  with radius <radius> meters'''
  return Obstacle(Konstant(radius), 0 * eul(0))

def Rectangle(shape_length, shape_width):
  '''<Length of rectangle> in direction flush with self.orientation, <Width of rectangle> in direction normal to self.orientation
  Returns Obstacle object representing a rectangle
  of length <shape_length> meters, width <shape_width> meters'''
  def dummy(theta):
    th_corn = arctan(shape_width / shape_length)
    th = theta % PI
    if th_corn < th < (PI - th_corn):
      return math.hypot( 0.5*shape_width, 0.5*shape_length*cos(th))
    else:
      return math.hypot(0.5*shape_length, 0.5*shape_width*sin(th))
  fu = Funxion(dummy, 0)
  obst = Obstacle(fu, 0*eul(0))
  return obst


      

### CAR CLASS

class Car:

  def __init__(self, mass, pwr_drive, energy):
    '''<Mass of car>, <Intentional power on car>, <Kinetic energy of car>'''
    self.name = "car" #default name for this Car object
    self.mass, self.pwr_drive, self.energy = mass, pwr_drive, energy
    self.pos, self.delta_pos = 0*eul(0), 0*eul(0) #Position, change in position, stored as a phasor
    self.stop = 0 #boolean true if car is intentionally stopped
    self.brake = 0 #boolean true if driver is intentionally braking
    self.rev = 0 #boolean true if car is in reverse
    self.tilt = 0+0.0 #downward tilt of the ground, in radians, 0 for level ground
    self.speed = (2.0 * self.energy /self.mass) ** (0.5) #speed, from kinetic energy
    self.pwr_grav = GRAV * sin(self.tilt) * self.mass * self.speed #mechanical power due to gravity
    self.pwr_other = 0 + 0.0 #mechanical power from sources other than gravity, the engine, or the brakes
    self.pwr_net = self.pwr_drive + self.pwr_grav + self.pwr_other
    # INSTANCE VARIABLES FOR DISTANCE CALCULATIONS
    self.delta_path = 0 + 0.0 #Distance travelled along path in a short time
    self.path = 0 + 0.0  #Total distance travelled since instantiation
    self.helm = 0 + 0.0 #Horizontal driving angle, in Radians per Meter rightward
    self.is_curved = abs(self.helm) > MIN_HELM #Boolean true if car's path is curved
    self.delta_tau = 0 #Short time delay for actual time delays like time.sleep()
    self.t_poll = 0.1 #Short time delay for predicting the car‘s motion
    self.t_now = 0.0 #Amount of time car has been in motion
    self.shape = Rectangle(5,3)
    self.frix_coeff = 0 #Coefficient of kinetic friction 
    
    
    


  
  def sweep(self):
    self.speed = (2.0 * self.energy / self.mass) ** (0.5)
    self.brake = self.pwr_drive < 0 #Set brake flag if driver intentionally reducing energy
    self.pwr_grav = GRAV * sin(self.tilt) * self.mass * self.speed
    self.pwr_net = self.pwr_grav + self.pwr_drive + self.pwr_other
    self.is_curved = abs(self.helm) > MIN_HELM
    if (self.speed <= 0) and (sgn(self.pwr_drive) < 0) :
      self.stop = 1
      self.energy = 0.0
    if (self.pwr_drive > 0 ):
      self.stop = 0
      self.brake = 0
    self.shape.position = self.pos 
    #CODE NOTE: The instance variables self.path and self.delta_path are NOT affected
    #by the sweep() method, because sweep method does not handle the motion itself
    #and also because recalculating the displacement too often can lead to inaccurate results

  def reverse(self):
    '''Puts this Car in reverse'''
    if not self.rev:
      self.rev = 1
      self.helm *= -1
      self.tilt *= -1
      self.pwr_grav *= -1
    self.sweep() 
  
  def forward(self):
    '''Takes this Car out of reverse'''
    if self.rev:
      self.rev = 0
      self.helm *= -1
      self.tilt *= -1
      self.pwr_grav *= -1
    self.sweep() 

  def friction(self, m_k):
    '''<Coefficient of Kinetic Friction>
    Calculates kinetic friction, in SI units, of this car
    using the simple model taught in beginner-level physics'''
    weight = GRAV * self.mass 
    weight *= cos( self.tilt )
    return -1 * abs(m_k * weight * self.speed  )

  def frict_encaps(self, kin_energy):
    '''<Kinetic Energy> in SI units   
    Car.friction method, but encapsulated 
    as a function for use in Funxion objects
    NOTE : <kin_eneergy> is a dummy parameter, 
    this function actually does all calculations using instance variables'''
    return self.friction(self.frix_coeff)

  def respeed(self, new_speed):
    '''<New speed> of Car in SI units   
    Changes the speed of the Car to <new_speed> 
    and adjusts instance variables to match.
    WARNING: Does not simulate any motion,
    just changes the speed'''
    self.energy = kinetic(self.mass, new_speed)
    self.sweep()

    
  def other_fxn(self):
    #Returns default Funxion object  
    #for use in Land objects
    #to calculate self.pwr_other
    #for this Car object   
    return Funxion(self.frict_encaps, self.energy)
    
  def drive_fxn(self):
    #Returns default Funxion object   
    #for use in Land objects
    #to output self.pwr_drive   
    #for this Car object   
    def dummy(t):
      return self.pwr_drive  
    return Funxion(dummy, self.t_now)



  

  def go(self, delta_t):
    '''<Short amount of time> in secondss
    Changes the instance variables to model a small amount of motion in the car'''
    energy_initial = 0.0 + self.energy
    d_t = 0.0 + delta_t
    self.sweep()
    if (self.pwr_net < 0) and ( abs(self.pwr_net * delta_t) > energy_initial):
      d_t = energy_initial / self.pwr_net #amount of time it will take car to stop at this rate
    time.sleep(self.delta_tau) #delay for threading purposes
    self.energy += self.pwr_net * d_t
    #
    self.t_now += d_t
    #
    self.sweep()
    energy_avg = avg([energy_initial, self.energy]) #average energy during the short time, 
    #assuming constant mechanical power throughout
    speed_avg = (2.0 * energy_avg / self.mass) ** 0.5 #average speed during the time elapsed 
    self.delta_path = speed_avg * d_t
    curv_theta = self.helm * self.delta_path
    curv_rad = xdiv(0, 1.0, self.helm)
    if not self.is_curved:
      self.delta_pos = self.delta_path * eul(0)
    else:
      self.delta_pos = curve_ray(curv_theta, curv_rad)
    if not self.rev:
      self.path += self.delta_path
      self.pos += self.delta_pos
      #EXPERIMENTAL CODE
      self.shape.move(self.delta_pos)
      self.shape.turn(2 * angle(self.delta_pos))
      #END EXPERIMENTAL CODE
    else:
      self.path -= self.delta_path
      self.pos -= self.delta_pos
      #EXPERIMENTAL CODE
      self.shape.move(self.delta_pos)
      self.shape.turn(2 * angle(self.delta_pos))
      #END EXPERIMENTAL CODE
    

  def go_drive(self, delta_t, drivefunxion):
    '''<Amount of time> in seconds for motion, <Funxion objet> representing self.pwr_drive as a  function of time
    Uses self.go( <delta_t> ) but decides self.pwr_driveby the function encapsulated in <drivefunxion>'''
    drivefunxion.feed(self.t_now)
    self.pwr_drive = drivefunxion.out_val
    self.sweep() 
    self.go(delta_t)
    
  def go_other(self, drivefunxion, otherfunxion): 
    '''#@param   drivefunxion   is Funxion object
    #representing self.pwr_drive as a function of Time
    #@param   otherfunxion   is Funxion object
    #representing self.pwr_other as a function of Energy'''   
    otherfunxion.feed(self.energy)
    self.pwr_other = otherfunxion.out_val
    self.sweep()
    self.go_drive(self.t_poll, drivefunxion)

    
  def travel_raw(self, dur, delta_t): #NOT RECOMMENDED
    '''<Duration> in seconds, <Polling Period> in seconds
    Iterates the method Car.go(delta_t) until <dur> seconds have supposedly elapsed'''
    n = int(dur // delta_t)
    for i in range(n):
      self.go(delta_t)

  def travel(self, dur, delta_t): #NOT RECOMMENDED
    '''<Duration> in seconds, <Polling Period> in seconds
    Iterates the method Car.go(float), but in a threaded manner,
    so other code can execute at the same time'''
    thr = threading.Thread(target = self.travel_raw, name = "", args = (dur, delta_t))
    thr.start()
    
  def to_speed_lin(self, speed_des, accel_time): #NOT RECOMMENDED
    '''#@param   speed_des   is desired speed
    #@param   accel_time   is time it should take to speed up/slow down
    #Gets this Car object to desired speed
    #by iterating Car.go( float ) at constant self.pwrdrive
    '''
    energy_des = 0.5 * self.mass * (speed_des ** 2)
    energy_delta = energy_des - self.energy
    pwr_net_des = energy_delta / accel_time
    self.pwr_drive = pwr_net_des - self.pwr_grav - self.pwr_other
    self.sweep()
    self.travel(accel_time, self.t_poll)
    

  def str_motion(self):
    '''Returns string with attributes directly relevant to the motion'''
    strega = ""
    strega+="self.speed = "+str(self.speed)+"  meters per second\n"
    strega+="self.pwr_drive = "+str(self.pwr_drive)+"  watts\n"
    strega+="self.pwr_net = "+str(self.pwr_net)+"  watts\n"
    strega+="self.energy = "+str(self.energy)+"  watt-seconds\n"
    strega+="self.path = "+str(self.path)+"  meters \n"
    strega+="self.delta_path = "+str(self.delta_path)+"  meters\n"
    strega+="self.pos = "+str(self.pos)+"  meters\n"
    strega+="self.delta_pos = "+str(self.delta_pos)+"  meters\n"
    strega += "self.helm = "+str(self.helm)+"  radians per meter\n"
    strega += "self.t_now = "+str(self.t_now)+"  seconds\n"
    strega+="\n\n\n"
    return strega

  def str_power(self):
    '''Returns string containig all the 
    power-related attributes of this Car object'''
    strega = "self.pwr_net = "+str(self.pwr_net)+"  watts\n"
    strega += "self.pwr_drive = "+str(self.pwr_drive)+"  watts\n"
    strega += "self.pwr_other = "+str(self.pwr_other)+"  watts\n"
    strega += "self.pwr_grav = "+str(self.pwr_grav) +"  watts\n"
    strega += "self.tilt = "+str(self.tilt)+"  radians downward\n"
    strega += "self.energy = "+str(self.energy)+" watt-seconds\n\n\n"
    return strega 

  def str_time(self):
    '''Returns string containing all time-related attributes
    of this Car object'''
    strega = "self.delta_tau = "+str(self.delta_tau)+"  seconds of processor time\n"
    strega += "self.t_poll = "+str(self.t_poll)+"  seconds of simulated time\n"
    strega += "self.t_now = "+str(self.t_now)+"  seconds of simulated time\n\n\n"
    return strega

  def str_flags(self):
    '''Returns string cpontaining all boolean attributes
    of this Car object'''
    strega = "self.brake = "+str(int(self.brake))+"  boolean\n"
    strega += "self.stop = "+str(int(self.stop)) +"  boolean\n"
    strega += "self.rev = "+str(int(self.rev)) +"  boolean\n"
    strega += "self.is_curved = "+str(int(self.is_curved))+"  boolean\n\n\n"
    return strega


  def tell_motion_raw(self, dur, d_tau):
    '''<Duration> in seconds you want to run this method, <Polling Period> in seconds between print statements
    Prints attrributes related to the motion of the car every <d_tau> seconds, 
    until <dur> seconds have elapsed'''
    n = int(dur // d_tau)
    for i in range(n):
      print(self.str_motion())
      time.sleep(d_tau)

  def tell_motion(self, dur, d_tau):
    '''<Duration> in seconds this method will run, <Polling Period>
    Prints attributes related to this Car's motion continually every <d_tau> seconds,
    but in a threaded manner'''
    thr = threading.Thread(target = self.tell_motion_raw, name="", args = (dur, d_tau))
    thr.start()





    
  def rename(self, new_name):
    self.name = new_name
    self.sweep() 
  
  def copy(self):
    car = Car(self.mass, self.pwr_drive, self.energy)
    car.stop, car.brake, car.rev = self.stop, self.brake, self.rev
    car.tilt, car.helm = self.tilt, self.helm
    car.pwr_grav, car.pwr_other, car.pwr_net = self.pwr_grav, self.pwr_other, self.pwr_net
    car.speed, car.path, car.delta_path = self.speed, self.path, self.delta_path
    car.pos, car.delta_pos = self.pos, self.delta_pos
    car.is_curved, car.delta_tau = self.is_curved, self.delta_tau
    car.t_poll, car.t_now = self.t_poll, self.t_now
    car.frix_coeff, car.shape = car.frix_coeff, car.shape.copy()
    return car

  def copy_name(self, new_name):
    '''<String> representing name of new object
    Copies this Car object and gives resulting object a new name'''
    car = self.copy() 
    car.rename(new_name)
    return car

  
    
  
    

    
  
  def tostring(self):
    strega = "Car   "+self.name+"   has these attributes:\n"
    strega += "self.energy = "+str(self.energy)+"  watt-seconds\n"
    strega += "self.pwr_drive = "+str(self.pwr_drive)+"  watts\n"
    strega += "self.pwr_grav = "+str(self.pwr_grav)+"  watts\n"
    strega += "self.pwr_other = "+str(self.pwr_other)+"  watts\n"
    strega += "self.pwr_net = "+str(self.pwr_net)+"  watts\n"
    strega += "self.tilt = "+str(self.tilt)+"  radians downhill\n"
    strega += "self.mass = "+str(self.mass)+"  kilograms\n"
    strega += "self.frix_coeff = "+str(self.frix_coeff)+"  unitless\n"
    strega += "self.stop = "+str(self.stop)+"  boolean\n"
    strega += "self.brake = "+str(self.brake)+"  boolean\n"
    strega += "self.rev = "+str(self.rev)+"  boolean\n"
    strega += "self.is_curved = "+str(self.is_curved)+"  boolean\n"
    strega += "self.helm = "+str(self.helm)+"  radians per meter rightward\n"
    strega += "self.path = "+str(self.path)+"  meters\n"
    strega += "self.pos = "+str(self.pos)+" meters\n"
    strega += "self.speed = "+str(self.speed)+"  meters per second\n"
    strega += "self.delta_path = "+str(self.delta_path)+"  meters\n"
    strega += "self.delta_pos = "+str(self.delta_pos)+"  meters\n"
    strega += "self.delta_tau = "+str(self.delta_tau)+"  seconds\n"
    strega += "self.t_poll = "+str(self.t_poll)+  "  seconds\n"
    strega += "self.t_now = "+str(self.t_now)+"  seconds\n"
    strega += "self.shape = "+str(self.shape)+"  Obstacle object\n"
    strega += "END Car   "+self.name+"   attributes\n\n\n"
    return strega

### CAR - RETURNING FUNCTIONS

def Blankcar():
  '''Returns Car with default instance variables'''
  return Car(1000,0,0)

def named_cars(name_vec):
  '''<List of Strings>
  Returns list of Car objects 
  with those names'''
  x = []
  for i in range(len(name_vec)):
    x.append(Blankcar())
    x[i].rename(name_vec[i])
  return x

### LAND OBJECT 

class Land:
  
  def __init__(self, tiltvec, helmvec, delta_space):
    self.name = "land"
    self.tiltvec, self.helmvec = tiltvec, helmvec
    self.tilt_avg = avg(self.tiltvec) #average tilt
    self.helm_avg = avg(self.helmvec) #average helm
    #self.tiltvec holds the hill angles, self.helmvec holds the curvature of the road in the same units as Car.helm and Car.tilt
    self.delta_space =  delta_space
    self.size = min(len(self.tiltvec), len(self.helmvec)) # vec
    self.frixvec = podvec(0, self.size) #holds coefficients of friction
    #in case I ever want to model that
    #EXPERIMENTAL CODE    
    self.frix_avg = avg(self.frixvec)  
    #END EXPERIMENTAL CODE     
    self.compass = NORTH #orientation of the land AKA driver’s direction upon entering, as a unit phasor, assumed NORTH by default
    self.carvec = [] #list of Car objects inteacting with this Land object
    #assumed empty by default
    self.phasorvec = []
    self.carspacevec = [] 
    self.othervec, self.drivevec = [],[]
    #EXPERIMENTAL CODE
    self.land_position = 0*eul(0)  #Absolute position of this Land's startpoint,
    #in same units as Obstacle.position
    #END EXPERIMENTAL CODE
    for i in range(self.size):
      ang = self.helmvec[i]*self.delta_space
      if abs(ang) < MIN_HELM:
        self.phasorvec.append(self.delta_space * eul(0))
      else:
        curv_r = 1.0 / abs(self.helmvec[i])
        self.phasorvec.append(curve_ray(ang, curv_r))
    self.ray = to_ray(self.size*self.delta_space, avg(self.helmvec)) 
    #Phasor representing net displacement of path, in same units as Car.pos
    self.tau_max = 0 #Processor time constant
    self.t_land = 1/50 #Polling rate, in simulated seconds, for Car objects in this Land 
    #EXPERIMENTAL INSTANCE VARIABLES   
    self.outlandvec=[]
    self.outangvec=[]
    self.inlandvec=[]
    self.inangvec=[]
    self.inspacevec=[]
    self.has_in = 0 #boolean, only true if other Land objects frow into this one
    self.has_out = 0 #boolean, only true if this Land object flows into others

  def rename(self, name):
    self.name = name
    
    
  def ndx_point(self, path):
    '''<Distance travelled> along path represented by this Land object, in same units as Car.path
    Returns index of self.tiltvec and self.helmvec corresponding to that point along the path
    '''
    n = int(path // self.delta_space)
    n = min(n, self.size - 1)
    return n 

  def sweep(self):
    for car in self.carvec:
      car.delta_tau = self.tau_max
      car.t_poll = self.t_land
      ndx_pt = self.ndx_point( self.carspacevec[ self.carvec.index(car) ])
      car.helm = self.helmvec[ndx_pt]
      car.tilt = self.tiltvec[ndx_pt]
      car.frix_coeff = self.frixvec[ndx_pt]
      car.sweep()
    while(len(self.tiltvec) < len(self.helmvec)):
      self.tiltvec.append(self.tiltvec[-1])
    while(len(self.helmvec) < len(self.tiltvec)):
      self.helmvec.append(self.helmvec[-1])
    self.size = min(len(self.helmvec) , len(self.tiltvec))
    while( len(self.frixvec) < self.size):
      self.frixvec.append(self.frixvec[-1])
    for i in range(self.size):
      ang = self.delta_space * self.helmvec[i]
      if ang < MIN_HELM:
        self.phasorvec[i] = self.delta_space * eul(ang) 
      else:
        curv_r = 1.0 / self.helmvec[i]
        self.phasorvec[i] = curve_ray(ang, curv_r)
    while( len(self.phasorvec) > self.size):
      self.phasorvec.pop(-1)
    self.ray = to_ray(self.size*self.delta_space, avg(self.helmvec))
    self.tilt_avg = avg(self.tiltvec)
    self.helm_avg = avg(self.helmvec)
    self.frix_avg = avg(self.frixvec)
    #EXPERIMENTAL CODE
    self.has_in = len(self.inlandvec) > 0
    self.has_out = len(self.outlandvec) > 0

  def point_ray(self, space):
    '''<Distance travelled> along path represnted by theis Land object
    Returns phasor representing net displacement between beginning of path
    and that point in space''' 
    ndx = self.ndx_point(space)
    d_upto = ndx * self.delta_space
    sh = 0
    for i in range(ndx):
      sh += self.helmvec[i]
    h_avg = xdiv(0, float(sh), ndx)
    r = to_ray(d_upto, h_avg)
    r += to_ray(space - d_upto, self.helmvec[ndx])
    return r

    
  def drxn_ray(self, space):
    #@param   space   is a float
    #   representing how far along the path of this Land object you are
    #Returns direction phasor representing
    #absolute direction your car is pointing,
    #in radians rightward of due NORTH
    drxn = eul(0) * self.compass
    ndx = self.ndx_point(space)
    for i in range(ndx):
      drxn *= eul(self.delta_space * self.helmvec[i])
    ds_last = space - (ndx * self.delta_space)
    drxn *= eul(ds_last * self.helmvec[ndx])
    return drxn

  def recompass(self, new_compass):
    '''<Direction phasor> representing new value for self.new_compass
    Reorients this Land object and the cars in it'''
    redrxn = new_compass / self.compass 
    self.compass = new_compass    
    for car in self.carvec:
      car.shape.orientation *= redrxn
    self.sweep() 

  def endspace(self):
    '''Returns toal amount of distance travelled
    on road represented by this Land, in Car.path units'''
    return self.delta_space * self.size

  def reposition(self, entrance_pos):
    '''<Absolute position> of this Land's startpoint, in same units as Obstacle.position
    CHanges the startpoint of this Land object,
    and adjusts Car.pos and Car.obst.position for every Car in here'''
    self.land_position = entrance_pos
    for car in self.carvec:
      ndx_car = self.carvec.index(car)
      z = self.point_ray(self.carspacevec[ndx_car])
      car.pos = z + self.land_position
      car.sweep()
    self.sweep()

    
  def extend(self, tilt_v, helm_v ):
    #Adds more points to this Land object
    #@param   tilt_v   is list of values
    #   to be added to self.tiltvec
    #@param   helm_v   is list of values 
    #   to be added to self.helmvec
    n = min(len(tilt_v), len(helm_v))
    for i in range(n):
      self.tiltvec.append(tilt_v[i])
      self.helmvec.append(helm_v[i])
      self.frixvec.append(self.frixvec[-1])
      self.phasorvec.append(to_ray(self.delta_space, helm_v[i]))
    self.sweep()
    
  def concat(self, land):
    print("Land.concat active\n\n")
    size_og = 0+self.size
    self.extend(land.tiltvec, land.helmvec)
    self.sweep()
    for i in range(size_og, self.size):
      ix = i - size_og
      fr = land.frixvec[ix]
      self.frixvec.append(fr)
    self.sweep()
    print("Land.concat finished\n\n\n")
    
  def conform(self, prec):
    #@param   prec   is new value for self.delta_space
    #Makes it so this Land object has self.delta_space = <prec>
    #but it represents the same Physical characteristics of this Land object.
    totalspace = self.delta_space * self.size   
    pasado = 0 + 0.0   
    tv, hv, frv = [],[],[]
    while( pasado < totalspace):
      ndx = self.ndx_point(pasado)
      tv.append(self.tiltvec[ndx])
      hv.append(self.helmvec[ndx])
      frv.append(self.frixvec[ndx])
      pasado += prec   
    self.delta_space = prec   
    self.tiltvec = tv   
    self.helmvec = hv   
    self.frixvec = frv  
    #EXPERIMENTAL CODE to stop it from crashing
    while( len(self.phasorvec) < len(self.helmvec)):
      self.phasorvec.append(0)
    #END EXPERIMENTAL CODE
    self.sweep()
    
  def constant_tilt(self, tilt):
    #Sets this Land objects tilt vakues
    #to <tilt> at every point   
    for i in range(self.size):
      self.tiltvec[i] = tilt
    self.sweep()
    
  def constant_helm(self, helm):
    #Sets this Land object’s selm values
    #to <helm> at every point   
    for i in range(self.size):
      self.helmvec[i] = helm
    self.sweep()
    
  def constant_frix(self, frix):
    #Sets this Land object’s friction coefficients
    #to <frix> at every piint   
    for i in range(self.size):
      self.frixvec[i] = frix
    self.sweep()
    
      
    

    
  def enter_car(self, car, space):
    '''<Car object>, <Place on path> in Car.path units
    Puts a Car object to interact with this Land object'''
    self.carvec.append(car)
    self.carspacevec.append(space )
    #EXPERIMENTAL CODE   
    car.pos = self.land_position + (self.point_ray(space)*self.compass)
    self.sweep()

    
  def exit_car(self, car):
    '''<Car object>
    Removes that Car from this Land object if it's in here'''
    if car in self.carvec:
      ndx = self.carvec.index( car  )
      self.carvec.pop(ndx)
      self.carspacevec.pop(ndx)
      #self.drivevec.pop(ndx)
      #self.othervec.pop(ndx)
    self.sweep()

  
      
    
  def ir_raw(self, car,  space):
    '''<Car object>,  <Point where it is>  on the road in Car.path units
    Induces Car.go method in <car>,
    simulating self.t_land seconds of motion,
    starting <space> meters into the path, 
    then updates this Land object's instance variables'''
    if car not in self.carvec:
      self.enter_car(car, space)
    ndx_car = self.carvec.index(car)
    ndx_space = self.ndx_point(space)
    land_tilt, land_helm = self.tiltvec[ndx_space], self.helmvec[ndx_space]
    car.tilt = land_tilt
    car.helm = land_helm
    car.sweep()
    car.go(self.t_land)
    if not car.rev:
      self.carspacevec[ndx_car] += car.delta_path
    else:
      self.carspacevec[ndx_car] -= car.delta_path

  def ir(self, car,  space):
    '''<Car object>,  <Place> on the road, in car.path units 
    Threaded version of   Land.ir_raw( Car,  float )'''
    thr = threading.Thread(target = self.ir_raw, name = "", args = (car,  space))
    thr.start()

  def ir_todos(self):
    '''Does the method self.ir( Car, float )
    for every Car in this Land object,
    for the same amount of time'''
    #No params needed 
    #because we take those from this Land object's instance variables.
    #delta_t_vec = podvec(self.t_land, len(self.carvec))
    v = funxvec_2(self.ir, self.carvec,  self.carspacevec)
    time.sleep(self.tau_max)

  def ir_other_raw(self, drivefunxion, otherfunxion, ndx_car):
    #@param   drivefunxion   is Funxion object 
    #   representing Car.pwr_drive vs. time
    #@param   otherfunxion   is Funxion object 
    #   representing Car.pwr_other vs. kinetic energy
    #@param   ndx_car   is an int representing 
    #   the index of self.carvec
    #Induces the method Car.go_other( Funxion, Funxion )
    #in a Car object which is already in this Land object
    car = self.carvec[ndx_car]
    self.sweep() 
    car.go_other(drivefunxion, otherfunxion)
    signum = 1
    if car.rev:
      signum *= -1
    #car.path += signum * car.delta_path
    self.carspacevec[ndx_car] += signum * car.delta_path 

  def ir_other(self, drivefunxion, otherfunxion, ndx_car):
    thr = threading.Thread(target = self.ir_other_raw, name="", args = (drivefunxion, otherfunxion, ndx_car))
    thr.start()
    
  def ir_otros(self, drivevec, othervec):
    #@param   drivevec   is a List of Funxion objects 
    #   to be used to calculate Car.pwr_drive  
    #@param   othervec   ... calculate Car.pwr_other
    #Lists correspond to self.carvec
    #Does method self.ir_other,
    #but on every Car in this Land object  
    for i in range(len(self.carvec)):
      drive, other = drivevec[i], othervec[i]
      self.ir_other(drive, other, i)
    time.sleep(self.tau_max) 
    
  def ir_mios(self):
    #Does self.ir_otros
    #but using only the instance variables
    #in this Land object
    self.ir_otros(self.drivevec, self.othervec)
    
  def viaje(self, dur_t):
    '''#@param   dur_t   is number of seconds
    #of time you wish to simulate.  
    #Iterates the method self.ir_mios
    #for <dur_t> simulated seconds''' 
    n = int(dur_t // self.t_land)
    for i in range(n):
      self.ir_mios()
    self.sweep()

  def repitar(self, dur_t):
    '''<Number> of seconds of simulated time
    Repeats the method self.ir_todos for
    <dur_t> simulated seconds'''
    n = 1 + int(dur_t // self.t_land)
    for i in range(n):
      self.ir_todos()   
    self.sweep()
    
  def inflow(self, land, space):
    '''#@param   land   is Land object representing
    #   one of the roads flowing into the road
    #   represented by this Land object.
    #Makes this Land object recognize 
    #<land> as one of the roads joining into this one'''
    self.has_in, land.has_out = 1,1
    drxn_in=land.drxn_ray(land.endspace())
    drxn_me = self.drxn_ray(space)
    ang = angle(drxn_me / drxn_in) 
    land.outangvec.append(ang)
    self.inangvec.append(-1*ang)
    land.outlandvec.append(self)
    self.inlandvec.append(land)
    self.inspacevec.append(space)
    land.sweep()
    self.sweep()
    
  def outflow(self, land):
    '''#@param   land   is another Land object 
    #   that you want to make flow out of this Land.  
    #Makes this Land object recognize <land> as representing
    #one of the forks at the end of the road'''
    land.inflow(self, 0)

  def has_fanout(self):
    '''Boolean, returns True
    only if the outflows have outflows'''
    v = []
    for x in self.outlandvec:
      v.append(len(self.outlandvec))
    return (sum(v) > 0)

  def has_fanin(self):
    '''Boolean, returns True
    only if the inflows have inflows'''
    v = []
    for x in self.inlandvec:
      v.append(len(x.inlandvec))
    return (sum(v) > 0)

  def connectvec_1(self):
    v = [self]
    for inland in self.inlandvec:
      if inland not in v:
        v.append(inland)
    for outland in self.outlandvec:
      if outland not in v:
        v.append(outland)
    return v

  def connectvec_2(self):
    v = self.connectvec_1()
    for outland in self.outlandvec:
      v_out = self.connectvec_1()   
      for x_out in v_out:
        if x_out not in v:
          v.append(x_out)
    for inland in self.inlandvec:
      v_in = inland.connectvec_1()
      for x_in in v_in:
        if x_in not in v:
          v.append(x_in)
    return v     

  def connectvec_3(self):
    v = self.connectvec_2()     
    v_land = copy(v)     
    for land in v:
      for outland in land.outlandvec:
        if outland not in v_land:
          v_land.append(land)
      for inland in land.inlandvec:
        if inland not in v_land:
          v_land.append(inland)
    connectmat = []
    for x in v_land:
      connectmat.append(x.connectvec_2())
    for row in connectmat:
      for i in range(len(row)):
        if row[i] not in v_land:
          v_land.append(row[i])
    return v_land
    
  def connectvec(self):
    '''Returns list of all the Land objects connected to this one, connected to their connections, etc.'''
    v_con = [self] 
    v_sizes = [len(v_con)]
    for x_out in self.outlandvec:
      if x_out not in v_con:
        v_con.append(x_out)
    for x_in in self.inlandvec:
      if x_in not in v_con:
        v_con.append(x_in)
    v_sizes.append(len(v_con))
    while(v_sizes[-1] != v_sizes[-2]):
      for con in v_con:
        for xo in con.outlandvec:
          if xo not in v_con:
            v_con.append(xo)
        for xi in con.inlandvec:
          if xi not in v_con:
            v_con.append(xi)
        v_sizes.append(len(v_con))
    return v_con   

  def connectvec_out(self):
    v_out = [self]
    v_sizes = [len(v_out)]
    for x1 in v_out:
      for x2 in x1.outlandvec:
        if x2 not in v_out:
          v_out.append(x2)
    v_sizes.append(len(v_out))
    while( v_sizes[-1] != v_sizes[-2]):
      for x in v_out:
        for y in x.outlandvec:
          if y not in v_out:
            v_out.append(y)   
      v_sizes.append(len(v_out))
    return v_out

  def connectvec_in(self):
    v_in = [self] 
    v_sizes = [len(v_in)]
    for x1 in v_in:
      for x2 in x1.inlandvec:
        if x2 not in v_in:
          v_in.append(x2)
    v_sizes.append(len(v_in))
    while(v_sizes[-1] != v_sizes[-2]):
      for x in v_in:
        for y in x.inlandvec:
          if y not in v_in:
            v_in.apend(y)
      v_sizes.append(len(v_in))
    return v_in

  def repos_in(self):
    for i in range(len(self.inlandvec)):
      inland = self.inlandvec[i]
      z = self.land_position + self.point_ray(self.inspacevec[i])
      z -= inland.ray * inland.compass
      if z != inland.land_position:
        inland.reposition(z)
    self.sweep()   

  def repos_out(self):
    for i in range(len(self.outlandvec)):
      outland = self.outlandvec[i]
      ndx = outland.inlandvec.index(self)
      z = self.land_position + self.ray
      z -= self.point_ray(outland.inspacevec[ndx])*outland.compass
      if z != outland.land_position:
        outland.reposition(z)
    self.sweep()

  def repos_downstream(self):
    '''Like Land.repos_out except it does All downsrteam connections,
    not just immediate connections'''
    v = self.connectvec_out()
    for i in range(len(v)):
      v[i].repos_out()
    self.sweep() 

  def repos_upstream(self):
    '''Like Land.repos_downstream except it deals with inflows'''
    v = self.connectvec_in() 
    for i in range(len(v)):
      v[i].repos_in()  
    self.sweep()   

  def repos_all(self):
    self.repos_upstream()    
    self.repos_downstream()

  def inlands_before(self, space):
    '''<NUmber> in Car.path units
    Returns list of Land objects which
    Inflow to this one at some point Upstream from <space>'''
    v=[]
    for i in range(len(self.inlandvec)):
      if self.inspacevec[i] <= space:
        v.append(self.inlandvec[i])
    return v

  def inlands_after(self, space):
    '''<Number> in Car.path units   
    Returns list of Land objects which 
    Inflow to this one Downstream from <space>'''
    v = []
    for i in range(len(self.inlandvec)):
      if self.inspacevec[i] >= space:
        v.append(self.inlandvec[i])
    return v

  def inlands_near(self, space, distance):
    '''<Number> in Car.path units, <Number> in Car.path units
    Returns list of inflows within <distance> from <space>'''
    v = []
    for i in range(len(self.inlandvec)):
      d_land = abs(space - self.inspacevec[i])
      if d_land <= distance:
        v.append(self.inlandvec[i])
    return v

  def lands_after(self, space):
    '''<Number> in Car.path units
    Returns list of Lands crossing this one
    Downstream from <space>, whether inflow or outflow'''
    v = self.inlands_after(space)
    for x in self.outlandvec:
      v.append(x)
    return v

  def lands_near(self, space, distance):
    '''<Numbr> in Car.path units, <Number> in Car.path units
    Like Land.inlands_near, except it also includes nearby outflows'''
    v = self.inlands_near(space, distance)
    d_out = abs(self.endspace() - space)
    if d_out <= distance:
      for x in self.outlandvec:
        v.append(x)    
    return v

  def sorted_cars(self):
    '''<Returns list of Car onjects in this Land, 
    sorted from most upstream to most downstream'''
    space_v = copy(self.carspacevec)
    space_v.sort()
    car_v = []
    bookmark = 0
    for i in range(len(space_v)):
      if i >= bookmark:
        s = space_v[i]
        ndx_v = val_index_vec(s, self.carspacevec)
        for ndx in ndx_v:
          car_v.append(self.carvec[ndx])
        bookmark += len(ndx_v)
    return car_v



    


    
    
    


      




    
    
    
    


  
  def tostring(self):
    strega = "Land   "+self.name+"   has these instance variables: \n"
    strega += "self.delta_space = "+str(self.delta_space)+"\n"
    strega += "self.tilt_avg = "+str(self.tilt_avg)+"\n"
    strega += "self.tiltvec = "+str(self.tiltvec)+"\n"
    strega += "self.helm_avg = "+str(self.helm_avg)+"\n"
    strega += "self.helmvec = "+str(self.helmvec)+"\n"
    strega += "self.frix_avg = "+str(self.frixvec)+"\n"
    strega += "self.frixvec = "+str(self.frixvec)+"\n"
    strega += "self.ray = "+str(self.ray)+"\n"
    strega += "self.phasorvec = "+str(self.phasorvec)+"\n"
    strega += "self.compass = "+str(self.compass)+"\n"
    strega += "self.carvec = "+str(self.carvec)+"\n"
    strega += "self.tau_max = "+str(self.tau_max)+"\n"
    strega += "self.t_land = "+str(self.t_land)+"\n"
    strega += "self.inlandvec = "+str(self.inlandvec)+"\n"
    strega += "self.inspacevec = "+str(self.inspacevec)+"\n"
    strega += "self.inangvec = "+str(self.inangvec)+"\n"
    strega += "self.outlandvec = "+str(self.outlandvec)+"\n"
    strega += "self.outangvec = "+str(self.outangvec)+"\n"
    strega += "self.has_out = "+str(self.has_out)+"\n"
    strega += "self.has_in = "+str(self.has_in)+"\n"
    strega += "END Land   "+self.name+"   instance variables\n\n\n"
    return strega
    
### FUNCTIONS THAT RETURN LAND OBJECTS  

def Flatland(distance, prec):
  #@param   distance   is path length in Car.path units
  #@param   prec   is Land.delta_space for this Land object
  #Returns Land object representing perfectly flat, straight,
  #and frictionless road
  n = 1 + int(distance // prec)
  v, u = podvec(0, n), podvec(0,n)
  kans = Land(v, u, prec)
  kans.sweep()
  return kans
        

def Funxionland(tiltfunxion, helmfunxion, mileage):
  '''#@param   tiltfunxion   is Funxion object representing
  #   tilt vs. path   
  #@param   helmfunxion   is Funxion object representing
  #   helm vs. path
  #@param   mileage   is lengthof the path,
  #   in meters, represented by this Land object   
  #Returns Land object representing oath <mikeage> meters ling,
  #where the tilt and helm are determined by 
  #<tiltfunxion> and <helmfunxion>, respectively'''
  n = 1 + int(mileage // LAND_PREC)
  xv = []
  for i in range(n):
    xv.append( i * LAND_PREC)
  tv = funxvec(tiltfunxion.funx, xv)
  hv= funxvec(helmfunxion.funx, xv)
  return Land(tv,  hv, LAND_PREC)

def Uniformland(tilt, helm, mileage):
  '''<<Tilt>, <Helm>, <Distance covered> in Car.path units   
  Returns Land object with the same tilt and helm throughout, 
  for <mileage> meters'''
  return Funxionland(Konstant(tilt), Konstant(helm), mileage)

def Blankland():
  '''Returns Land object representing 
  one mile of straight, flat, frictionless road'''
  return Uniformland(0,0,MILE)

def Straighthill(tiltvec):
  '''<List of tilt values>
  Returns Land object representing a perfectly straight road
  on sloped land'''
  helmvec = podvec(0, len(tiltvec))
  return Land(tiltvec, helmvec, LAND_PREC)

def Flatcurve(helmvec):
  '''<List of helm values>
  Returns Land object representing a curved road
  on perfectly flat ground'''
  tiltvec = podvec(0, len(helmvec))
  return Land(tiltvec, helmvec, LAND_PREC)

### TERRAIN CLASS     
class Terrain:

  def __init__(self, land_center):
    '''<Land object> 
    <land_center> is Land object which will be used as a reference for all Land objects 
    and their connections in this Terrain object'''
    self.name = "terrain"
    self.land_center = land_center
    self.direction = self.land_center.compass
    self.place_center = self.land_center.land_position
    self.landvec = self.land_center.connectvec()
    self.tau_center = self.land_center.tau_max #processor time constant
    self.t_terrain = self.land_center.t_land #simulation time constant
    self.t_life = 0 + 0.0 #Equivalent to Car.t_now
    for land in self.landvec:
      land.tau_max = self.tau_center   
      land.t_land = self.t_terrain
      land.conform(self.land_center.delta_space)
      land.sweep()
    self.land_center.repos_all()
    self.landplacevec = []
    for i in range(len(self.landvec)):
      self.landplacevec.append(self.landvec[i].land_position)

    
  def sweep(self):
    self.direction = self.land_center.compass
    self.place_center = self.land_center.land_position
    self.landvec = self.land_center.connectvec()   
    self.land_center.repos_all()   
    while(len(self.landvec) > len(self.landplacevec)):
      self.landplacevec.append(0) #to avoid crashing
    for i in range(len(self.landvec)):
      land = self.landvec[i]
      self.landplacevec[i] = land.land_position
      if land.delta_space != self.land_center.delta_space:
        land.conform(self.land_center.delta_space)
      land.t_land = self.t_terrain
      land.tau_max = self.tau_center
      land.sweep

    
  def rotate(self, ang):
    '''<Angle>
    Rotates every Land object in this Terrain by
    <ang> degrees rightward'''
    redrxn = eul(ang)
    for land in self.landvec:
      new_compass = land.compass * redrxn
      land.recompass(new_compass)
    self.sweep()   

  def landsearch(self, car):
    ret = self.land_center
    for land in self.landvec:
      if car in land.carvec:
        ret = land
    return ret

  def has_car(self, car):
    '''<Car object>
    Returns boolean telling whether <car> is in this Terrain object'''
    la = self.landsearch(car)   
    return (car in la.carvec)

  def changeland(self, car, land):
    land_init = self.landsearch(car)
    space_over = land_init.carspacevec[land_init.carvec.index(car)] - land_init.endspace()
    if car.rev:
      space_over = land.endspace() + land_init.carspacevec[land_init.carvec.index(car)]
    land_init.exit_car(car)
    land.enter_car(car, space_over)
    self.sweep()

  def rep(self, dur_t):
    for land in self.landvec:
      land.repitar(dur_t)
    self.sweep()

  def is_route(self, land_array):
    ret = 1
    for i in range(len(land_array) - 1):
      la_now, la_next = land_array[i], land_array[i+1]
      b=1
      b *= la_next in la_now.outlandvec 
      b += la_next in la_now.inlandvec
      ret *= b 
    return ret

  def outroute(self, land_contact, land_array):
    land_contact.append(land_array[0])
    for i in range(len(land_array) - 1):
      land_array[i].outflow(land_array[i+1])
    self.sweep()  

  def inroute(self, land_contact, land_array, inpoint_array):
    land_arr, in_arr = copy(land_array), copy(inpoint_array)
    land_arr.reverse()
    in_arr.reverse()   
    land_contact.inflow(land_arr[0], in_arr[0])
    for i in range(len(land_arr) - 1):
      land_arr[i].inflow(land_arr[i+1, in_arr[i+1]])
    self.sweep()   

  def outpll(self, land_contact, land_array):
    for la in land_array:
      land_contact.outflow(la)
    self.sweep()    

  def inpll(self, land_contact, land_array, inpoint):
    for la in land_array:
      land_contact.inflow(la, inpoint)
    self.sweep()    

### TERRAIN - RETURNING FUNCTIONS
def Blankterrain():
  return Terrain(Blankland())

def Carterrain( car ):
  '''<Car object>
  Returns Terrain object 
  based on instance variables in <car>'''
  lactr = Uniformland(car.tilt, car.helm, MILE)
  lactr.constant_frix(car.frix_coeff)
  lactr.recompass(car.shape.orientation)
  lactr.land_position = car.shape.position
  terra = Terrain(lactr)
  lactr.enter_car(car, 0)
  terra.sweep()   
  return terra

def Fleurterrain(fork_helm):
  '''<helm> for the curved branches 
  Returrns Terrain object in the shape of a fleur-de-lis
  '''
  out_v, in_v = [],[]
  for i in range(3):
    h = fork_helm * (i - 1)
    b1, b2 = Blankland(), Blankland()
    b1.constant_helm(h)   
    b2.constant_helm(h)    
    out_v.append(b1)
    in_v.append(b2)
  ter = Blankterrain()   
  ter.outpll(ter.land_center, out_v) 
  ter.inpll(ter.land_center, in_v, 0)
  return ter
      


### DRIVER CLASS 
class Driver:

  def __init__(self, car, terrain):
    '''<Car object>, <Terrain object>
    '''
    self.name = "driver"
    self.car = car
    self.terrain = terrain
    self.land = self.terrain.landsearch(self.car)
    self.pwr_out = self.car.pwr_drive 
    #self.xfer_funxion = Funxion(unity, self.pwr_out)
    #self.pwr_in = self.xfer_funxion.out_val
    self.t_num = 10 #integer to be used as multiplying factor for time constants
    self.next_land = terrain.land_center 
    self.drive_funxion = self.car.drive_fxn()
    self.other_funxion = self.car.other_fxn()
    self.ndx_car = self.land.carvec.index(self.car)

  def sweep(self):
    #self.land = self.terrain.landsearch(self.car)
    self.car.pwr_drive = self.pwr_out
    space = self.land.carspacevec[self.ndx_car]
    if self.car.rev and (space < 0):
      if self.next_land != self.land:
        self.terrain.changeland(self.car, self.next_land)
    elif (not self.car.rev) and (space > self.land.endspace()):
      if  self.next_land != self.land:
        self.terrain.changeland(self.car, self.next_land)
    self.terrain.sweep
    self.land = self.terrain.landsearch(self.car)
    self.ndx_car = self.land.carvec.index(self.car)

  def is_behind(self, other_car):
    '''<Car object>
    Returns True if and only if this Driver's Car
    is Behind/Upstream from <other_car>'''
    ret = 0
    ndx_self = self.land.carvec.index(self.car)
    point_self = self.land.carspacevec[ndx_self]
    if other_car in self.land.carvec:
      v=self.land.sorted_cars()  
      ndx_o = v.index(other_car)
      ret += ndx_self < ndx_o
    for la in self.land.inlands_after(point_self):
      ret += other_car in la.carvec
    for la in self.land.outlandvec:
      ret += other_car in la.carvec
    return bool(ret)

  def is_infront(self, other_car):
    '''<Car object>
    Returns True if and only if Driver.car
    is in Front of / Downstream from <other_car>'''
    ret = 0
    ndx_self = self.land.carvec.index(self.car)
    point_self = self.land.carspacevec[ndx_self]
    if other_car in self.land.carvec:
      point_car = self.land.carspacevec[self.land.carvec.index(other_car)]
      ret += point_car < point_self
    for la in self.land.inlands_before(point_self):
      ret += other_car in la.carvec
    return bool(ret)

  def distance_to(self, other_car):
    '''<Car object>
    Returns bumper-to-bumper distance between this Driver's car and that Car,
    in Car.path units'''
    #print("\n\nDriver.distance_to engaged\n")
    ret = INF
    other_land = self.terrain.landsearch(other_car)
    point_self = self.land.carspacevec[ self.land.carvec.index(self.car) ]
    #print("point_self = "+str(point_self)+"\n")
    if self.is_behind(other_car):
      #print("self.is_behind( "+other_car.name+" )\n")
      ret = 0 - self.car.shape.ray_length(0)
      #print("ret = "+str(ret)+"\n") 
      if other_car in self.land.carvec:
        #print("other_car in self.land.carvec\n")
        pt_car = self.land.carspacevec[ self.land.carvec.index(other_car) ]
        #print("pt_car = "+str(pt_car)+"\n")
        pt_car -= other_car.shape.ray_length(PI)
        #print("pt_car = "+str(pt_car)+"\n")
        pt_car -= point_self
        #print("pt_car = "+str(pt_car)+"\n")
        ret += pt_car
        #print("ret += pt_car\nret = "+str(ret)+"\n")
      elif other_land in self.land.outlandvec:
        #print("other_land in self.land.outlandvec")
        pt_car = other_land.carspacevec[ other_land.carvec.index(other_car) ]
        #print("pt_car = "+str(pt_car)+"\n")
        pt_car += self.land.endspace() - point_self
        #print("pt_car = "+str(pt_car)+"\n")
        pt_car -= other_car.shape.ray_length(PI/2)
        #print("pt_car = "+str(pt_car)+"\n")
        ret += pt_car
        #print("ret += pt_car\nret = "+str(ret)+"\n")
      elif other_land in self.land.inlands_after(point_self):
        #print("other_land in inland_after\n")
        pt_car = other_land.endspace() 
        #print("pt_car = "+str(pt_car)+"\n")
        pt_car -= other_land.carspacevec[other_land.carvec.index(other_car)]
        #print("pt_car = "+str(pt_car)+"\n")
        pt_car += self.land.inspacevec[self.land.inlandvec.index(other_land)]
        #print("pt_car = "+str(pt_car)+"\n")
        pt_car -= point_self
        #print("pt_car = "+str(pt_car)+"\n")
        pt_car -= other_car.shape.ray_length(0)
        #print("pt_car = "+str(pt_car)+"\n")
        ret += pt_car
        #print("ret += pt_car\nret = "+str(ret)+"\n")
    elif self.is_infront(other_car):
      #print("self.is_infront( "+other_car.name+")\n")
      ret = 0 - self.car.shape.ray_length(PI)  
      #print("ret = "+str(ret)+"\n")
      if other_car in self.land.carvec :
        #print("other_car in self.land.carvec\n")
        pt_car = self.land.carspacevec[self.land.carvec.index(other_car)]
        #print("pt_car = "+str(pt_car)+"\n")
        ret += point_self - pt_car
        #print("ret += point_self - pt_car\nret = "+str(ret)+"\n")
        ret -= other_car.shape.ray_length(0)
        #print("ret = "+str(ret)+"\n")
      elif other_land in self.land.inlands_before(point_self):
        #print("other_land in self.land.inlands_before(point_self)\n")
        ret += point_self
        #print("ret = "+str(ret)+"\n")
        pt_car = other_land.carspacevec[other_land.carvec.index(other_car)]
        #print("pt_car = "+str(pt_car)+"\n")
        ret += other_land.endspace() - pt_car
        #print("other_land.endspace() = "+str(other_land.endspace())+"\n")
        #print("ret += other_land.endspace() - pt_car\n")
        #print("ret = "+str(ret)+"\n")
        ret -= other_car.shape.ray_length(0)
        #print("ret = "+str(ret)+"\n")
    #print("Final value: ret = "+str(ret)+"\n\n")
    return ret

  def t_collision(self, other_car):
    '''<Car object>
    Returns time, in seconds, 
    until Driver.car would collide with <car>
    at their present speeds'''
    print("\n\nDriver.t_collision engaged\n")
    distance = self.distance_to(other_car)
    print("Driver.distance_to( "+other_car.name+" )  =  "+str(distance)+"\n")
    point_self = self.land.carspacevec[self.land.carvec.index(self.car)]
    print("point_self = "+str(point_self)+"\n")
    speed_twd = 0
    my_speed = 0 + self.car.speed  
    yo_speed = 0 + other_car.speed
    sgn_self, sgn_other = -1,1
    is_beh = self.is_behind(other_car)
    casa = 0
    casa += int(other_car.rev)
    casa += 2 * int(self.car.rev)
    casa += 4 * int(is_beh)
    print("is_beh = "+str(int(is_beh))+"\n")
    print("self.car.rev = "+str(self.car.rev)+"\n")
    print("other_car.rev = "+str(other_car.rev)+"\n")
    print("casa = "+str(casa)+"\n")
    if casa in [2,3,4,5]:
      sgn_self *= -1
    if casa in [1,3,4,6]:
      sgn_other *= -1
    my_speed *= sgn_self
    yo_speed *= sgn_other
    print("my_speed = "+str(my_speed)+"\n")
    print("yo_speed = "+str(yo_speed)+"\n")
    speed_twd += my_speed + yo_speed
    if speed_twd < 0:
      speed_twd = 0
    print("speed_twd = "+str(speed_twd)+"\n")
    ret =  xdiv(0, distance, speed_twd)
    print("returning "+str(ret)+"\n\n\n")
    return ret
    






  def front_car(self):
    '''Returns Car object closest in front 
    of this Driver object's Car'''
    cv = self.land.sorted_cars()
    if cv.index(self.car) < (len(cv)-1):
      return cv[ 1 + cv.index(self.car)]
    #Add other if/else branches
    #to deal with nearby cars on inflows and outflows
    else:
      return self.car #default fir when there are no cars
    
  def back_car(self):
    cv = self.land.sorted_cars()   
    if cv.index(self.car) > 0:
      return cv[ cv.index(self.car) - 1]
    #Add other if/else branches 
    #for cars on nearby inflows and outflows
    else:
      return self.car #default answer for when no cars behind you

  def reverse(self):
    '''Encapsulation of Car.reverse for this Driver object'''
    self.car.reverse()    

  def forward(self):
    '''Encapsulation of Car.forward for this Driver object'''
    self.car.forward()

  def position(self):
    '''Returns position phasor returning centerpointof this Driver object's Car'''
    return self.car.shape.position



  def pwr_equilibrium(self):
    '''Returns mechanical power necessary to counteract both gravity and friction'''
    return -1*(self.car.pwr_grav + self.car.pwr_other)

  def steady_fxn(self):
    '''Driver.pwr_equilibrium encapsulated in Funxion object'''
    def dummy(t):
      return self.pwr_equilibrium()
    return Funxion(dummy, self.car.t_now)

  def slow_fxn(self, speed, dur_accel):
    delta_kin = kinetic(self.car.mass, speed) - self.car.energy
    watt_xtra = float(delta_kin) / dur_accel
    def dummy(t):
      ret = 0
      if self.car.speed > speed:
        ret = watt_xtra
      return ret
    return Funxion(dummy, self.terrain.t_terrain)

  def fast_fxn(self, speed, dur_accel):
    delta_kin = kinetic(self.car.mass, speed) - self.car.energy
    watt_xtra = float(delta_kin) / dur_accel
    def dummy(t):
      ret = 0
      if self.car.speed < speed:
        ret = watt_xtra
      return ret
    return Funxion(dummy, self.terrain.t_life)


  def tospeed_fxn(self, speed, dur_accel):
    if speed < self.car.speed:
      return SumFunxion( [self.slow_fxn(speed, dur_accel), self.steady_fxn()], self.terrain.t_life)
    else:
      return SumFunxion([self.fast_fxn(speed, dur_accel), self.steady_fxn()] , self.terrain.t_life)


  def choose_next(self, downland):
    self.next_land = downland

  def choose_ang_out(self, ang):
    thv = copy(self.land.outangvec)
    lav = copy(self.land.outlandvec)
    if len(lav) == 0:
      self.choose_next(self.land)
    else:
      th_close = closest(ang, thv)
      ndx_land = thv.index(th_close)
      #ndx_v = val_index_vec(th_close, thv)
      if num_in_list(th_close, thv) > 1:
        ndx_v = val_index_vec(th_close, thv)
        helmv = []
        for ndx in ndx_v:
          helmv.append(lav[ndx].helmvec[0])
        th_close = closest(ang, helmv)
        ndx_land = helmv.index(th_close)
      l_next = lav[ndx_land]
      self.choose_next(l_next)


  def choose_straight(self):
    self.choose_ang(0)

  def choose_left(self):
    self.choose_ang(-1*PI/2)

  def choose_right(self):
    self.choose_ang(PI/2)

  def adjust_pwr(self):
    #self.other_funxion(self.car.energy)
    self.car.pwr_other = self.other_funxion.funx(self.car.energy)
    #self.drive_funx(self.terrain.t_life)
    self.pwr_out = self.drive_funxion.funx(self.terrain.t_life)
    #self.car.pwr_drive = self.pwr_out
    self.sweep()

  def rep_driver(self, dur_t):
    self.adjust_pwr()    
    self.terrain.rep(dur_t)
    self.sweep()      

  def rep_short(self):
    self.rep_driver(self.t_num * self.terrain.t_terrain)

  def rep_long(self, dur_sim):
    t_unit = self.t_num * self.terrain.t_terrain
    n = int(dur_sim // t_unit)
    for i in range(n):
      self.rep_short()


### DRIVER - RETURNING FUNCTIONS

def Terraindriver(terra):
  '''<Terrain object>
  Returns default Driver object
  made to drive in that Terrain'''
  car = Car(1000, 0, 0)
  terra.land_center.enter_car(car, 0)
  terra.sweep()
  dr = Driver(car, terra)
  return dr

def Cardriver(car):
  '''<Car object>
  Returns default Driver object
  meant to drive that Car'''
  lactr = Blankland()
  lactr.land_position = car.shape.position
  lactr.enter_car(car, 0) 
  terra = Terrain(lactr)
  terra.sweep()   
  return Driver(car, terra)


  

    


### TESTING SECTION

print("Iguana module running\n\n")

ctr, land_left, land_right = Blankland(), Blankland(), Blankland()
ctr.rename("ctr")
land_left.rename("land_left")
land_left.recompass(O_CLOCK ** 4)
land_right.rename("land_right")
land_right.recompass(O_CLOCK ** 2)
nc = named_cars(["camry", "cam_fwd", "cam_bk"])
camry, cam_fwd, cam_bk = nc[0], nc[1], nc[2]
pt_cam = MILE/2

ctr.enter_car(camry, pt_cam)
ctr.enter_car(cam_bk, pt_cam - 105)
camry.respeed(35)
cam_bk.respeed(35)



ter=Terrain(ctr)
my_rev = [0,0,1,1]
yo_rev = [0,1,0,1]

adam = Driver(camry, ter)
strega = ""
for i in range(4):
  if(my_rev[i]):
    camry.reverse()   
  else:
    camry.forward()    
  if(yo_rev[i]):
    cam_bk.reverse()   
  else:
    cam_bk.forward()   
  strega += "adam.t_collision( cam_bk) = "+str(adam.t_collision(cam_bk))+"\n"
print(strega)

#print(strega)





