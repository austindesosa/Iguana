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

MIN_HELM = 0.000001 * DEG # curves WITH LESS HELM THAN THIS
# can be considered straight lines

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
  '''《lHospial limit》, 《scalar》, 《other scalar》
  lhop is winner of l'Hospital's
  rule
  Multiplies the two POSSIBLY INFINITE scalars 《x1》 and 《x2》 in accoRdance with whatever lHospital's rule
  says 0*inf would be (see MATH NOTE above)
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

  def scale(self, k):
    '''<Scaling factor>
    Returns Funxion object just like this one,
    except the encapsulated function is scaled by <k>'''
    def dummy(x):
      return k * self.funx(x)
    fu = Funxion(dummy, self.in_val)
    return fu




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
      funxion.feed( x )
      yvec.append(funxion.out_val)
    return sum(yvec)
  fu = Funxion(dummy, 0)
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
      funxion.feed( x )
      yvec.append(funxion.out_val)
    return prod(yvec)
  fu = Funxion(dummy, 0)
  return fu


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
    self.delta_tau = 0.0000000001 #Short time delay for actual time delays like time.sleep()
    self.t_poll = 0.1 #Short time delay for predicting the car‘s motion
    self.t_now = 0.0 #Amount of time car has been in motion
    self.shape = Rectangle(5,3)
    
    
    


  
  def sweep(self):
    self.speed = (2.0 * self.energy / self.mass) ** (0.5)
    self.brake = self.pwr_drive < 0 #Set brake flag if driver intentionally reducing energy
    self.pwr_grav = GRAV * sin(self.tilt) * self.mass * self.speed
    self.pwr_net = self.pwr_grav + self.pwr_drive + self.pwr_other
    self.is_curved = abs(self.helm) > (0.000001 * DEG)
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
    #EXPERIMENTAL CODE 
    self.t_now += d_t
    #END EXPERIMENTAL CODE
    self.sweep()
    energy_avg = avg([energy_initial, self.energy]) #average energy duting the short time, 
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
    else:
      self.path -= self.delta_path
      self.pos -= self.delta_pos
    

  def go_drive(self, delta_t, drivefunxion):
    '''<Amount of time> in seconds for motion, <Funxion objet> representing self.pwr_drive as a  function of time
    Uses self.go( <delta_t> ) but decides self.pwr_driveby the function encapsulated in <drivefunxion>'''
    drivefunxion.feed(self.t_now)
    self.pwr_drive = drivefunxion.out_val
    self.sweep() 
    self.go(delta_t)
    
  def go_other(self, drivefunxion, otherfunxion):
    #@param   drivefunxion   is Funxion object
    #representing self.pwr_drive as a function of Time
    #@param   otherfunxion   is Funxion object
    #representing self.pwr_other as a function of Energy   
    otherfunxion.feed(self.energy)
    self.pwr_other = otherfunxion.out_val
    self.sweep()
    self.go_drive(self.t_poll, drivefunxion)

    
  def travel_raw(self, dur, delta_t):
    '''<Duration> in seconds, <Polling Period> in seconds
    Iterates the method Car.go(delta_t) until <dur> seconds have supposedly elapsed'''
    n = int(dur // delta_t)
    for i in range(n):
      self.go(delta_t)

  def travel(self, dur, delta_t):
    '''<Duration> in seconds, <Polling Period> in seconds
    Iterates the method Car.go(float), but in a threaded manner,
    so other code can execute at the same time'''
    thr = threading.Thread(target = self.travel_raw, name = "", args = (dur, delta_t))
    thr.start()
    
  def to_speed_lin(self, speed_des, accel_time):
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
    #UNTESTED
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
    strega += "END Car   "+self.name+"   attributes\n\n\n"
    return strega



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
    self.compass = NORTH #orientation of the land AKA driver’s direction upon entering, as a unit phasor, assumed NORTH by default
    self.carvec = [] #list of Car objects inteacting with this Land object
    #assumed empty by default
    self.phasorvec = []
    self.carspacevec = []
    for i in range(self.size):
      ang = self.helmvec[i]*self.delta_space
      if abs(ang) > MIN_HELM:
        self.phasorvec.append(self.delta_space * eul(0))
      else:
        curv_r = 1.0 / abs(self.helmvec[i])
        self.phasorvec.append(curve_ray(ang, curv_r))
    self.ray = sum(self.phasorvec) #Phasor representing net displacement of path, in same units as Car.pos
    self.tau_max = .00001 #
    self.t_land = .01 #Polling rate, in simulated seconds, for Car objects in this Land 
    
    
  def ndx_point(self, path):
    '''<Distance travelled> along path represented by this Land object, in same units as Car.path
    Returns index of self.tiltvec and self.helmvec corresponding to that point along the path
    '''
    n = int(path // self.delta_space)
    n = min(n, self.size)
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
      self.phasorvec.pop[-1]
    self.ray = sum(self.phasorvec)

  def point_ray(self, space):
    '''<Distance travelled> along path represnted by theis Land object
    Returns phasor representing net displacement between beginning of path
    and that point in space''' 
    ndx = self.ndx_point(space)
    ray = 0 * eul(0)
    for i in range(ndx):
      ray += self.phasorvec[i]
    last_space = space - (ndx * self.delta_space)
    h = self.helmvec[ndx]
    if (h < MIN_HELM):
      ray += last_space * eul(0)
    else:
      ray += curve_ray(last_space, last_space * h)
    return ray 

    
  def enter_car(self, car, space):
    '''<Car object>, <Place on path> in Car.path units
    Puts a Car object to interact with this Land object'''
    self.carvec.append(car)
    self.carspacevec.append(space )
    self.sweep()

    
  def exit_car(self, car):
    '''<Car object>
    Removes that Car from this Land object if it's in here'''
    if car in self.carvec:
      ndx = self.carvec.ndx( car  )
      
    
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
    print("Land.ir_other_raw invoked\n")
    car = self.carvec[ndx_car]
    self.sweep() 
    car.go_other(drivefunxion, otherfunxion)
    signum = 1
    if car.rev:
      signum *= -1
    car.path += signum * car.delta_path 

  def ir_other(self, drivefunxion, otherfunxion, ndx_car):
    thr = threading.Thread(target = self.ir_other_raw, name="", args = (drivefunxion, otherfunxion, ndx_car))
    thr.start()
    


  
  def tostring(self):
    strega = "Land   "+self.name+"   has these instance variables \n"
    strega += "self.delta_space = "+str(self.delta_space)+"\n"
    strega += "self.tilt_avg = "+str(self.tilt_avg)+"\n"
    strega += "self.tiltvec = "+str(self.tiltvec)+"\n"
    strega += "self.helm_avg = "+str(self.helm_avg)+"\n"
    strega += "self.helmvec = "+str(self.helmvec)+"\n"
    strega += "self.ray = "+str(self.ray)+"\n"
    strega += "self.phasorvec = "+str(self.phasorvec)+"\n"
    strega += "self.compass = "+str(self.compass)+"\n"
    strega += "self.carvec = "+str(self.carvec)+"\n"
    strega += "self.tau_max = "+str(self.tau_max)+"\n"
    strega += "self.t_land = "+str(self.t_land)+"\n"
    strega += "END Land   "+self.name+"   instance variables\n\n\n"
    return strega
        

### TESTING SECTION

#Instantiate Car object
camry_mass = 1000.0
camry_energy = (camry_mass / 2)*((30*MPH)**2)
camry_pwr = camry_energy * 0.15
camry = Car(camry_mass, camry_pwr, camry_energy)

camry.t_poll=0.2
camry.sweep()

#Declare Funxion objects
manejar = Linear(0.3*camry_pwr, camry_pwr)
otro = Ramp(-1 * camry_mass / 10)

print(camry.str_power)
for i in range(7):
  camry.go_other(manejar, otro)
  print(camry.str_power()) 
  







