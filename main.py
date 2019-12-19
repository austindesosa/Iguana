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

INCH=0.0254
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

O_CLOCK = eul(30 * DEG)
#People sometimes say "o'clock" as a unit of helm
#for example, "1 o'clock" means 30 degrees rightward
#O_CLOCK is that unit as a direction phasor

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

### PHYSICS MATH FUNCTIONS

#WARNING: Doesn't work for perfectly straight dlines
#only works for curved paths
def curve_ray(theta, kurv_r):
  '''<Angle in radians> subtended by a curve, <Curvature Radius> of that curve
  Returns phasor (complex number) 
  whose Magnitude is the distance from startpoint to endpoint of the curve as the crow flies,
  and whose Angle is the angle of that ray wrt driver's direction 
  at beginning of driving along that curve'''
  mag_ang = 0.5 * (PI - theta)
  mag = 2 * kurv_r * cos(mag_ang)
  return mag * eul(theta / 2)

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



  def tostring(self):
    strega="Funxion   "+self.name+"  has these instance variables:\n"
    strega+="self.in_val = "+str(self.in_val)+"\nself.out_val = "+str(self.out_val)
    strega+="\nEND Funxion   "+self.name+"   instance variables\n\n"
    return strega


class Car:

  def __init__(self, mass, pwr_drive, energy):
    '''<Mass of car>, <Intentional power on car>, <Kinetic energy of car>'''
    self.name = "car" #default name for this Car object
    self.mass, self.pwr_drive, self.energy = mass, pwr_drive, energy
    self.pos, self.delta_pos = 0*eul(0), 0*eul(0) #Position, changge in position, sored as a phasor
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
    #EXPERIMENTAL INSTANCE VARIABLES
    self.is_curved = abs(self.helm) < (0.000001 * DEG) #Boolean true if car's path is curved
    self.delta_tau = 0.1 #Short time delay for actual time delays like time.sleep()
    


  
  def sweep(self):
    self.speed = (2.0 * self.energy / self.mass) ** (0.5)
    self.brake = self.pwr_drive < 0 #Set brake flag if driver intentionally reducing energy
    self.pwr_grav = GRAV * sin(self.tilt) * self.mass * self.speed
    self.is_curved = abs(self.helm) < (0.000001 * DEG)
    if (self.speed <= 0) and (sgn(self.pwr_drive) < 0) :
      self.stop = 1
      self.energy = 0.0
    if (self.pwr_drive > 0 ):
      self.stop = 0
      self.brake = 0
    #CODE NOTE: THe instance variables self.path and self.delta_path are NOT affected
    #by the sweep() method, because sweep method does not handle the motion itself
    #and also because recalculating the displacement too often can lead to inaccurate results

  def go(self, delta_t):
    '''<Short amount of time> in secondss
    Changes the instance variables to model a small amount of motion in the car'''
    energy_initial = 0.0 + self.energy
    d_t = 0.0 + delta_t
    if (self.pwr_net < 0) and ( abs(self.pwr_net * delta_t) > energy_initial):
      d_t = energy_initial / self.pwr_net #amount of time it will take car to stop at this rate
    time.sleep(self.delta_tau) #experimental code, delay for threading purposes
    self.energy += self.pwr_net * d_t
    self.sweep()
    energy_avg = avg([energy_initial, self.energy]) #average energy duting the short time, 
    #assuming constant mechanical power throughout
    speed_avg = (2.0 * energy_avg / self.mass) ** 0.5 #average speed duting the time elapsed
    self.delta_path = speed_avg * d_t
    #Mostly experimental/untested code the rest of the method
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
    #the rest of the code in this method is for testing purposes only

    
  def travel_raw(self, dur, delta_t):
    '''<Duration> in seconds, <Polling Period> in seconds
    Iterates the method Car.go(delta_t) until <dur> seconds have supposedly elapsed'''
    n = dur // delta_t
    for i in range(n):
      self.go(delta_t)

  def travel(self, dur, delta_t):
    '''<Duration> in seconds, <Polling Period> in seconds
    Iterates the method Car.go(float), but in a threaded manner,
    so other code can execute at the same time'''
    thr = threading.Thread(target = self.travel_raw, name = "", args = (dur, delta_t))
    thr.start()

  def str_motion(self):
    '''Returns string with attributes directly relevant to the motion'''
    strega = ""
    strega+="self.speed = "+str(self.path)+"  meters per second\n"
    strega+="self.pwr_drive = "+str(self.path)+"  watts\n"
    strega+="self.pwr_net = "+(self.path)+"  watts\n"
    strega+="self.energy = "+(self.path)+"  watt-seconds\n"
    strega+="self.path = "+(self.path)+"  meters \n"
    strega+="self.delta_path = "+(self.delta_path)+"  meters\n"
    strega+="self.pos = "+(self.pos)+"  meters\n"
    strega+="self.delta_pos = "+(self.delta_pos)+"  meters\n"
    strega+="\n\n\n"
    return strega

  def tell_motion_raw(self, dur, d_tau):
    #UNTESTED
    '''<Duration> in seconds you want to run this method, <Polling Period> in seconds between print statements
    Prints attrributes related to the motion of the car every <d_tau> seconds, 
    until ther have been <'''
    n = dur // d_tau
    for i in range(n):
      print(self.str_motion)
      self.sleep(d_tau)

  def tell_motion(self, dur, d_tau):
    #UNTESTED
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
    return car

  def copy_name(self, new_name):
    '''<String> representing name of new object
    Copies this Car object and gives resulting object a new name'''
    car = self.copy() 
    car.rename(new_name)
    return car

  
    
  
    

    
  
  def tostring(self):
    strega = "Car   "+self.name+"   has these attributes:\n"
    strega += "self.energy = "+str(self.energy)+"  watt-secondss\n"
    strega += "self.pwr_drive = "+str(self.pwr_drive)+"  watts\n"
    strega += "self.pwr_grav = "+str(self.pwr_grav)+"  watts\n"
    strega += "self.pwr_other = "+str(self.pwr_other)+"  watts\n"
    strega += "self.pwr_net = "+str(self.pwr_net)+"  watts\n"
    strega += "self.tilt = "+str(self.tilt)+"  radians downhill\n"
    strega += "self.mass = "+str(self.mass)+"  kilograms\n"
    strega += "self.stop = "+str(self.stop)+"  boolean\n"
    strega += "self.brake = "+str(self.brake)+"  boolean\n"
    strega += "self.rev = "+str(self.rev)+"  boolean\n"
    strega += "self.is_curved"+str(self.is_curved)+"  boolean\n"
    strega += "self.helm = "+str(self.helm)+"  radians per meter rightward\n"
    strega += "self.path = "+str(self.path)+"  meters\n"
    strega += "self.pos = "+str(self.pos)+" meters\n"
    strega += "self.speed = "+str(self.speed)+"  meters per second\n"
    strega += "self.delta_path = "+str(self.delta_path)+"  meters\n"
    strega += "self.delta_pos = "+str(self.delta_path)+"  meters\n"
    strega += "self.delta_tau = "+str(self.delta_path)+"  seconds\n"
    strega += "END Car   "+self.name+"   attributes\n\n\n"
    return strega


### TESTING SECTION







