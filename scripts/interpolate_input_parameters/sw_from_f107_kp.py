import math

def hpi_from_gw(gw):
  if gw <= 2.5:
    return '1'
  elif gw <= 3.94:
    return '2'
  elif gw <= 6.22:
    return '3'
  elif gw <= 9.82:
    return '4'
  elif gw <= 15.49:
    return '5'
  elif gw <= 24.44:
    return '6'
  elif gw <= 38.56:
    return '7'
  elif gw <= 60.85:
    return '8'
  elif gw <= 96.0:
    return '9'
  else:
    return '10'

def swbt_calc(swbz,swby):
  return math.sqrt(swbz**2+swby**2)

def swden_calc(): # this isn't used?
  return 5.0

def swvel_calc(kp):
  return 317.0+55.84*kp-2.71*kp**2

def swang_calc(by,bz):
  return math.atan2(by,bz)/math.pi*180

def swby_calc():
  return 0.0

def swbz_calc(kp,f107):
  return -0.085*kp**2 - 0.08104*kp + 0.4337 + 0.00794 * f107 - 0.00219 * kp * f107

def hemi_pow_calc(kp):
  return 1.29 + 15.60*kp - 4.93*kp**2 + 0.64*kp**3

def calc_solar_data(kp,f107):
  swbt     = []
  swangle  = []
  swvel    = []
  swbz     = []
  hemi_pow = []
  hemi_pow_idx = []
  for i in range(len(f107)):
    swbz.append(    swbz_calc(kp[i],f107[i]))
    swbt.append(    swbt_calc(swbz[i],swby_calc()))
    hemi_pow.append(hemi_pow_calc(kp[i]))
    swangle.append( swang_calc(swby_calc(),swbz[i]))
    swvel.append(   swvel_calc(kp[i]))
    hemi_pow_idx.append(hpi_from_gw(hemi_pow[i]))
  return swbt,swangle,swvel,swbz,hemi_pow,hemi_pow_idx


