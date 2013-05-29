import numpy as np
import argparse
import math
from pylab import *
import matplotlib.pyplot as plt


# Furty 
c = 299792458 #
deg2rad = 180.0/np.pi # degrees to radians conversion
r = 6371009.0      # Earth radius in meters

# site geographic coordinates in degrees
# UTEF Praha
#long1o = 14.4229437
#lat1o = 50.0673815
#long1 = long1o/deg2rad
#lat1 = lat1o/deg2rad

#Edmonton 
long1o = -113.523
lat1o = 53.47949
long1 = long1o/deg2rad
lat1 = lat1o/deg2rad



# detector coordinates for 
#x1 = -6.403
#y1 = -7.907
#z1 = 0
#x2 = 3.602
#y2 = -9.286
#z2 = 0

x0 = 11.549
y0 = 0.341
z0 = 0
x1 = 1.499
y1 = 0.253
z1 = 0
x2 = 6.566
y2 = 8.971
z2 = 0





def get_data(path):
  print "Using ",path
  data = []
  f=open(path,"r")
  for line in f.readlines():
    evn=line.split()
    for i in range(1,len(evn)):
      evn[i]=float(evn[i])
    data.append(evn)

  f.close()
  return(data)
 
def get_jd(data):
  y=data[1]
  m=data[2]
  d=data[3]
  hr=data[4]
  mi=data[5]
  sec=data[6]
  if m < 3:
    f = m + 12
    g = y - 1
  else:
    f = m
    g = y
  JD_date = d + trunc((153*f - 457)/5.0) + 365*g + floor(g/4.0) - floor(g/100.0) + floor(g/400.0) + 1721118.5
  JD_rest = hr/24.0 + mi/(24.0*60) + sec/(24.0*3600)
  JD = JD_date + JD_rest
  return (JD)

def get_lst(jd):
  lst = 6.697374558 + 2400.051337*((jd - 2451545)/36525.0) + 24*(jd + 0.5 - trunc(jd + 0.5)) + long1o/15.0
  lst = lst%24.0
  return(lst)

def get_dec(az,alt,lat1):
  dec=arcsin(sin(alt)*sin(lat1) - cos(alt)*cos(az)*cos(lat1))
  return(dec)

def get_h_angl(az,alt,lat1,dec):
    cos_hang = (sin(alt)*cos(lat1) + cos(az)*cos(alt)*sin(lat1))/cos(dec)
    sin_hang = (cos(alt)*sin(az))/cos(dec)
    if cos_hang < 0:
        hang = pi - arcsin(sin_hang)
    else:
        hang = arcsin(sin_hang)
    return(hang)

def get_ra(lst,h_angl):
  ra=(15.0*lst/deg2rad) - h_angl
  return(ra)

def get_gal_lat(ra,dec):
  gal_lat = arcsin(cos(dec)*cos((27.4/deg2rad))*cos(ra - (192.25/deg2rad)) + sin(dec)*sin((27.4/deg2rad)))
  return(gal_lat)

def get_gal_long(ra,dec,gal_lat):
  gal_long = math.atan2((sin(dec) - sin(gal_lat)*sin((27.4)/deg2rad)), (cos(dec)*sin(ra - (192.25/deg2rad))*cos((27.4/deg2rad)))) + (33.0/deg2rad)
  return(gal_long)

def get_coord(data):
  gal_lat_lst=[]
  gal_long_lst=[]
  source = []
  azim = []
  altit = []
  zaltit=[]
  ra_lst=[]
  dec_lst=[]
  for i in range(1,len(data)):
    if data[i][0] == "a":
     # 1=> (B - A) | 2=> (C - A)
 
     d1 = c*(data[i][9]-data[i][8])*25e-12 #s
     d2 = c*(data[i][10]-data[i][8])*25e-12  #s


     A = c*(d2 - d1*((y2-y0)/(y1-y0)))/((x2 - x0) - (x1 - x0)*((y2-y0)/(y1-y0)))
     B = c*(d2 - d1*((x2-x0)/(x1-x0)))/((y2 - y0) - (y1 - y0)*((x2-x0)/(x1-x0)))
     az = np.arctan2(A,B) #* deg2rad
     azz = (az*deg2rad)+180.0 
       


     alt = np.arccos(np.sqrt((A*A + B*B)/(c**2)))# * deg2rad
     zalt=90-(alt*deg2rad)  
     if math.isnan(az)==False and math.isnan(alt)==False:
       jd = get_jd(data[i])
       lst = get_lst(jd) 
       dec = get_dec(az,alt,lat1) 
       h_angl = get_h_angl(az,alt,lat1,dec)
       ra = get_ra(lst,h_angl) 
       gal_lat=get_gal_lat(ra,dec)
       gal_long=get_gal_long(ra,dec,gal_lat)

       coord=(az,alt,jd,lst,dec,h_angl,ra,gal_lat,gal_long) # IN RADIANS ALL  
       source.append(coord)
       azim.append(azz)
       zaltit.append(zalt)

       gal_lat_lst.append(gal_lat*deg2rad)
       gal_long_lst.append(gal_long*deg2rad)  
       ra_lst.append(ra*deg2rad)
       dec_lst.append(dec*deg2rad) 
       #print az,alt,jd,lst,dec,h_angl,ra, gal_lat,gal_long
      # print ra,dec,gal_lat, gal_long
  return(source,azim,zaltit,ra_lst,dec_lst,gal_lat_lst,gal_long_lst)




def altaz_plot(az,zalt):
#  print len(az),len(zalt), min(az), max(az), min(alt),max(alt)
  plt.figure(figsize=(10,10))
  plt.subplot(111, polar=True)
  plt.plot(az,zalt,',' ,markersize=1)
#  plt.plot.set_yticks(range(90, 0, 10))
 # ax.set_rmax(100.0)
  plt.grid(True)
 # plt.show()
 # plt.savefig("out_file.png",dpi=60)
  return()

def althist_plot(alt,glat):
  fig = plt.figure(figsize=(20,10))
  ax=fig.add_subplot(121)
  ax.hist(alt,45,color="green",alpha=0.75)
  ax.set_xlabel('Zenit distance [deg]')
  ax.set_ylabel('Event Count')
  ax.grid(True)
# plt.savefig("out_file1.png",dpi=60)
 # plt.show()
  ax1=fig.add_subplot(122)
  ax1.hist(glat,45,color="blue",alpha=0.75)
  ax1.set_xlabel('Galactic latitude [deg]')
  ax1.set_ylabel('Event Count')
  ax1.grid(True)



  return()

def gal_coord_plot(glat,glong):
  fig = plt.figure(figsize=(20,10))
  ax=fig.add_subplot(111,projection='mollweide')
  ax.set_xlabel('long')
  ax.set_ylabel('lat')
 # ax.set_xticklabels(np.arange(30,331,30))
  hist,xedges,yedges = np.histogram2d(glat,glong,bins=[100,100],range=[[-90,90],[-180,180]])
  X,Y = np.meshgrid(np.radians(yedges),np.radians(xedges))
  image = ax.pcolormesh(X,Y,hist,shading='flat')
  cb = fig.colorbar(image, orientation='horizontal')
  ax.grid(True)
  

  return()


if __name__ == "__main__":

  parser = argparse.ArgumentParser(description="Czelta event position calculator")
  parser.add_argument('--data',type=str, nargs=1, required=True, action='store',
                     help='data for analysis')
  parser.add_argument('--site',type=str,nargs=1, required=True, action='store',
                     help='site name, suported: laurent (will be more in futere)')
  args = parser.parse_args()
  path = args.data[0]
  site = args.site[0]
  if site=='laurent' :
    ev=get_data(path)
    print "  "
    print "data loaded from laurent site"
    res,az,alt,ra_l,dec_l,gal_lat_lst,gal_long_lst=get_coord(ev)  
    altaz_plot(az,alt)
    althist_plot(alt,gal_lat_lst)
    gal_coord_plot(gal_lat_lst,gal_long_lst)
    plt.show()











