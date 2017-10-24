import numpy as np
from math import cos, sin, atan2, pi
from random import random
# from datetime import datetime
import datetime
def fcc(xlo, xhi, ylo, yhi, zlo, zhi, lattice):
    '''numpy array of fcc'''
    num = (xhi-xlo)*(yhi-ylo)*(zhi-zlo)
    data = np.zeros([num*4,3])
    index = 0
    cell = np.array([[0, 0, 0], [0.5*lattice, 0.5*lattice, 0], [0.5*lattice, 0, 0.5*lattice], [0, 0.5*lattice, 0.5*lattice]])
    for i in range(xlo, xhi):
        for j in range(ylo, yhi):
            for k in range(zlo, zhi):
                data[index:index+4, 0] = cell[:,0] + i*lattice
                data[index:index+4, 1] = cell[:,1] + j*lattice
                data[index:index+4, 2] = cell[:,2] + k*lattice
                index += 4
    return data

def center(xlo, xhi, ylo, yhi, zlo, zhi, lattice):
    '''numpy array of fcc'''
    num = (xhi-xlo)*(yhi-ylo)*(zhi-zlo)
    data = np.zeros([num,3])
    index = 0
    cell = np.array([0.5*lattice, 0.5*lattice, 0.5*lattice])
    for i in range(xlo, xhi):
        for j in range(ylo, yhi):
            for k in range(zlo, zhi):
                data[index, 0] = cell[0] + i*lattice
                data[index, 1] = cell[1] + j*lattice
                data[index, 2] = cell[2] + k*lattice
                index += 1
    return data


data = fcc(0,1,0,1,0,1,2)

def write_lammps(filename, data, lattice):
    '''write the simplest type of the data file'''
    xlo = min(data[:,0])
    xhi = max(data[:,0]) + 0.5*lattice
    ylo = min(data[:,1])
    yhi = max(data[:,1]) + 0.5*lattice
    zlo = min(data[:,2])
    zhi = max(data[:,2]) + 0.5*lattice

    filename = filename + '.data'
    file = open(filename, 'w')

    now = datetime.datetime.now()
    file.write('lammpsdata made by wzh in ' + now.strftime('%Y-%m-%d')  + '\n')

    file.write('%d\t %s\n' %(data.shape[0],'atoms'))
    atom_type = 1
    file.write('%d\t %s\n\n' %(atom_type, 'atom types'))

    file.write('%f %f \t%s\n' %(xlo, xhi, 'xlo xhi'))
    file.write('%f %f \t%s\n' %(ylo, yhi, 'ylo yhi'))
    file.write('%f %f \t%s\n\n' %(zlo, zhi, 'zlo zhi'))

    mass = [63.5]
    file.write('%s\n\n'%'Masses')
    file.write('%d %f\t\n'% (1, mass[0]))

    file.write('\n%s\n\n'%'Atoms')


    index = 1

    for i in range(data.shape[0]):
        file.write('%d %d %d %f %f %f %f \n'%(index, index, 1, 0, data[i,0], data[i,1], data[i,2]))
        index +=1

    file.close()

def write_lammps_common(filename, data, lattice):
    '''write the simplest type of the data file'''
    xlo = min(data[:,0])
    xhi = max(data[:,0]) + 0.5*lattice
    ylo = min(data[:,1])
    yhi = max(data[:,1]) + 0.5*lattice
    zlo = min(data[:,2])
    zhi = max(data[:,2]) + 0.5*lattice

    filename = filename + '.data'
    file = open(filename, 'w')

    now = datetime.datetime.now()
    file.write('lammpsdata made by wzh in ' + now.strftime('%Y-%m-%d')  +'\n')

    file.write('%d\t %s\n' %(data.shape[0],'atoms'))
    atom_type = 1
    file.write('%d\t %s\n\n' %(len(np.unique(data[:,3])), 'atom types'))

    file.write('%f %f \t%s\n' %(xlo, xhi, 'xlo xhi'))
    file.write('%f %f \t%s\n' %(ylo, yhi, 'ylo yhi'))
    file.write('%f %f \t%s\n\n' %(zlo, zhi, 'zlo zhi'))

    mass = [4.003, 63.5]
    file.write('%s\n\n'%'Masses')
    file.write('%d %f\t\n'% (1, mass[0]))
    file.write('%d %f\t\n'% (1, mass[0]))

    file.write('\n%s\n\n'%'Atoms')


    index = 1

    for i in range(data.shape[0]):
        file.write('%d %d %d %f %f %f %f \n'%(index, index, data[i,3], 0, data[i,0], data[i,1], data[i,2]))
        index +=1

    file.close()

def write_lammps_eff_he(filename, data, lattice, k=1, r=1):
    '''write the simplest type of the data file'''
    eradius_he = 0.6
    xlo = min(data[:,0])
    xhi = max(data[:,0]) + 0.5*lattice
    ylo = min(data[:,1])
    yhi = max(data[:,1]) + 0.5*lattice
    zlo = min(data[:,2])
    zhi = max(data[:,2]) + 0.5*lattice

    filename = filename + '_eff.data'
    file = open(filename, 'w')

    now = datetime.datetime.now()
    file.write('lammpsdata made by wzh in ' + now.strftime('%Y-%m-%d')  + 'k = %.2f lat*r = %f' %(k, lat*r) +'\n')

    file.write('%d\t %s\n' %(data.shape[0]*3,'atoms'))
    atom_type = 2
    file.write('%d\t %s\n\n' %(atom_type, 'atom types'))

    file.write('%f %f \t%s\n' %(xlo, xhi, 'xlo xhi'))
    file.write('%f %f \t%s\n' %(ylo, yhi, 'ylo yhi'))
    file.write('%f %f \t%s\n\n' %(zlo, zhi, 'zlo zhi'))

    mass = [4.003, 1.]
    file.write('%s\n\n'%'Masses')
    file.write('%d %f\t\n'% (1, mass[0]))
    file.write('%d %f\t\n'% (2, mass[1]))

    file.write('\n%s\n\n'%'Atoms')


    index = 1

    for i in range(data.shape[0]):
        #index type q spin eradius x y z
        file.write('%d %d %f %d %f %f %f %f \n'%(index, 1, 2.0, 0, 0., data[i,0], data[i,1], data[i,2]))
        index +=1

    for i in range(data.shape[0]):
        file.write('%d %d %f %d %f %f %f %f \n'%(index, 2, 0, 1, eradius_he, data[i,0], data[i,1], data[i,2]))
        file.write('%d %d %f %d %f %f %f %f \n'%(index+1, 2, 0, -1, eradius_he, data[i,0], data[i,1], data[i,2]))
        index +=2
    file.close()

def write_lammps_eff_h2(filename, data, lattice, k, r):
    '''write the simplest type of the data file'''
    eradius_h = 0.98
    xlo = min(data[:,0])-1
    xhi = max(data[:,0]) + 0.5*lattice+1
    ylo = min(data[:,1])-1
    yhi = max(data[:,1]) + 0.5*lattice+1
    zlo = min(data[:,2])
    zhi = max(data[:,2]) + 0.5*lattice

    filename = filename + '_eff.data'
    file = open(filename, 'w')

    now = datetime.datetime.now()
    file.write('lammpsdata made by wzh in ' + now.strftime('%Y-%m-%d')  + 'k = %.2f lat*r = %f' %(k, lat*r) +'\n')

    file.write('%d\t %s\n' %(data.shape[0]*4,'atoms'))
    atom_type = 2
    file.write('%d\t %s\n\n' %(atom_type, 'atom types'))

    file.write('%f %f \t%s\n' %(xlo, xhi, 'xlo xhi'))
    file.write('%f %f \t%s\n' %(ylo, yhi, 'ylo yhi'))
    file.write('%f %f \t%s\n\n' %(zlo, zhi, 'zlo zhi'))

    mass = [1.00794, 1.]
    file.write('%s\n\n'%'Masses')
    file.write('%d %f\t\n'% (1, mass[0]))
    file.write('%d %f\t\n'% (2, mass[1]))

    file.write('\n%s\n\n'%'Atoms')


    index = 1

    for i in range(data.shape[0]):
        #index type q spin eradius x y z
        theta = random()
        phi = random()
        r = 0.37
        delta_x = r*sin(theta)*cos(phi)
        delta_y = r*sin(theta)*cos(phi)
        delta_z = r*cos(theta)
        file.write('%d %d %f %d %f %f %f %f \n'%(index, 1, 1.0, 0, 0., data[i,0]+delta_x, data[i,1]+delta_y, data[i,2]+delta_z))
        file.write('%d %d %f %d %f %f %f %f \n'%(index+1, 1, 1.0, 0, 0., data[i,0]-delta_x, data[i,1]-delta_y, data[i,2]-delta_z))
        index +=2

    for i in range(data.shape[0]):
        file.write('%d %d %f %d %f %f %f %f \n'%(index, 2, 0, 1, eradius_h, data[i,0], data[i,1], data[i,2]))
        file.write('%d %d %f %d %f %f %f %f \n'%(index+1, 2, 0, -1, eradius_h, data[i,0], data[i,1], data[i,2]))
        index +=2
    file.close()

def write_lammps_eff_h2_li(filename, data, data2, lattice):
    '''write the simplest type of the data file'''
    eradius_h = 0.98
    eradius_li_s = 0.37
    eradius_li_p = 1.75
    #xlo = min(data2[:,0])-1
    xhi = max(data2[:,0]) + 0.5*lattice+1
    xlo =-xhi
    #ylo = min(data2[:,1])-1
    yhi = max(data2[:,1]) + 0.5*lattice+1
    ylo =-yhi
    zlo = min(data2[:,2])
    zhi = max(data2[:,2]) + 0.5*lattice

    filename = filename + '_eff.data'
    file = open(filename, 'w')

    now = datetime.datetime.now()
    file.write('lammpsdata made by wzh in ' + now.strftime('%Y-%m-%d') + '\n')

    file.write('%d\t %s\n' %(data.shape[0]*4+ data2.shape[0]*4,'atoms'))
    atom_type = 3
    file.write('%d\t %s\n\n' %(atom_type, 'atom types'))

    file.write('%f %f \t%s\n' %(xlo, xhi, 'xlo xhi'))
    file.write('%f %f \t%s\n' %(ylo, yhi, 'ylo yhi'))
    file.write('%f %f \t%s\n\n' %(zlo, zhi, 'zlo zhi'))

    mass = [1.00794, 6.941, 1.]
    file.write('%s\n\n'%'Masses')
    file.write('%d %f\t\n'% (1, mass[0]))
    file.write('%d %f\t\n'% (2, mass[1]))
    file.write('%d %f\t\n'% (3, mass[2]))

    file.write('\n%s\n\n'%'Atoms')


    index = 1

    for i in range(data.shape[0]):
        #index type q spin eradius x y z
        theta = random()
        phi = random()
        r = 0.37
        delta_x = r*sin(theta)*cos(phi)
        delta_y = r*sin(theta)*sin(phi)
        delta_z = r*cos(theta)
        file.write('%d %d %f %d %f %f %f %f \n'%(index, 1, 1.0, 0, 0., data[i,0]+delta_x, data[i,1]+delta_y, data[i,2]+delta_z))
        file.write('%d %d %f %d %f %f %f %f \n'%(index+1, 1, 1.0, 0, 0., data[i,0]-delta_x, data[i,1]-delta_y, data[i,2]-delta_z))
        index +=2

    for i in range(data.shape[0]):
        file.write('%d %d %f %d %f %f %f %f \n'%(index, 3, 0, 1, eradius_h, data[i,0], data[i,1], data[i,2]))
        file.write('%d %d %f %d %f %f %f %f \n'%(index+1, 3, 0, -1, eradius_h, data[i,0], data[i,1], data[i,2]))
        index +=2

    for i in range(data2.shape[0]):
        #index type q spin eradius x y z
        file.write('%d %d %f %d %f %f %f %f \n'%(index, 2, 3.0, 0, 0., data2[i,0], data2[i,1], data2[i,2]))
        index +=1

    s = 0
    for i in range(data2.shape[0]):
        file.write('%d %d %f %d %f %f %f %f \n'%(index, 3, 0, 1, eradius_li_s, data2[i,0], data2[i,1], data2[i,2]))
        file.write('%d %d %f %d %f %f %f %f \n'%(index+1, 3, 0, -1, eradius_li_s, data2[i,0], data2[i,1], data2[i,2]))
        file.write('%d %d %f %d %f %f %f %f \n'%(index+2, 3, 0, pow(-1, s), eradius_li_p, data2[i,0], data2[i,1] + 0.5*lattice, data2[i,2]))
        s = s + 1
        index +=3

    file.close()

def write_lammps_eff_he_li_faker(filename, data, data2, lattice):
    '''write the simplest type of the data file'''
    eradius_he = 0.6
    eradius_li_s = 0.37
    eradius_li_p = 1.75
    #xlo = min(data2[:,0])-1
    xhi = max(data2[:,0]) + 0.5*lattice+1
    #ylo = min(data2[:,1])-1
    yhi = max(data2[:,1]) + 0.5*lattice+1
    xlo = -xhi
    ylo = -yhi
    zlo = min(data2[:,2])
    zhi = max(data2[:,2]) + 0.5*lattice

    filename = filename + '_eff.data'
    file = open(filename, 'w')

    now = datetime.datetime.now()
    file.write('lammpsdata made by wzh in ' + now.strftime('%Y-%m-%d') + '\n')

    file.write('%d\t %s\n' %(data.shape[0]*3+ data2.shape[0]*4,'atoms'))
    atom_type = 3
    file.write('%d\t %s\n\n' %(atom_type, 'atom types'))

    file.write('%f %f \t%s\n' %(xlo, xhi, 'xlo xhi'))
    file.write('%f %f \t%s\n' %(ylo, yhi, 'ylo yhi'))
    file.write('%f %f \t%s\n\n' %(zlo, zhi, 'zlo zhi'))

    mass = [4.003, 6.941, 1.]
    file.write('%s\n\n'%'Masses')
    file.write('%d %f\t\n'% (1, mass[0]))
    file.write('%d %f\t\n'% (2, mass[1]))
    file.write('%d %f\t\n'% (3, mass[2]))

    file.write('\n%s\n\n'%'Atoms')


    index = 1

    for i in range(data.shape[0]):
        #index type q spin eradius x y z
        file.write('%d %d %d %f %f %f %f \n'%(index,index, 1, 0, data[i,0], data[i,1], data[i,2]))
        index +=1

    for i in range(data.shape[0]):
        file.write('%d %d %d %f %f %f %f \n'%(index, index, 3, 0 ,data[i,0], data[i,1], data[i,2]))
        file.write('%d %d %d %f %f %f %f \n'%(index+1, index+1, 3, 0, data[i,0], data[i,1], data[i,2]))
        index +=2

    for i in range(data2.shape[0]):
        #index type q spin eradius x y z
        file.write('%d %d %d %f %f %f %f \n'%(index, index, 2, 0., data2[i,0], data2[i,1], data2[i,2]))
        index +=1

    s = 0
    for i in range(data2.shape[0]):
        file.write('%d %d %d %f %f %f %f \n'%(index, index, 3, 0, data2[i,0], data2[i,1], data2[i,2]))
        file.write('%d %d %d %f %f %f %f \n'%(index+1, index+1, 3, 0, data2[i,0], data2[i,1], data2[i,2]))
        file.write('%d %d %d %f %f %f %f \n'%(index+2, index+2, 3, 0, data2[i,0], data2[i,1], data2[i,2]+0.5*lattice))
        s = s + 1
        index +=3



    file.close()


def write_lammps_eff_he_li(filename, data, data2, lattice):
    '''write the simplest type of the data file'''
    eradius_he = 0.6
    eradius_li_s = 0.37
    eradius_li_p = 1.75
    #xlo = min(data2[:,0])-1
    xhi = max(data2[:,0]) + 0.5*lattice+1
    #ylo = min(data2[:,1])-1
    yhi = max(data2[:,1]) + 0.5*lattice+1
    xlo = -xhi
    ylo = -yhi
    zlo = min(data2[:,2])
    zhi = max(data2[:,2]) + 0.5*lattice

    filename = filename + '_eff.data'
    file = open(filename, 'w')

    now = datetime.datetime.now()
    file.write('lammpsdata made by wzh in ' + now.strftime('%Y-%m-%d') + '\n')

    file.write('%d\t %s\n' %(data.shape[0]*3+ data2.shape[0]*4,'atoms'))
    atom_type = 3
    file.write('%d\t %s\n\n' %(atom_type, 'atom types'))

    file.write('%f %f \t%s\n' %(xlo, xhi, 'xlo xhi'))
    file.write('%f %f \t%s\n' %(ylo, yhi, 'ylo yhi'))
    file.write('%f %f \t%s\n\n' %(zlo, zhi, 'zlo zhi'))

    mass = [4.003, 6.941, 1.]
    file.write('%s\n\n'%'Masses')
    file.write('%d %f\t\n'% (1, mass[0]))
    file.write('%d %f\t\n'% (2, mass[1]))
    file.write('%d %f\t\n'% (3, mass[2]))

    file.write('\n%s\n\n'%'Atoms')


    index = 1

    for i in range(data.shape[0]):
        #index type q spin eradius x y z
        file.write('%d %d %f %d %f %f %f %f \n'%(index, 1, 2.0, 0, 0., data[i,0], data[i,1], data[i,2]))
        index +=1

    for i in range(data.shape[0]):
        file.write('%d %d %f %d %f %f %f %f \n'%(index, 3, 0, 1, eradius_he, data[i,0], data[i,1], data[i,2]))
        file.write('%d %d %f %d %f %f %f %f \n'%(index+1, 3, 0, -1, eradius_he, data[i,0], data[i,1], data[i,2]))
        index +=2

    for i in range(data2.shape[0]):
        #index type q spin eradius x y z
        file.write('%d %d %f %d %f %f %f %f \n'%(index, 2, 3.0, 0, 0., data2[i,0], data2[i,1], data2[i,2]))
        index +=1

    s = 0
    for i in range(data2.shape[0]):
        file.write('%d %d %f %d %f %f %f %f \n'%(index, 3, 0, 1, eradius_li_s, data2[i,0], data2[i,1], data2[i,2]))
        file.write('%d %d %f %d %f %f %f %f \n'%(index+1, 3, 0, -1, eradius_li_s, data2[i,0], data2[i,1], data2[i,2]))
        file.write('%d %d %f %d %f %f %f %f \n'%(index+2, 3, 0, pow(-1, s), eradius_li_p, data2[i,0], data2[i,1], data2[i,2]+0.5*lattice))
        s = s + 1
        index +=3



    file.close()




def split(data, k, r):
    '''split based on the origin and the x-y plane'''
    data = data[ (data[:,1] < data[:,0]*k) & (data[:,1] > -data[:,0]*k), :]
    data = data[(data[:,1]**2 + data[:,0]**2 < r**2)]

    return data

def split_circle(data, r):
    data = data[(data[:,1]**2 + data[:,0]**2 < r**2)]

    return data

def conv_split_circle(data, r):
    data = data[(data[:,1]**2 + data[:,0]**2 > r**2)]

    return data

def sin_data(data, r, A, n):
    data = data[(data[:,0]**2+data[:,1]**2)**(1./2) < (r + A*np.sin(n*np.arctan2(data[:,1], data[:,0])))]
    return data

def conv_sin_data(data, r, A, n):
    data = data[(data[:,0]**2+data[:,1]**2)**(1./2) > (r + A*np.sin(n*np.arctan2(data[:,1], data[:,0])) + 1)]
    return data

def lattice(rho, mass, fcc = 4):
    '''density unit ---- g/cm^3, mass ---- g/mol, fcc ---- lattice number'''
    return (fcc*mass*10/6.02/rho)**(1./3)

def rho(lattice, mass, fcc = 4):
    return fcc*mass*10/6.02/lattice**3;

def lattice_center(rho, mass, num=1):
    return (num*mass*10/6.02/rho)**(1./3)


def angel_to_k(angel):
    return tan(angel/2./180*pi)







r = 100
rs = int(r*0.6)
angel = 15
rho = 0.12
lat = lattice(rho, 4.003, 4)

r = int(600/lat)
r = 1000
print r*lat
rs = int(r*0.6)
data = fcc(0,r,-r,r,-4,4, lat)
k = angel_to_k(angel)
data = split(data, k, r*lat)
filename = 'den%.2f_r%d' % (rho, r)
write_lammps(filename, data, lat, k, r)
print k,lat, r*lat, data.shape


for i in range(5,11):
    i = i*0.5
    lat = lattice(i, 4.003, 4)
    data = fcc(0, 5, 0, 5, 0, 5, lat)
    filename = 'rho%0.1f' % i
    write_lammps(filename, data, lat, k=1, r=5)

angel = 15
rho = 0.12
lat = lattice(rho, 4.003, 4)
k = angel_to_k(angel)
for r in range(1,4):
    r = r*100
    data = fcc(0,r, -r, r, -4, 4, lat)
    data = split(data, k , r*lat)
    filename = 'ang30_r%d' % ( r)
    write_lammps(filename, data, lat, k, r)


rho = 0.12
lat = lattice(rho, 4.003, 4)
r = 15
data = fcc(-r, r, -r, r, -2, 2, lat)
data = split_circle(data, r*lat)
filename = 'small20%d' %(r)
write_lammps_eff_he(filename, data, lat, k, r)


rho = 0.12
lat = lattice(rho, 1.00794*2, 4)
r = 20
data = fcc(-r, r, -r, r, -2, 2, lat)
data = split_circle(data, r*lat)
filename = 'h2_20%d' %(r)
write_lammps_eff_h2(filename, data, lat, k, r)

rho = 0.12
lat1 = 4.42
r = 60
data = center(-r, r, -r, r, -2, 2, lat1)
data = split_circle(data, r*lat1)
lat = 4.42
r = 80
rs = 60
data2 = fcc(-r, r, -r, r, -2 , 2, lat)
data2 = split_circle(data2, r*lat)
data2 = conv_split_circle(data2, rs*lat+3)

all = np.concatenate([data, data2])

one1 = np.zeros(data.shape[0])+2
one2 = np.zeros(data2.shape[0])+1

test1 = np.c_[data, one1]
test2 = np.c_[data2, one2]

test = np.concatenate([test1, test2])

filename = '1017_Li_H2_eff'
write_lammps_common(filename,test, lat1)
write_lammps_eff_he_li(filename, data, data2, lat)

# Cu-He
Cu_lat = 3.615
R = 100+1
r = 70
n = 10
k = 0.05
A = k/n*2*pi*r
#A=0
Cu = fcc(-R, R, -R, R, -2, 2, Cu_lat)
R = 100
Cu = split_circle(Cu, R*Cu_lat)
Cu = conv_sin_data(Cu, r*Cu_lat, A*Cu_lat, 8)

He_lat =5
He = fcc(-R, R, -R, R, -2, 2, Cu_lat)
He = sin_data(He, r*Cu_lat, A*Cu_lat, 8)

one1 = np.zeros(Cu.shape[0])+2
one2 = np.zeros(He.shape[0])+1

test1 = np.c_[Cu, one1]
test2 = np.c_[He, one2]

test = np.concatenate([test1, test2])

filename = 'Cu_He_classic_R100_r70_mol8'
write_lammps_common(filename,test, Cu_lat)

all = np.concatenate([Cu, He])
filename = 'Cu_He'
write_lammps(filename, all, Cu_lat)


filename = 'Cu100'
write_lammps(filename, Cu, Cu_lat)
