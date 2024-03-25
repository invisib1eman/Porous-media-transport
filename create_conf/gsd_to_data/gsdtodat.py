import sys
import gsd
import gsd.hoomd
import numpy as np


class Bond(object):
    def __init__(self, btype, i, j):
        self.btype = btype
        self.i = i
        self.j = j

class SystemEnity(object):
    """Defines a general ligand
    """

    def __init__(self, n_particles, n_bonds=0, n_angles=0, n_dihedrals=0, n_impropers=0):
        """

        :param n_particle: the number of particles in a System
        """

        # number of total particles
        self.par_num = n_particles
        self.par_types = 0
        self.par_position = np.zeros((self.par_num, 3))
        self.par_image = np.zeros((self.par_num, 3), dtype=int)
        self.par_typeid = np.zeros(self.par_num, dtype=int)
        self.par_charge = np.zeros(self.par_num)
        self.par_velocity = np.zeros((self.par_num, 3))
        # number of bonds
        self.bond_num = n_bonds
        self.bond_types = 0
        self.bond_typeid = np.zeros(self.bond_num, dtype=int)
        self.bond_group = np.zeros((self.bond_num, 2), dtype=int)
        # number of angles
        self.angle_num = n_angles
        self.bond_types = 0
        self.angle_typeid = np.zeros(self.angle_num, dtype=int)
        self.angle_group = np.zeros((self.angle_num, 3), dtype=int)
        self.box = 0


# the input file name
input_fname = sys.argv[1]
# the output file name
output_fname = sys.argv[2]
# nchain, chain_length
length1 = 40 + 6
nchain1 = 870
length2 = 1 + 2
nchain2 = 3*100*8

config = gsd.hoomd.open(input_fname, mode='rb').read_frame(0)
wrap_position = config.particles.position
image = config.particles.image
box = config.configuration.box


unwrap_position = np.zeros_like(wrap_position)
temp_image = np.zeros_like(image)
for tag in range(config.particles.N):
    unwrap_position[tag] = config.particles.position[tag] + \
        np.multiply(config.particles.image[tag], config.configuration.box[:-3])
    temp_image[tag] = np.array([0, 0, 0])

Sys = SystemEnity(config.particles.N,  config.bonds.N)
temp_box = box[:3]/2.0
Sys.par_position = unwrap_position + temp_box
Sys.par_typeid = config.particles.typeid + 1
Sys.par_charge = [0.0]*config.particles.N
Sys.par_velocity = config.particles.velocity
Sys.par_image = temp_image
Sys.par_types = len(config.particles.types)
# add bonds between the monomers
Sys.bond_types = len(config.bonds.types)
Sys.bond_group = config.bonds.group + 1
Sys.bond_typeid = config.bonds.typeid + 1

# Sys.configuration.box = box*2*np.amax(np.abs(image))
Sys.box = box

datafile = open(output_fname, 'w')
mol_id = open("mol_id.txt", 'w')
mol_id.write("# id mol type\n")

datafile.write("LAMMPS data file\n")
datafile.write("%d atoms\n" % Sys.par_num)
datafile.write("%d bonds\n\n" % Sys.bond_num)
# datafile.write("%d angles\n\n" % Sys.angleself.angle_num)

datafile.write("%d atom types\n" % Sys.par_types)
#datafile.write("%d bond types\n\n" % Sys.bond_types)
datafile.write("%d bond types\n\n" % 1)

datafile.write(str(0) + " " + str(Sys.box[0]) + " xlo xhi\n")
datafile.write(str(0) + " " + str(Sys.box[1]) + " ylo yhi\n")
datafile.write(str(0) + " " + str(Sys.box[2]) + " zlo zhi\n\n")

datafile.write("Atoms # full\n\n")

for i in range(Sys.par_num):
    s = str(i+1) + " "
    if i < (length1*nchain1):
        mol = int(i/length1) + 1
    else:
        mol = int((i-length1*nchain1)/length2) + nchain1 + 1
    ss = str(i+1) + " " + str(mol) + " " + str(Sys.par_typeid[i]) + '\n'
    mol_id.write(ss)
    mol = 1
    s += str(mol) + " "
    s += str(Sys.par_typeid[i]) + " "
    s += str(Sys.par_charge[i]) + " "
    s += str(Sys.par_position[i][0]) + " " + str(Sys.par_position[i]
                                                 [1]) + " " + str(Sys.par_position[i][2]) + " "
    s += str(Sys.par_image[i][0]) + " " + str(Sys.par_image[i]
                                              [1]) + " " + str(Sys.par_image[i][2]) + "\n"
    datafile.write(s)


datafile.write("\nVelocities\n\n")
for i in range(Sys.par_num):
    s = str(i+1) + " "
    s += str(Sys.par_velocity[i][0]) + " " + str(Sys.par_velocity[i]
                                                 [1]) + " " + str(Sys.par_velocity[i][2]) + "\n"
    datafile.write(s)

datafile.write("\nBonds\n\n")
for i in range(Sys.bond_num):
    #s = str(i+1) + " " + str(Sys.bond_typeid[i]) + " "
    s = str(i+1) + " " + str(1) + " "
    s += str(Sys.bond_group[i][0]) + " " + str(Sys.bond_group[i][1]) + "\n"
    datafile.write(s)

'''
datafile.write("\nAngles\n\n")

for i in range(Sys.angle_num):
    s = str(i+1) + " " + str(Sys.angle_typeid[i]) + " "
    s += str(Sys.angle_group[i][0]) + " " + str(Sys.angle_group[i]
                                              [1]) + " " + str(Sys.angle_group[i][2]) + "\n"
	datafile.write(s)
'''

datafile.close()
mol_id.close()
