import warnings
warnings.filterwarnings('ignore')

import numpy as np
from numpy import linalg as la
import mbuild as mb
import sys
import random

#import freud
#from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield
#from mbuild.formats.gsdwriter import write_gsd
#import foyer
#import unyt as u
#import gsd

class Quaternion(object):
    """
        quaternion class
    """
    def __init__(self, orientation, hoomd=False):
        """
        :param hoomd_orientation: quaternion in hoomd format [real,x,y,z] 
        """

        if (hoomd):
            self.q = (np.array((orientation[1], orientation[2], orientation[3], orientation[0]), dtype=np.float64))
        else:
            self.q = (np.array((orientation[0], orientation[1], orientation[2], orientation[3]), dtype=np.float64))

    def multiply(self, quaternion2):
        """
        
        :param quaternion2: quaternion to multipply this quaternion by 
        :return: another quaternion which is the result of the multiplication
        """
        x0, y0, z0, w0 = quaternion2.q
        x1, y1, z1, w1 = self.q
        new = np.array((
            x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
            -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
            x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0,
            -x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0), dtype=np.float64)
        return Quaternion(new)

    def conjugate(self):
        """
        
        :return: the quaternion's conjugate as a quaternion object 
        """
        return Quaternion(np.array((-self.q[0], -self.q[1], -self.q[2], self.q[3]), dtype=np.float64))

    def orient(self, vector):
        """
        :param vector: the vector to be rotated [x, y, z] 
        :return: the rotated vector [x, y, z]
        """
        q2 = [vector[0], vector[1], vector[2], 0.0]
        v_quat = Quaternion(q2)
        return self.multiply(v_quat).multiply(self.conjugate()).q[:3]

    def inverse(self):
        """
        :return: the inverse quaternion 
        """
        q0, q1, q2, q3 = self.q
        bottom = q0 ** 2 + q1 ** 2 + q2 ** 2 + q3 ** 2
        q4 = np.divide([-q0, -q1, -q2, q3], bottom)
        return Quaternion(q4)

    def de_orient(self, position):
        """
        :param position: position vector to be deoriented [x, y, z]
        :return: the "unrotated" position
        """

        return self.inverse().orient(position)

class QuaternionBetween(Quaternion):
    """
    calculates Quaternion between 2 Vectors
    """

    def __init__(self, vector1, vector2, hoomd=False):
        """
        
        :param vector1: vector the quaternion goes from
        :param vector2: vector the quaternion goes to
        :param hoomd: set to true if you want hoomd style Quaternion
        """

        cross = np.cross(vector1, vector2)
        w = la.norm(vector1) * la.norm(vector2) + np.dot(vector1, vector2)
        length = la.norm([cross[0], cross[1], cross[2], w])
        unit = [cross[0], cross[1], cross[2], w] / length

        super(QuaternionBetween, self).__init__(unit, hoomd=hoomd)

def points_on_unit_sphere(n):
        """calculates n points distributed on unit sphere

        :param n: number of points to distribute
        :return: positions of points on unit sphere
        """

        pts = []
        inc = np.pi * (3 - np.sqrt(5))
        off = 2 / float(n)
        for k in range(n):
            y = k * off - 1 + (off/2)
            r = np.sqrt(1 - y*y)
            phi = k * inc
            pts.append([np.cos(phi)*r, y, np.sin(phi)*r])
        return np.array(pts)

class Bead(mb.Compound):
    def __init__(self):
        super(Bead, self).__init__()
        # Add bead
        self.add(mb.Particle(name='C', pos=[0,0,0]), label='m[$]')
        self.add(mb.Port(anchor=self[0], orientation=[1, 0, 0], separation=bond_length), label='right')
        self.add(mb.Port(anchor=self[0], orientation=[-1, 0, 0], separation=bond_length), label='left')

class Center_bead(mb.Compound):
    def __init__(self, num=4):
        super(Center_bead, self).__init__()
        # Add bead
        self.add(mb.Particle(name='C'), label='c[$]')
        #Add ports and rotate the direction of the ports
        direction=points_on_unit_sphere(num)
        for i in range(num):
            Label = 'label' + str(i)
            self.add(mb.Port(anchor=self[0], orientation=list(direction[i]), separation=bond_length), label=Label)

class End_bead(mb.Compound):
    def __init__(self):
        super(End_bead,self).__init__()
        # Add the end bead
        self.add(mb.Particle(name='A'), label='a[$]')
        # Add ports
        self.add(mb.Port(anchor=self[0], orientation=[1, 0, 0], separation=bond_length), label='right')

class Linear_chain(mb.Compound):
    def __init__(self, chain_length=1, cap_front=True, cap_end=False):
        super(Linear_chain, self).__init__()
        end_bead = End_bead()
        last_monomer = Bead()
        if cap_front:
            mb.force_overlap(move_this=end_bead,
                            from_positions=end_bead['right'],
                            to_positions=last_monomer['right'])
            self.add(end_bead, label='right-cap')
        self.add(last_monomer, label='ma[$]')
        for i in range (chain_length-1):
            current_monomer = Bead()
            mb.force_overlap(move_this=current_monomer,
                             from_positions=current_monomer['right'],
                             to_positions=last_monomer['left'])
            self.add(current_monomer)
            last_monomer=current_monomer
        if cap_end:
            end_bead = End_bead()
            mb.force_overlap(move_this=end_bead,
                            from_positions=end_bead['right'],
                            to_positions=last_monomer['left'])
            self.add(end_bead, label='left-cap')

class Cross_linker(mb.Compound):
    def __init__(self, chain_length=1):
        super(Cross_linker, self).__init__()
        end_bead1 = End_bead()
        end_bead1[0].name = 'B'
        self.add(end_bead1, label='right-cap')
        if chain_length > 0:
            last_monomer = Bead()
            mb.force_overlap(move_this=end_bead1,
                            from_positions=end_bead1['right'],
                            to_positions=last_monomer['right'])
            self.add(last_monomer, label='ma[$]')
            for i in range (chain_length-1):
                current_monomer = Bead()
                mb.force_overlap(move_this=current_monomer,
                                from_positions=current_monomer['right'],
                                to_positions=last_monomer['left'])
                self.add(current_monomer)
                last_monomer=current_monomer
            end_bead2 = End_bead()
            end_bead2[0].name = 'B'
            mb.force_overlap(move_this=end_bead2,
                            from_positions=end_bead2['right'],
                            to_positions=last_monomer['left'])
            self.add(end_bead2, label='left-cap')
        else:
            end_bead2 = End_bead()
            end_bead2[0].name = 'B'
            mb.force_overlap(move_this=end_bead1,
                            from_positions=end_bead1['right'],
                            to_positions=end_bead2['right'])
            self.add(end_bead2, label='left-cap')



class Star_chain(mb.Compound):
    def __init__(self, arm_num=1, arm_len=3):
        super(Star_chain, self).__init__()
        # Determine the number of port according to the value of arm_num (arm_num:the number of arm)
        direction = points_on_unit_sphere(arm_num)
        cent = Center_bead(arm_num)
        # Build the linear arm accroding to the arm length
        self.add(cent)
        for i in range(arm_num):
            linear = Linear_chain(arm_len)
            vec = linear[-1].pos - linear[0].pos
            q = QuaternionBetween(vec, direction[i])
            for x in range(arm_len+1):
                linear[x].pos = q.orient(linear[x].pos)
            linear.add(mb.Port(
                anchor=linear[-1], orientation=direction[i], separation=bond_length), label='left')
            mb.force_overlap(move_this=linear,
                             from_positions=linear['left'],
                             to_positions=cent['label'+str(i)])
            self.add(linear)


class Branch_chain(mb.Compound):
    def __init__(self, back_len=6, branch_num=3, branch_len=2, seed=1245, disorder=0):
        super(Branch_chain, self).__init__()

        # Build the backbone accroding to the backbone length
        backbone = Linear_chain(back_len, cap_front=False)
        self.add(backbone)
        # Determine the number of port according to the value of arm_num (arm_num:the number of arm)
        
        '''
        if disorder == 0:
            np.random.seed(seed)
        rand_num = []
        while len(rand_num) < branch_num:
            val_rand = np.random.randint(0, back_len)
            if val_rand not in rand_num:
                rand_num.append(val_rand)
        rand_num.sort()
        '''
        rand_array = []
        for i in range(back_len):
            if i < branch_num:
                rand_array.append(1)
            else:
                rand_array.append(3)
        for j in range(1000):
            random.shuffle(rand_array)
        rand_array = np.array(rand_array)

        rand_num = list(np.nonzero(rand_array == 1)[0])
        rand_num.sort()
        
        for i in range(branch_num):
            Label = 'label' + str(i)
            backbone.add(mb.Port(anchor=backbone[int(rand_num[i])], orientation=[0, 1, 0], separation=bond_length), label=Label)

        for i in range(branch_num):
            if branch_len > 1:
                linear = Linear_chain(branch_len-1)
                linear.add(mb.Port(anchor=linear[-1], orientation=[-1, 0, 0], separation=bond_length), label='right')
                linear.rotate(-np.pi/2, around=[0, 1, 0])
            else:
                linear = End_bead()
                linear['right'].rotate(np.pi/2, around=[0, 1, 0])
            mb.force_overlap(move_this=linear,
                            from_positions=linear['right'],
                            to_positions=backbone['label'+str(i)])
            self.add(linear)



class Copolymer(mb.Compound):
    def __init__(self, num_C=6, num_A=3,  seed=1245, disorder=0):
        super(Copolymer, self).__init__()

        # Determine the number of port according to the value of arm_num (arm_num:the number of arm)
        rand_num = []
        if disorder == 0:
            np.random.seed(seed)
        while len(rand_num) < num_A:
            val_rand = np.random.randint(0, num_A+num_C)
            if val_rand not in rand_num:
                rand_num.append(val_rand)
        
        rand_num.sort()
        #print(rand_num)
        
        if  0 in rand_num:
            temp_bead = Bead()
            temp_bead[0].name = 'A'
            last_monomer= temp_bead
            self.add(last_monomer)
        else:
            last_monomer= Bead()
            self.add(last_monomer)
        for i in range(1, num_A+num_C):
            if i in rand_num:
                temp_bead = Bead()
                temp_bead[0].name = 'A'
                current_monomer = temp_bead
            else:
                current_monomer = Bead()
            mb.force_overlap(move_this=current_monomer,
                             from_positions=current_monomer['right'],
                             to_positions=last_monomer['left'])
            self.add(current_monomer)
            last_monomer=current_monomer


bond_length = 0.05
num_poly = 782 # the number of the chain
num_clinker = 1876 # the number of the cross-linker
A_perchanin = 6
mainchain_len = 40
################################################################
rng = np.random.default_rng()
NAvaluesi = np.zeros(num_poly, dtype=int)
for ii in range(10000):
    NAvaluesi = rng.binomial(mainchain_len, A_perchanin/mainchain_len, num_poly)
    if np.sum(NAvaluesi) == A_perchanin*num_poly:
        break
    else:
        continue

print(np.sum(NAvaluesi))

s = 0
s2 = 0
for i in range(0, num_poly):
  s = s + NAvaluesi[i]
  s2 = s2 + NAvaluesi[i]*NAvaluesi[i]
s = s / float(num_poly)
stdevABS = np.sqrt(s2/float(num_poly) - s**2)
print("Absolute standard deviation = %g" % (stdevABS))

fmol = open( "MolId.dat", 'w')
####################################################################

cubic_box_length = 40
box = mb.Box([cubic_box_length, cubic_box_length, cubic_box_length])
#freud_box = freud.box.Box.from_box(box.lengths)

disorder = 1
# 1) disorder = 1 which indicates each chain has different random squences for A beads
# 2) disorder = 0 which indicates each chain has the same random squences for A beads
coply = []
if disorder == 1:
    ncomp=[1]*num_poly + [num_clinker]
    for i in range(num_poly):
        # 17 = the number of C bead, 6 = the number of A bead
        #poly=Copolymer(17, 6, 123, disorder=0)
        # 40 = the length of the backbone, 6 = the number of branch point, 1 = the length of the branch, 
        poly = Branch_chain(40, NAvaluesi[i], 1, 23513, disorder=disorder)
        total_chainLen = mainchain_len + NAvaluesi[i]
        fmol.write(str(i+1) + '\t' + str(total_chainLen) + '\n')
        coply.append(poly)

elif disorder == 0:
    ncomp=[num_poly] + [num_clinker]
    # 17 = the number of C bead, 6 = the number of A bead
    #poly=Copolymer(17, 6, 123, disorder=0)
    # 40 = the length of the backbone, 6 = the number of branch point, 1 = the length of the branch, 
    poly = Branch_chain(40, 6, 1, 23513, disorder=disorder)
    coply.append(poly)

#for star polymer
#4=the number of arm, 3+1=the length of arm
#poly = Star_chain(4, 3)
#coply.append(poly)
# for copolymer 
#cross_linker = Cross_linker(0)
# for branch chain
cross_linker = Cross_linker(1)
coply.append(cross_linker)
filled_box = mb.fill_box(compound=coply,
                      n_compounds=ncomp,
                      box=box, overlap=0.5,
                      )

#filled_box.xyz = freud_box.wrap(filled_box.xyz)

# Apply foyer force field
#ff = foyer.Forcefield('ff.xml')
#structure = ff.apply(filled_box)

#filled_box.visualize()
fname = sys.argv[1] + '.gsd'
filled_box.save(filename=fname, overwrite=True)
