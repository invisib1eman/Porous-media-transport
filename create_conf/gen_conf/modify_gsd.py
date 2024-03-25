import gsd
import gsd.hoomd
import numpy as np
#import freud
import os
import glob
import sys

# the input file name
input_fname = sys.argv[1]
# the output file name 
output_fname = sys.argv[2]

config = gsd.hoomd.open(input_fname + '.gsd', mode='rb').read_frame(0)
wrap_position = config.particles.position
image = config.particles.image
box = config.configuration.box

#freud_box = freud.box.Box.from_box(box)
#unwrap_position = freud_box.unwrap(wrap_position, image)
"""
unwrap_position = []
for tag in range(config.particles.N):
    temp = config.particles.position[tag] + np.multiply(config.particles.image[tag], config.configuration.box[:-3])
    unwrap_position.append(temp)
"""
# create an unwrap position configuration of the chain 
snap = gsd.hoomd.Snapshot()

# add the monomers to the configuration
snap.particles.N = config.particles.N 
snap.particles.image = image
snap.particles.position = wrap_position
snap.particles.typeid = config.particles.typeid
snap.particles.mass = [1.0]*config.particles.N 
snap.particles.charge = [0.0]*config.particles.N 
snap.particles.types = config.particles.types  

# add bonds between the monomers
snap.bonds.types = config.bonds.types
snap.bonds.N = config.bonds.N
snap.bonds.group = config.bonds.group
snap.bonds.typeid = config.bonds.typeid

#snap.configuration.box = box*2*np.amax(np.abs(image))
snap.configuration.box = box

with gsd.hoomd.open(name=output_fname + '.gsd', mode='wb') as f:
    f.append(snap)
