import gsd
import gsd.hoomd
import numpy as np
import freud
import os
import glob

names = [os.path.basename(x)[:] for x in sorted(glob.glob('b*.gsd'))]
# the input file name
input_fname = names[0]
# the output file name 
output_fname = 'unwrap.gsd'

config = gsd.hoomd.open(input_fname, mode='rb').read_frame(0)
wrap_position = config.particles.position
image = config.particles.image
box = config.configuration.box

freud_box = freud.box.Box.from_box(box)

unwrap_position = freud_box.unwrap(wrap_position, image)
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
snap.particles.position = unwrap_position
snap.particles.typeid = config.particles.typeid
snap.particles.mass = [1.0]*config.particles.N 
snap.particles.charge = [0.0]*config.particles.N 

# give the particle a type name to refer to later
snap.particles.types = config.particles.types  

# add bonds between the monomers
snap.bonds.types = config.bonds.types
snap.bonds.N = config.bonds.N
snap.bonds.group = config.bonds.group
snap.bonds.typeid = config.bonds.typeid

snap.configuration.box = box


with gsd.hoomd.open(name=output_fname, mode='wb') as f:
    f.append(snap)
