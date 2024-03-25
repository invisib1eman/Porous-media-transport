import sys
import numpy as np

data = None
# the input file name
if len(sys.argv) > 3:
    input_fname = sys.argv[3]
    data=(np.genfromtxt(fname=input_fname, skip_header=0).T).astype(int)

# nchain, chain_length
nchain1 = 800
nchain2 = 1600

length1 = 40 + 6
length2 = 1 + 2
num_AC = nchain1*length1
num_BC = nchain2*length2
par_num = int(sys.argv[1])

mol_id = open(sys.argv[2], 'w')
mol_id.write("# id mol\n")

id_par = 0
if data is None:
    for i in range(par_num):
        if i < num_AC:
            mol = int(i/length1) + 1
        else:
            mol = int((i-num_AC)/length2) + nchain1 + 1
        ss = str(i+1) + " " + str(mol) + '\n'
        mol_id.write(ss)
else:
    for i in range(len(data[1])):
        mol = i+1
        for j in range(data[1][i]):
            id_par += 1
            #if (j+1) > 40:
            #    type_par = 1
            #else:
            #    type_par = 3
            #ss = str(id_par) + " " + str(mol) + " " + str(type_par) + '\n'
            ss = str(id_par) + " " + str(mol) +  '\n'
            mol_id.write(ss)
    
    num_AC = np.sum(data[1])
    nchain1 = len(data[1])
    for i in range(par_num-num_AC):
        id_par = num_AC + i + 1
        #if i%3 == 1:
        #    type_par = 3
        #else:
        #    type_par = 2
        mol = int(i/length2) + nchain1 + 1
        #ss = str(id_par) + " " + str(mol) + " " + str(type_par) + '\n'
        ss = str(id_par) + " " + str(mol) + '\n'
        mol_id.write(ss)


mol_id.close()
