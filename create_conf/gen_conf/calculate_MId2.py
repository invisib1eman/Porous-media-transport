import sys
import numpy as np
import warnings

data = None
bond_data = None
# the input file name
if len(sys.argv) > 1:
    input_fname = sys.argv[1]
    data=(np.genfromtxt(fname=input_fname, skip_header=0).T).astype(int)
    bond_data=(np.genfromtxt(fname=sys.argv[2], skip_header=9).T).astype(int)
    bond_data=bond_data[0:3]
    arg = np.argsort(bond_data[1])
    bond_data[1] = bond_data[1][arg]
    bond_data[2] = bond_data[2][arg]

# nchain, chain_length
length1 = 40 + 6
nchain1 = 800
length2 = 1 + 2
nchain2 = 2*800
num_AC = 36800
num_BC = 4800
par_num = num_AC + num_BC

mol_id = open("mol_id_type.txt", 'w')
mol_id.write("# id mol type\n")

id_par = 0
cid = 0
sum_val = 0
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
            if (j+1) > 40:
                type_par = 1
                cid = bond_data[2][sum_val + j-1]
            else:
                type_par = 3
                cid = 0
            ss = str(id_par) + "\t" + str(mol) + "\t" + str(type_par)  + "\t" + str(cid) + '\n'
            #ss = str(id_par) + " " + str(mol) +  '\n'
            mol_id.write(ss)
        sum_val += data[1][i] - 1

    cid = 0
    num_AC = np.sum(data[1])
    nchain1 = len(data[1])
    
    """
    for i in range(par_num-num_AC):
        id_par = num_AC + i + 1
        if i%3 == 1:
            type_par = 4
        else:
            type_par = 2
        mol = int(i/length2) + nchain1 + 1
        ss = str(id_par) + " " + str(mol) + " " + str(type_par)  + '\n'
        #ss = str(id_par) + " " + str(mol) + '\n'
        mol_id.write(ss)
    """

mol_id.close()
