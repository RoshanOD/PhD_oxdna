import argparse
import os
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-N","--iteration",help="index of iteration")
parser.parse_args()
args = parser.parse_args()

# print(f"{args.directory}")
# print(type(args.iteration))
# print('a string args.iteration')
#----Reading reaction coordinate (num_bonds, distance bin, weights)
data_set = []
react_coord, weights = [],[]
with open("/rds/general/user/rjo20/home/oxdna_work/MC/18mer_trap/RUN"+str(int(args.iteration)-1)+"/energy.dat") as f:
    print(f)
    for line in f:
        useful = line.split(' ')[-4:]
        react_coord.append(float(useful[0])-float(useful[2])); weights.append(float(useful[-1]))
data_set.append([react_coord,weights])
# print(np.shape(data_set))

min_val, max_val = -3.0, 18.0 # of reaction coordinate
freq,bins = np.histogram(data_set[0], bins=np.arange(min_val-0.5,max_val+1.5,1))
avg_freq = len(data_set[0][0])/(max_val-min_val+1)
print("Freq per bin"+str(int(avg_freq)))
# Careful; if run was cut short/didn't end, len(data_set[i][0]) will be too short
resid = freq-avg_freq
norm_res = resid/np.std(resid)

#----Updating new weights for next round
old_w = []
with open("/rds/general/user/rjo20/home/oxdna_work/MC/18mer_trap/RUN" + str(int(args.iteration)-1) + "/wfile.txt") as wf:
    print(wf)
    for line in wf:
        old_w.append(float(line.split(' ')[-1]))
new_w = np.round(old_w*np.exp(-norm_res),6)
new_w = [w if w!=0. else 1e-6 for w in new_w]

out_file = "/rds/general/user/rjo20/home/oxdna_work/MC/18mer_trap/RUN" + args.iteration + "/wfile.txt"

col1 = [0,0,0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
col2 = [3,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
new_wfile = open(out_file, 'w')
for i in range(len(col1)):
    write_str = "%d %d %f\n"%(col1[i],col2[i],new_w[i]) # string formatting
    new_wfile.write(write_str)
new_wfile.close()
