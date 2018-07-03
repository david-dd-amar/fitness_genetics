# This script receives a bfile path
# It first transforms the bed to ped and keeps it in a temp dir
# It then goes over all the snps and recodes indels:
# For each SNP that has at least one allele size 2 or more (I):
# smaller allele is recoded as D and the longer as I.
# Map file stays the same.
# We then recode back to bed and delete the temp folder.
import re
import random
import os

def read_bim(p):
    f = open(p)
    d = {}
    for l in f:
        l = l.rstrip()
        arr = re.split("\\s",l)
        #if arr[4]=="0" or arr[5]=="0":print(l)
        rid = arr[0]+';'+arr[3]
        d[rid] = arr
    f.close()
    return d
def read_bim_ids_as_array(p):
    f = open(p)
    ids = []
    for l in f:
        l = l.rstrip()
        arr = re.split("\\s",l)
        rid = arr[0]+';'+arr[3]
        ids.append(rid)
    f.close()
    return ids

bfile = "/Users/David/Desktop/repos/fitness_genetics/python_sh/merged_control_geno"
bim = bfile +".bim"
out_path = "/Users/David/Desktop/repos/fitness_genetics/python_sh/merged_control_geno_recoded"
sep = " "
plink_cmd = "/Users/David/Desktop/elite/analysis/plink_mac/plink"
tmp = os.path.dirname(os.path.realpath(__file__)) + "/tmp"+str(random.randint(1, 1000))
os.system("mkdir " + tmp)

print("reading the bim, searching for the indels")
d = read_bim(bim)
d1 = {}
indel_keys = {}
for k in d:
    arr = d[k]
    s1 = arr[4]
    s2 = arr[5]
    if len(s1)>1:
        newarr = [arr[0],arr[1],arr[2],arr[3],"I","D"]
        d1[k]=newarr
        indel_keys[k]=1
        continue
    if len(s2)>1:
        newarr = [arr[0],arr[1],arr[2],arr[3],"D","I"]
        d1[k]=newarr
        indel_keys[k]=1
        continue
    d1[k]=arr

bim_ids = read_bim_ids_as_array(bim)
indel_key_inds = {}
i=0
for k in bim_ids:
    if k in indel_keys:
        indel_key_inds[i]=i
    i+=1
print ("found " + str(len(indel_key_inds)) + " indels in the bim file")

print("creating the ped file")
# Create the ped
import subprocess
command = plink_cmd +" --bfile " + bfile + " --recode --out " + tmp + "/ped"
print("the command is:")
print(command)
process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
process.wait()

print("recoding the ped")
# Read the ped and recode
ped_path = tmp+"/ped.ped"
new_ped_path = tmp+"/new.ped"
f = open(ped_path)
o = open(new_ped_path,"w")
for l in f:
    l = l.rstrip()
    arr = re.split("\\s",l)
    for j in range(6,len(arr),2):
        j_ind = (j-6)/2
        if j_ind in indel_key_inds:
            s1 = arr[j]
            s2 = arr[j+1]
            if len(s1)>1 and s1 != "0" and s1!="NA":
                arr[j] = "I"
            if len(s2)>1 and s2 != "0" and s2!="NA":
                arr[j+1] = "I"
            if len(s1)<=1 and s1 != "0" and s1!="NA":
                arr[j] = "D"
            if len(s2)<=1 and s2 != "0" and s2!="NA":
                arr[j+1] = "D"
    newl = " ".join(arr)
    o.write(newl+"\n")
o.close()

print("transforming the new ped to bed")
map_path = tmp+"/ped.map"
command = plink_cmd +" --ped " + new_ped_path +" --map " + map_path  +  " --make-bed --out " + out_path
print("the command is:")
print(command)
process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
process.wait()

os.system("rm -r " + tmp)






