import os, shutil

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

path = "/nobackup/cm16jf/databases/protein_faa_reps/bacteriaseqs/protein_faa_reps"
namesFile = "/nobackup/cm16jf/databases/protein_faa_reps/bacteriaseqs/protein_faa_reps/bacteria_gc90_codes.txt"

#Reads each line and adds it to a list

with open(namesFile, "r") as f:
    data = f.read().splitlines()

#Trys to make the directory if it does not already exist

try:

    os.mkdir(path + "/bacteria_gc90_seqs")
except OSError:
    print ("Creation of the directory %s failed" % path)
else:
    print ("Successfully created the directory %s " % path)

    
for fileName in data:
    print("Looking for: " + fileName + "_protein.faa ...")
    x = find(fileName + "_protein.faa", path + "/bacteria")
    if x != None:
        print(x + " found!")
        shutil.copy(x, path + "/bacteria_gc90_seqs" )
    else:
        print("Not found")

