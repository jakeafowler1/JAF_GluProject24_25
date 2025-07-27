#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# In[ ]:


hits = "/nobackup/cm16jf/databases/protein_faa_reps/cyanobacteria_seqs/tabloutcyano.txt"
searchForHits = "/nobackup/cm16jf/databases/protein_faa_reps/cyanobacteria_seqs/nrcyano_seqs.faa"
writefile = "/nobackup/cm16jf/databases/protein_faa_reps/cyanobacteria_seqs/cyanoHitsseq.faa"

f = open(hits, "r") # Open the file in read mode that contains the hits
foundHits = [] # Predefine list
ind = 0;

# Loop through each line in the file, the contents of each line will output to 'l'
for l in f:
    currLine = l
    endString = currLine.find("-") #Finds the first occurance of '-' which defines the protien code
    isolatedString = currLine[0:endString].strip() #Isolates the string and removes the spaces
    foundHits.append(isolatedString) # Populate the list with each line in the file
    
# Breaks the for loop when a '#' is read indicating the end of the file, but also ignores the first 3 lines (the title lines)

    if (isolatedString == "#") and (ind > 3):
        break

    
    ind = ind + 1
    
    

del foundHits[0:3] # Removes the first 3 lines to just give the protiens codes
del foundHits[-1] # Removes the hash at the end

# Opening the file to find the hits in
f2 = open(searchForHits, "r")

isolatedHits = []
foundLine = False
ind = 0

# loop through each line in the file one at a time
for l2 in f2:
    if l2[0] == ">":
        if foundLine == True:
            foundLine = False

        
        protCode = l2[1:-1]
        if protCode in foundHits:
            foundLine = True

    if foundLine == True:
        isolatedHits.append(l2)

f2.close()

f3 = open(writefile, "a")
f3.writelines(isolatedHits)
f3.close()
