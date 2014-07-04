#!/usr/bin/env python

""" Print the sequence of temperatures that a given structure
visits in the replica exchange simulation. """

import sys
import commands

istart=0
iold=istart-1

# must match T values in config file ...
temp=[293,297,301,305,309,313,317,321,325,329,334,339,344,349,354,359,364,369,374,379,385,391,397,403,409,415,421,427,433,439]

t0=temp[0]

# A dictionary of the structures at the various temperatures in the current step
curr={}
currtemp={}
# Initialize the current structures
for t in temp:
	curr[t]=t
	currtemp[t]=t

# A list of the structures at the temperatures along the rex simulation
structure=[]
# A list of the temperatures a given structure has during the simulation
tempstruc=[]

# Read the input from the rex output files.
input=commands.getoutput("grep ^the rexswap.log | awk '{print $8, $6, $7;}'").split('\n')

# For every line in the input
for line in input:
        (i,t1,t2)=map(int,line.split())
# if we have a new index
	if (i!=iold):
 		# append a series of structure arrays
		for ii in range(iold, i):
			structure.append(curr.copy())
			for t in temp:
				currtemp[curr[t]]=t
			tempstruc.append(currtemp.copy())
		# and set the old index to the current one
		iold=i
	# Update the structures
	tmp=curr[t1]
	curr[t1]=curr[t2]
	curr[t2]=tmp

# Add the current structures to the end
structure.append(curr.copy())
for t in temp:
	currtemp[curr[t]]=t
tempstruc.append(currtemp.copy())

# Print out all the structures:
for t in temp:
	for i in range(len(structure)):
		print tempstruc[i][t]
		#	print structure[i][t]
	print "&"




