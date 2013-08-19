#!/usr/bin/env python

import sys, math

def print_connect(bl):

    #sys.stdout.write('CONNECT')

    mydict = {}

    for elt in bl:
        if elt[0] in mydict.keys():
            mydict[elt[0]].append(elt[1])
        else:
            mydict[elt[0]] = [elt[1],]

    for k in mydict.keys():
        sys.stdout.write('CONNECT %d' % k)
        for x in mydict[k]:
            sys.stdout.write(' %d' % x)
        sys.stdout.write('\n')

if __name__ == '__main__':

    if len(sys.argv) < 2:
        sys.stderr.write('Usage: psf2connect.py <PSF file>\n')
        sys.exit(-1)

    bondlist = []
    psffp = open(sys.argv[1],'r')

    l1 = psffp.readline()
    if not l1.startswith('PSF'):
        sys.stderr.write('Malformed PSF!\n')
        sys.exit(-1)

    junk = psffp.readline()
    l2 = psffp.readline()
    ntitl = int(l2.split()[0])
    sys.stderr.write('Skipping %d title lines.\n' % ntitl)

    # also skip the blank line after the title
    for x in range(ntitl+1):
        junk = psffp.readline()

    latom = psffp.readline()
    natom = int(latom.split()[0])
    sys.stderr.write('Got %d atoms, but we don\'t care about them.\n' % natom)

    # again, skip over the blank line at the end
    for x in range(natom+1):
        junk = psffp.readline()
    
    lbond = psffp.readline()
    nbond = int(lbond.split()[0])
    if nbond <= 0:
        sys.stderr.write('HEY! There are no bonds here!\n')
        sys.exit(-1)
    sys.stderr.write('Getting the information for %d bonds.\n' % nbond)

    nbread = 0
    while nbread < nbond:
        line = psffp.readline()
        larr = line.split()
        if len(larr) % 2 != 0:
            sys.stderr.write('HEY! Incomplete bonds on line. Aborting.\n')
            sys.exit-1

        for i in range(0,len(larr),2):
            atom1, atom2 = int(larr[i]), int(larr[i+1])
            if atom1 > natom or atom1 < 0:
                sys.stderr.write('HEY! Bad atom number %d.\n' % atom1)
                sys.exit(-1)
            if atom2 > natom or atom2 < 0:
                sys.stderr.write('HEY! Bad atom number %d.\n' % atom1)
                sys.exit(-1)

            # in order
            if atom2 < atom1:
                atom1, atom2 = atom2, atom1
            bondlist.append((atom1,atom2))

        nbread += len(larr)/2

    psffp.close()
    print_connect(bondlist)
