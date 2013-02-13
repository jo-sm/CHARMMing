#!/usr/bin/env python

# usage example:
# ./rename_dupes.py -n antechamber.rtf -b antechamber_orig /usr/local/charmming/toppar/top_all36_prot.rtf /usr/local/charmming/toppar/top_all36_cgenff.rtf


import sys,os,shutil
from optparse import OptionParser

def modify_prm(d,f):

   fp = open(f,'r')
   lines = fp.readlines()
   fp.close()

   fp = open(f,'w')
   for line in lines:
       line = line.upper()
       if line.startswith('MASS'):
           continue
       for an in d.keys():
           if '%s ' % an in line:
               line = line.replace('%s ' % an,'%s ' % d[an])
       fp.write(line)
   fp.close()

def build_dict(alist,f):

    atmnum_dict = {}
    rename_dict = {}
    lines = []

    # make the rename dict
    fp = open(f,'r')
    if not fp:
        raise()

    for line in fp:
        line = line.upper()
        lines.append(line)
        line = line.strip()
        if not line.startswith('MASS'):
            continue
 
        mynum = int(line.split()[1])
        myatc = line.split()[2]
        orignum = mynum
        origatc = myatc

        # find a unique atom number
        good = False
        while not good:
            good = True
            for elt in alist:
                if mynum == elt[0]:
                    good=False
                    mynum+=1
                    break

        # find a unique atom type code
        good = False
        while not good:
            good = True
            for elt in alist:
                if myatc == elt[1]:
                    good = False
                    if myatc == origatc and len(myatc) < 5:
                        myatc+='A'
                    else:
                        newchar=ord(myatc[-1])+1
                        myatc[-1]=newchar
                    break

        if orignum != mynum:
            atmnum_dict[orignum]=mynum
        if origatc != myatc:
            rename_dict[origatc]=myatc

    fp.close()

    # if we found any atoms to rename or retype, overwrite the
    # existing file.
    if atmnum_dict or rename_dict:
        fp = open(f,'w')
        if not fp:
            raise

        for line in lines:
            if line.startswith('MASS'):
                if int(line.split()[1]) in atmnum_dict.keys():
                    line = line.replace(' %s ' % line.split()[1],' %s ' % atmnum_dict[int(line.split()[1])])
                if line.split()[2] in rename_dict.keys():
                    line = line.replace('%-5s' % line.split()[2],'%-5s' % rename_dict[line.split()[2]])
            elif line.startswith('ATOM'):
                if line.split()[2] in rename_dict.keys():
                    line = line.replace('%-5s' % line.split()[2],'%-5s' % rename_dict[line.split()[2]])

            fp.write(line)

        fp.close()

    return rename_dict

def build_list(d,f):

    fp = open(f,'r')
    if not fp:
        raise

    for line in fp:
        line = line.strip()
        line = line.upper()
        if not line.startswith('MASS'):
            continue
        mynum = int(line.split()[1])
        atmtuple = tuple((mynum,line.split()[2]))
        if atmtuple in d:
            sys.stderr.write('Multiple definitions of atom number %d with code %s.\n' % (atmtuple[0],atmtuple[1]))
        else:
            for tp in d:
                if atmtuple[0] == tp[0]:
                    sys.stderr.write('ERROR: atom number %d is already used!\n' % tp[0])
                    raise
                if atmtuple[1] == tp[1]:
                    sys.stderr.write('ERROR: atom type code %s is already used!\n' % tp[1])
                    raise

        d.append(atmtuple)

    fp.close()

if __name__ == '__main__':

    atmlist = []

    parser = OptionParser()
    parser.add_option('-v','--verbose',action='store_true',default=False,help='Unpleasantly verbose print-out')
    parser.add_option('-n','--new',metavar='FILE',help='New RTF file whose duplicate names should be changed (NB must end in .rtf)')
    parser.add_option('-b','--backup',metavar='STRING',help='Make a backup starting with STRING')
    (opts,args) = parser.parse_args(sys.argv[1:])

    if not opts.new:
        sys.stderr.write('You must provide a newly-generated RTF to work on!\n')
        sys.exit(-1)
    if not args:
        sys.stderr.write('No existing topology files provided. This is a no-op.\n')
        sys.exit(-1)
    if not opts.new.endswith('.rtf'):
        sys.stderr.write("Your RTF file must end with .rtf. Plus, you're ugly!\n")
        sys.exit(-1)

    basename = opts.new.replace('.rtf','')
    for file in args:
        if opts.verbose:
            print 'Looking for atom types in %s' % file
        ##try:
        build_list(atmlist,file)
        ##except:
        ##    sys.stderr.write('Problem with file %s ... is it readable?\n' % file)
        ##    sys.exit(-2)

    if opts.verbose:
        for elt in atmlist:
            print 'Atom with index %d has type code %s.' % (elt[0],elt[1])
    if opts.backup:
        shutil.copy(opts.new,opts.backup + '.rtf')
    rename_dict = build_dict(atmlist,opts.new)
    if not rename_dict:
        sys.stderr.write('Did not find any atoms that needed to be renamed. Finishing up.')
        sys.exit(0)
    if opts.verbose:
        for k in rename_dict.keys():
            print 'Renaming atom type %s to %s' % (k,rename_dict[k])

    prmfile = basename + '.prm'
    if opts.backup:
        shutil.copy(prmfile,opts.backup + '.prm')
    modify_prm(rename_dict,prmfile)

    pdbfile = basename + '.pdb'
    if opts.backup:
        shutil.copy(pdbfile,opts.backup + '.pdb')
    for k in rename_dict.keys():
        os.system("sed -i -e 's/ %s / %s /ig' %s" % (k,rename_dict[k],pdbfile))
