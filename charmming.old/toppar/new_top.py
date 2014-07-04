#!/usr/bin/env python
# fcp v0.2 06/01/2010
# req python v2.5

import lib.Etc as Etc

class Top(object):
    """."""
    def __init__(self,iterator,name=''):
        """Constructor"""
        self.name       = name
        self.text       = Etc.stripBare(iterator,'!')
        self.__split    = list( Etc.paragraphs(self.text,['RESI','PRES','END']) )
        # Split file into: Header/Residue/Patch Residue
        self.header     = self.__split[0]
        self.resi       = [ paragraph for paragraph in self.__split[1:] if paragraph[0][:4] == 'RESI' ]
        self.pres       = [ paragraph for paragraph in self.__split[1:] if paragraph[0][:4] == 'PRES' ]
        # Parse Header
        self.charmmHeader = [ line for line in self.header if line.startswith('*') and len( line.split() ) > 1 ]
        self.version    = self.__version()
        self.mass       = [ ' '.join(line.split()[:4]) for line in self.header if line.startswith('MASS') ]
        self.link       = [ ' '.join(line.split()[:2]) for line in self.header if line.startswith('DECL') ]
        self.defaultPatch = self.__defaultPatch()
        self.autoGenerate = self.__autoGenerate() 
        # Lists
        self.linkList   = self.get_linkList()
        # Dictionaries
        self.resiDict   = self.get_resiDict()
        self.presDict   = self.get_presDict()
        self.massDict   = self.get_massDict()
        # Regenerate 
        self.regenerate_text()

    def __len__(self):
        """Returns number of residues in top file."""
        return len(self.resiDict)

    def __version(self):
        """Called by the constructor to parse charmm version info."""
        try:
            version = [ line.split()[0] for line in self.header if line.startswith( tuple(map(str,range(10))) ) ]
            version.sort()
            return version[-1]
        except:
            return -1

    def __defaultPatch(self):
        """Called by the constructor to parse default patch info."""
        for line in self.header:
            if line.startswith('DEFA'):
                return line
        return None

    def __autoGenerate(self):
        """Called by the constructor to parse auto generate info."""
        for line in self.header:
            if line.startswith('AUTO'):
                return line
        return None

    def apply_defaultPatch(self):
        """Apply default patching options to all resi without explicit patching options, and remove default."""
        # Can only apply default patch if it exists!
        if not self.defaultPatch: return
        default = ' '.join( ['PATCHING'] + self.defaultPatch.split()[1:] )
        # Determine if each residue has explicit patching info
        for key in self.resiDict.keys():
            patching = False
            for line in self.resiDict[key]:
                if line.startswith('PATC'):
                    patching = True
                    break
            # If it doesnt, add the default patching line
            if not patching:
                self.resiDict[key] += [default]
        # Remove the default information, since all residues now have explicit patching info
        self.defaultPatch = None

    def get_linkList(self):
        """."""
        return list( set( [ line.split()[1] for line in self.link ] ) )

    def get_resiDict(self):
        """Build a dictionary of all residues present."""
        return dict(zip( [ resi[0].split()[1] for resi in self.resi ] ,self.resi))

    def get_presDict(self):
        """Build a dictionary of all patch residues present."""
        return dict(zip( [ pres[0].split()[1] for pres in self.pres ] ,self.pres))

    def get_massDict(self):
        """Build a dictionary of all masses present."""
        return dict(zip( [ mass.split()[2] for mass in self.mass ] ,self.mass))

    def regenerate_text(self):
        """Rebuild raw text from header, resi and pres info."""
        buffer = []
        self.regenerate_header()
        buffer += self.header
        self.regenerate_resi()
        buffer += self.resi
        self.regenerate_pres()
        buffer += self.pres
        buffer += ['END']
        self.text = Etc.flatten(buffer)

    def regenerate_header(self):
        """Rebuild header text from charmm header, version, mass, link and default sections"""
        buffer = []
        buffer += self.charmmHeader
        buffer += ['*']
        buffer += ['  %s  1' % self.version]
        buffer += [' ']
        self.regenerate_mass()
        buffer += self.mass
        buffer += [' ']
        self.regenerate_link()
        buffer += self.link
        buffer += [' ']
        if self.defaultPatch: buffer += [self.defaultPatch,' ']
        if self.autoGenerate: buffer += [self.autoGenerate,' ']
        self.header = buffer

    def regenerate_link(self):
        """Rebuild link text from link info."""
        # Sort by #/+/- and alpha
        buffer = []
        for symbol in ['#','+','-']:
            tmpLink = [ taco[1:] for taco in self.linkList if taco.startswith(symbol) ]
            tmpLink.sort()
            tmpLink = [ symbol + taco for taco in tmpLink ]
            buffer.append(tmpLink)
        buffer = Etc.flatten(buffer)
        # Formatting
        buffer2 = []
        for line in buffer:
            new_line = ['DECL'] + [line]
            new_line = '%-5s%6s' % tuple(new_line)
            buffer2.append(new_line)
        self.link = buffer2

    def regenerate_resi(self):
        """Rebuild resi text from residue info."""
        buffer = []
        for value in self.resiDict.values():
            buffer2 = value + [' ']
            buffer += [buffer2]
        # Sort by residue name by alpha
        def resi_index(value):
            return value[0].split()[1]
        buffer.sort(key=resi_index)
        self.resi = buffer

    def regenerate_pres(self):
        """Rebuild pres text from patch residue info."""
        buffer = []
        for value in self.presDict.values():
            buffer2 = value + [' ']
            buffer += [buffer2]
        # Sort by residue name by alpha
        def pres_index(value):
            return value[0].split()[1]
        buffer.sort(key=pres_index)
        self.pres = buffer

    def regenerate_mass(self):
        """Rebuild mass text from mass info. """
        buffer = self.massDict.values()
        # Sort by original indexing
        def mass_index(value):
            return int( value.split()[1] )
        buffer.sort(key=mass_index)
        # Formatting
        buffer2 = []
        for i,line in enumerate(buffer):
            args = tuple( [i+1] + line.split()[2:] )
            buffer2.append('MASS%6d%6s%12s' % args)
        self.mass = buffer2

    def cull(self,iterator):
        """."""
        self.apply_defaultPatch()
        self.cull_resi(iterator)
        self.cull_pres()
        self.cull_mass()
        self.cull_link()
        self.regenerate_text()

    def cull_resi(self,iterator):
        """Takes a list of residues and eliminates all residues save those explicitly specified."""
        top_set     = set( self.resiDict.keys() )   # residues in the topology file
        pdb_set     = set( iterator )               # residues in the pdb file
        cull_set    = top_set - pdb_set             # top .not. pdb
        # Throw an error if you have a residue in the pdb but not in the topology file
        if not pdb_set <= top_set:
            error_set = pdb_set - top_set
            error_string = ' '.join( list( error_set ) )
            print 'cull_resi: The following residues are to be culled, yet are not present in the topology file: %s\n' % error_string
            return
        # If culling would result in an empty top file, dont cull!
        if not cull_set:
            print 'cull_resi: This cull command would result in a null topology file, aborting\n'
            return
        # Cull & Regenerate
        for key in cull_set: del self.resiDict[key]
        self.regenerate_resi()

    def cull_pres(self):
        """Eliminates extraneous patching residues. Based upon residues present in pdb and patching info in top file."""
        # Figure out which patching residues should stay
        buffer = [ taco.split() for taco in Etc.flatten(self.resi) if taco.startswith('PATC') ]
        try:
            buffer += self.defaultPatch.split()
        except: pass
        buffer = list( set(Etc.flatten(buffer)) )
        buffer = [ taco for taco in buffer if not taco.startswith( ('NONE','FIRS','LAST','PATC','DEFA') ) ]
        ###
        top_set     = set( self.presDict.keys() )   # pResidues in the topology file
        pdb_set     = set( buffer )                 # pResidues required for the pdb file
        cull_set    = top_set - pdb_set             # top .not. pdb
        # Throw an error if you have a residue in the pdb but not in the topology file
        if not pdb_set <= top_set:
            error_set = pdb_set - top_set
            error_string = ' '.join( list( error_set ) )
            print 'cull_resi: The following residues are to be culled, yet are not present in the topology file: %s\n' % error_string
            return
        # Cull & Regenerate
        for key in cull_set: del self.presDict[key]
        self.regenerate_pres()

    def cull_mass(self):
        """Eliminates extraneous masses.  Based upon residues and patch residues present in pdb and top file."""
        # Figure out which masses should stay
        resi_buffer = [ line.split()[2] for line in Etc.flatten(self.resi) if line.startswith('ATOM') ]
        pres_buffer = [ line.split()[2] for line in Etc.flatten(self.pres) if line.startswith('ATOM') ]
        buffer = list( set(resi_buffer + pres_buffer) )
        ###
        top_set     = set( [ line.split()[2] for line in self.mass ] )  # Atom types from the top file
        pdb_set     = set( buffer )                 # Atom types required for the pdb file
        cull_set    = top_set - pdb_set             # top .not. pdb
        # Throw an error if you have a mass in the pdb but not in the topology file
        if not pdb_set <= top_set:
            error_set = pdb_set - top_set
            error_string = ' '.join( list( error_set ) )
            print 'cull_mass: The following masses are to be culled, yet are not present in the topology file: %s\n' % error_string
            return
        # Cull & Regenerate
        for key in cull_set: del self.massDict[key]
        self.regenerate_mass()

    def cull_link(self):
        """Eliminates extraneous DECL statements.  Based upon residues and patch residues present in pdb and top file."""
        # Figure out which links should stay
        resi_buffer = [ line.split()[1:] for line in Etc.flatten(self.resi) ]
        pres_buffer = [ line.split()[1:] for line in Etc.flatten(self.pres) ]
        buffer = [ taco for taco in Etc.flatten(resi_buffer + pres_buffer) if taco.startswith( ('#','+','-') ) ]
        def isnt_number(string):
            try:
                float(string)
                return False
            except ValueError:
                return True
        buffer = list( set( filter(isnt_number,buffer) ) )
        ###
        top_set     = set( self.linkList )          # Atom labels from the top file
        pdb_set     = set( buffer )                 # Atom labels required for the pdb file
        cull_set    = top_set - pdb_set             # top .not. pdb
        # Throw an error if you have a mass in the pdb but not in the topology file
        if not pdb_set <= top_set:
            error_set = pdb_set - top_set
            error_string = ' '.join( list( error_set ) )
            print 'cull_link: The following links are to be culled, yet are not present in the topology file: %s\n' % error_string
            return
        # Cull & Regenerate
        self.linkList = buffer
        self.regenerate_link()

    def Print(self):
        """Print topology file."""
        def unchomp(line):
            return '%s\n' % line
        return map(unchomp,self.text)

    def Write(self,file):
        """Writes topology file."""
        write_to = open(file,'w')
        for line in self.Print():
            write_to.write(line)
        write_to.close()

    def merge(self,*args):
        """Merge 2 or more topology files."""
        header  = self.header
        resi    = self.resi
        pres    = self.pres
        try:
            for top in args:
                assert type(top) == Top
                header += top.header
                resi   += top.resi
                pres   += top.pres
            buffer = Etc.flatten(header + resi + pres)
            return Top(buffer)
        except AssertionError:
            print 'merge: Only valid Top objects may be merged\n'
            print 'merge: "%s" is of the type: %s\n' % (top,type(top))

    def __add__(self,other):
        """Merge 2 topology files."""
        return self.merge(other)

    def export_atoms(self,wildCard = ['X']):
        """Main form of communication with Par module.  Returns a dictionary of lists which indicate atom groups for parameter culling."""
        pass
