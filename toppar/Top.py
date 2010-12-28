#!/usr/bin/env python
# fcp v0.1  01/11/2010

import copy,string

class Top(dict):
    """.
    """
    def __init__(self,iterator,name=''):
        """Constructor"""
        self.name       = name
        self.text       = [ line.lstrip() for line in iterator if ( not line.lstrip().startswith('!') and line.strip() ) ]
        self.check_me   = [ line.split('!')[0].lstrip().upper() for line in self.text ]
        # Initialize
        dict.__init__( self, self.split() )
        self.header     = self['header']
        del self['header']
        # Initialize
        self.charmm_header  = self.set_charmm_header()
        self.version        = self.set_version()
        self.links          = self.set_links()
        self.defaultPatch   = self.set_default_patch()
        self.patchRes       = self.set_patch_res()
        self.masses         = self.set_masses()
        if self.defaultPatch: self.set_patching()
#       del self.header

    def split(self):
        """Splits the file into n + 1 sections, where n is the number of residues in the file, and
        the +1 is the 'header', consisting of the charmm_header, version, masses and links.  This returns
        a tuple, with residue names as keys and 'header' as the +1 key.  This function has no buisness
        getting called by anything except the constructor."""
        keys    = []
        values  = []
        parsing_header  = True
        section_buffer  = []
        first_line      = 'header'
        for i,line in enumerate(self.check_me):
            section_buffer.append(self.text[i])
            if line.startswith('RESI') or line.startswith('PRES') or line.startswith('END'):
                if parsing_header:
                    parsing_header = False
                    keys.append('header')
                else:
                    keys.append(first_line.split()[1])
                values.append( [first_line] + section_buffer[:-1] )
                first_line      = self.text[i]
                section_buffer  = []
        return zip(keys,values)

    def set_charmm_header(self):
        """Return all lines that start with '*' __and__ have other text in them.
        The terminal header line is skipped!"""
        return [ line for line in self.header if ( line.startswith('*') and len( line.split() ) > 1 ) ]

    def set_version(self):
        """Return an integer denoting the version number of Charmm this topology file requires."""
        for line in self.header:
            try:
                return int( line.split()[0] )
            except ValueError:
                pass

    def set_masses(self):
        """Returns a list of lines of text with the mass information.  Duplicates are nuked,
        and masses are indexed properly."""
        masses = [ line for line in self.header if line.startswith('MASS') ]

        unique = {}                         # Nuke dupes
        for line in masses:
            unique[line.split()[2]] = line
        masses = [ unique[key] for key in unique.keys() ]

        def mass_index(string):             # Sort them by their original index
            return int( string.split()[1] )
        masses.sort(key=mass_index)

        buffer = []                         # Re-Index
        for i,mass in enumerate(masses):
            replace_me = ' %s ' % (mass.split()[1])
            buffer.append( mass.replace(replace_me, ' %i ' % (i+1)) )
        masses = buffer

        return masses

    def get_masses(self):
        """Returns a list of lines of text with the mass information.  Duplicates are nuked,
        and masses are indexed properly."""
        masses = [ line for line in self.masses if line.startswith('MASS') ]

        unique = {}                         # Nuke dupes
        for line in masses:
            unique[line.split()[2]] = line
        masses = [ unique[key] for key in unique ]

        def mass_index(string):             # Sort them by their original index
            return int( string.split()[1] )
        masses.sort(key=mass_index)

        buffer = []                         # Re-Index
        for i,mass in enumerate(masses):
            replace_me = ' %s ' % (mass.split()[1])
            buffer.append( mass.replace(replace_me, ' %i ' % (i+1)) )
        masses = buffer

        return masses

    def set_links(self):
        """Returns a list of lines of text which contains the link information.  Duplicates are nuked."""
        links = [ line.split()[1] + ' ' + str(i) for i,line in enumerate(self.header) if line.startswith('DECL') ]
        
        unique = {}                         # Nuke dupes
        for line in links:
            unique[line.split()[0]] = line
        links = [ unique[key] for key in unique ]

        def link_index(string):
            return int( string.split()[1] )

        links.sort(key=link_index)          # Sort
        links = [ 'DECL  %-s\n' % line.split()[0] for line in links ]

        return links

    def get_links(self):
        """Returns a list of lines of text which contains the link information.  Duplicates are nuked."""
        links = [ line.split()[1] + ' ' + str(i) for i,line in enumerate(self.links) if line.startswith('DECL') ]
        
        unique = {}                         # Nuke dupes
        for line in links:
            unique[line.split()[0]] = line
        links = [ unique[key] for key in unique ]

        def link_index(string):
            return int( string.split()[1] )

        links.sort(key=link_index)          # Sort
        links = [ 'DECL  %-s\n' % line.split()[0] for line in links ]

        return links

    def set_default_patch(self):
        """If default patching information exists, return it."""
        for line in self.header:
            if line.startswith('DEFA'):
                return line
        return None

    def merge(self,*args):
        """Takes two Top objects and merges them, returning another Top object.  Redundant residues, masses,
        links, etc are removed (quietly)."""
        new_top = copy.deepcopy(self)
        # Charmm header
        new_top.charmm_header = ['* This .rtf file was generated by CHARMMING, the original headers follow:\n']
        new_top.charmm_header += ['* name: %s\n' % self.name]
        new_top.charmm_header += self.charmm_header
        for arg in args:
            new_top.charmm_header += ['* name: %s\n' % arg.name]
            new_top.charmm_header += arg.charmm_header
        # Version
        for arg in args:
            if arg.version > new_top.version:
                new_top.version = arg.version
        # Masses
        for arg in args:
            new_top.masses += arg.masses
        new_top.masses = new_top.get_masses()
        # Links
        for arg in args:
            new_top.links += arg.links
        new_top.links = new_top.get_links()
        # Residues
        for arg in args:
            for key in arg:
                new_top[key] = arg[key]
        # Patch Res
        for arg in args:
            for key in arg.patchRes:
                new_top.patchRes[key] = arg.patchRes[key]

        new_top.text   = new_top.Print()
        new_top.check_me   = [ line.split('!')[0].lstrip().upper() for line in new_top.text ]
        return new_top

    def cgMerge(self,*args):
        """Takes two Top objects and merges them, returning another Top object.  Redundant residues, masses,
        links, etc are removed (quietly)."""
        new_top = copy.deepcopy(self)
        # Charmm header
        new_top.charmm_header = ['* This .rtf file was generated by CHARMMING, the original headers follow:\n']
        new_top.charmm_header += ['* name: %s\n' % self.name]
        new_top.charmm_header += self.charmm_header
        for arg in args:
            new_top.charmm_header += ['* name: %s\n' % arg.name]
            new_top.charmm_header += arg.charmm_header
        # Version
        for arg in args:
            if arg.version > new_top.version:
                new_top.version = arg.version
        # Masses
        for arg in args:
            new_top.masses += arg.masses
        new_top.masses = new_top.get_masses()
        # Links
        for arg in args:
            new_top.links += arg.links
        new_top.links = new_top.get_links()
        # Residues
        for arg in args:
            for key in arg:
                new_top[key] = arg[key]
        # Patch Res
        for arg in args:
            for key in arg.patchRes:
                new_top.patchRes[key] = arg.patchRes[key]

        new_top.text   = new_top.Print(cg=True)
        new_top.check_me   = [ line.split('!')[0].lstrip().upper() for line in new_top.text ]
        return new_top

    def cull(self,iterator):
        """Takes a list of residues, and eliminates all residues except those specified, then eliminates all
        masses except those in the specified residues."""
        new_top = copy.deepcopy(self)
        # Cull residues
        delete_keys = [ key for key in new_top if ( key not in iterator ) ]
        for key in delete_keys:
            del new_top[key]
        # TODO Need logic for patch culling -- talk to Rick
        # Cull masses from residues
        buffer_Res = [ line.split()[2] for values in new_top.values() for line in values if line.startswith('ATOM') ]
        # Cull masses from patches
        buffer_pRes = [ line.split()[2] for values in new_top.patchRes.values() for line in values if line.startswith('ATOM') ]
        buffer = list(set(buffer_Res + buffer_pRes))
        new_top.masses = [ line for line in new_top.masses if line.split()[2] in buffer ]
        new_top.masses = new_top.get_masses()
        new_top.text   = new_top.Print()
        new_top.check_me   = [ line.split('!')[0].lstrip().upper() for line in new_top.text ]
        return new_top

    def get_atom_list(self,res=''):
        """Prints a list of atom types that exist in the specified residues in the file.  Alternatively,
        if a residue is not specified, it prints a list of all atoms types in the entire file."""
        if not res:         # ie the default
            return [ line.split()[2] for line in self.masses ]
        else:
            try:
                return [ line.split()[2] for line in self[res] if line.startswith('ATOM') ]
            except KeyError:
                return [ line.split()[2] for line in self.patchRes[res] if line.startswith('ATOM') ]

    def get_atom_sets(self):
        """Returns a dictionary of sets.  One key/value pair per residue, where each key is the res name and
        each value is a set of atom types found within that specific residue.  This method is critical for 
        communication with the Par module."""
        atom_sets = {}
        wildcards = ['X']
        for res in self:
            atom_sets[res] = set( self.get_atom_list(res) + wildcards )
        for res in self.patchRes:
            atom_sets[res] = set( self.get_atom_list(res) + wildcards )
        return atom_sets

    def Print(self,cg=False):
        """Returns a list of strings that can be fed to a file handle to write a .rtf file."""
        char_score     = dict( [ (char,i+1) for i,char in enumerate(string.uppercase) ] )
        def res_score(arg):
            return ( char_score[ arg[0] ] * 1000 + int(arg[1:]) )

        print_me = []

        for line in self.charmm_header:
            print_me.append(line)
        print_me.append('*')
        print_me.append('\n')

        print_me.append('%i  1' % self.version)
        print_me.append('\n')

        for line in self.masses:
            print_me.append(line)
        print_me.append('\n')

        for line in self.links:
            print_me.append(line)
        print_me.append('\n')

        if cg:
            for key in sorted(self,key=res_score):
                for line in self[key]:
                    print_me.append(line)
                print_me.append('\n')
        else:
            for key in self:
                for line in self[key]:
                    print_me.append(line)
                print_me.append('\n')

        for key in self.patchRes:
            for line in self.patchRes[key]:
                print_me.append(line)
            print_me.append('\n')

        print_me.append('END')

        return print_me

    def set_patch_res(self):
        """Go through all residue information and partition RESI and PRES."""

        resi_buffer = []
        patch_buffer = {}

        for key in self:
            if self[key][0].split()[0] == 'RESI':
                resi_buffer.append(key)
            elif self[key][0].split()[0] == 'PRES':
                patch_buffer[key] = self[key]
            else:
                raise AssertionError
        
        buffer = self.keys()
        for res in buffer:
            if res not in resi_buffer:
                del self[res]
        
        return patch_buffer

    def set_patching(self):
        """If default patch info was specified in the header, apply this to all residues without explicit patch info."""
        patch_line = ' '.join( ['PATCHING '] + self.defaultPatch.split()[1:] + ['\n'] )
        for res in self:

            for line in self[res]: # If patch info is present, skip current residue
                if line.startswith('PATC'):
                    break
            
            self[res].append(patch_line)
        return
