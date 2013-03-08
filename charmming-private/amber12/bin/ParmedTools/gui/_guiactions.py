"""
This module contains all of the methods to invoke specific ParmedActions listed
in that module. Each method will establish the necessary variables (StringVar),
create whatever window it needs to get information from the user (or give
information to the user), wait for that window to close, then dispatch that
Action to the class in ParmedActions. 

Follow the general trend if you wish to add your method to the GUI. Note, any
method that you want accessible through the GUI must have an action method put
here with the same name as the class found in ParmedActions.
"""

from Tkinter import *
from tkMessageBox import askyesno, showinfo, showerror
from ParmedTools import ParmedActions
from ParmedTools.gui.guifiletools import file_chooser, save_file_chooser
from ParmedTools.gui import _guiwidgets

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def writefrcmod(amber_prmtop, messages):
   """ Dumps an frcmod file to a given filename """
   fname = save_file_chooser('frcmod', '.frcmod')
   if fname: 
      action = ParmedActions.writefrcmod(amber_prmtop, fname)
      action.execute()
      messages.add_line(str(action))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def loadrestrt(amber_prmtop, messages):
   """ Finds a file to load as a restart """
   type_list = [('Inpcrd', '*.inpcrd'), ('Inpcrd', '*.crd'),
                ('Restart', '*.restrt'), ('Restart', '*.rst7'),
                ('All Files', '*')]
   fname = file_chooser('Amber Coordinate File', type_list)
   if fname: 
      action = ParmedActions.loadrestrt(amber_prmtop, fname)
      action.execute()
      messages.add_line(str(action))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def writeoff(amber_prmtop, messages):
   """ Dumps an OFF library to a given filename """
   fname = save_file_chooser('OFF library', '.lib')
   if fname: 
      action = ParmedActions.writeoff(amber_prmtop, fname)
      action.execute()
      messages.add_line(str(action))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def changeradii(amber_prmtop, messages):
   """ Allows users to change the GB Radius set """
   title = 'Choose a Radius Set'
   desc = ('Select a radius set for implicit solvent\n' +
           'calculations. This has the same effect as\n' +
           '"set default PBRadii <value>" in tleap')
   radius_selection = StringVar()
   namelist = ['bondi', 'mbondi', 'mbondi2', 'mbondi3', 'amber6']
   cmd_window = _guiwidgets.RadioButtonWindow(
                  amber_prmtop, title, desc, radius_selection, namelist)
   # Wait until the window is destroyed, then get the variable and pass it
   # over to the class
   cmd_window.wait_window()
   sel = str(radius_selection.get())
   if sel:
      action = ParmedActions.changeradii(amber_prmtop, sel)
      action.execute()
      messages.add_line(str(action))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def changeljpair(amber_prmtop, messages):
   """ Changes a pair-wise LJ interaction for pre-combined epsilon/Rmin """
   # The variables we need for changeljpair
   widget_list = [('MaskEntry', 'Atom(s) Type 1 Mask'),
                  ('MaskEntry', 'Atom(s) Type 2 Mask'),
                  ('Entry', 'Combined Radius'),
                  ('Entry', 'Combined Well Depth')]
   # Variable list -- we need 2 masks and 2 floats
   var_list = [StringVar(), StringVar(), StringVar(), StringVar()]
   # description
   description = ' '.join(ParmedActions.changeljpair.__doc__.split())
   cmd_window = _guiwidgets.ActionWindow('changeLJPair', amber_prmtop,
                             widget_list, var_list, description)
   cmd_window.wait_window()
   # Make sure we didn't cancel (or just press OK with no input), or just leave
   vars_exist = True in [bool(v.get()) for v in var_list]
   if not vars_exist: return
   # Now that we did something, do it
   var1, var2 = var_list[0].get(), var_list[1].get()
   var3, var4 = var_list[2].get(), var_list[3].get()
   try:
      action = ParmedActions.changeljpair(amber_prmtop, var1, var2, var3, var4)
      action.execute()
      messages.add_line(str(action))
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def outparm(amber_prmtop, messages):
   """ Output a final topology file """
   fname = save_file_chooser('prmtop', '.prmtop')
   if fname: 
      action = ParmedActions.outparm(amber_prmtop, fname)
      action.execute()
      messages.add_line(str(action))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printflags(amber_prmtop, messages):
   """ Prints all of the flags in the topology file """
   window = Toplevel()
   window.resizable(width=False, height=False)
   window.title('%%FLAG list in %s' % amber_prmtop)
   text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                        pady=5, width=80, height=20)
   text.pack()
   text.grid(row=0, column=0)
   action = ParmedActions.printflags(amber_prmtop)
   text.text.insert(END, str(action))
   window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printpointers(amber_prmtop, messages):
   """ Prints all of the flags in the topology file """
   window = Toplevel()
   window.resizable(width=False, height=False)
   window.title('POINTER list in %s' % amber_prmtop)
   text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                        pady=5, width=80, height=20)
   text.pack()
   text.grid(row=0, column=0)
   action = ParmedActions.printpointers(amber_prmtop)
   text.text.insert(END, str(action))
   window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def changelj14pair(amber_prmtop, messages):
   """ Changes specific 1-4 Lennard Jones pairs """
   # Only good for chamber topologies
   if not amber_prmtop.chamber:
      showerror('Incompatible',
                'changeLJ14Pair is only valid for chamber topologies!')
      return
   # variables we need for changelj14pair
   widget_list = [('MaskEntry', 'Atom(s) Type 1 Mask'),
                  ('MaskEntry', 'Atom(s) Type 2 Mask'),
                  ('Entry', 'Combined Radius'),
                  ('Entry', 'Combined Well Depth')]
   # Variable list -- we need 2 masks and 2 floats
   var_list = [StringVar(), StringVar(), StringVar(), StringVar()]
   # description
   description = ' '.join(ParmedActions.changeljpair.__doc__.split())
   cmd_window = _guiwidgets.ActionWindow('changeLJ14Pair', amber_prmtop,
                             widget_list, var_list, description)
   cmd_window.wait_window()
   # Make sure we didn't cancel (or just press OK with no input), or just leave
   vars_exist = True in [bool(v.get()) for v in var_list]
   if not vars_exist: return
   # Now that we did something, do it
   var1, var2 = var_list[0].get(), var_list[1].get()
   var3, var4 = var_list[2].get(), var_list[3].get()
   try:
      action=ParmedActions.changelj14pair(amber_prmtop, var1, var2, var3, var4)
      action.execute()
      messages.add_line(str(action))
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def checkvalidity(amber_prmtop, messages):
   """ Basic validity checks """
   action = ParmedActions.checkvalidity(amber_prmtop)
   valid = action.execute()
   messages.add_line(str(action))

   if not valid:
      showerror('Bad Topology!', 'Problems found. Check terminal or your ' +
                'standard error output for details')
   else:
      showinfo('Check passed', '%s looks OK to me.' % amber_prmtop)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def change(amber_prmtop, messages):
   """ Allows us to change a specific atomic property """
   # The spinbox is sent with the Spinbox, label, and then a list of all of the
   # values to give to it
   widget_list = [('Spinbox', 'Property to change', 'CHARGE', 'MASS',
                   'RADII', 'SCREEN', 'ATOM_NAME', 'AMBER_ATOM_TYPE'),
                  ('MaskEntry', 'Atoms to change'),
                  ('Entry', 'New Value for Property')]
   # We need 3 string variables, then get the description
   var_list = [StringVar(), StringVar(), StringVar()]
   description = 'Changes the property of given atoms to a new value'
   # Create the window, open it, then wait for it to close
   cmd_window = _guiwidgets.ActionWindow('change', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   vars_found = True in [bool(v.get()) for v in var_list]
   if not vars_found: return
   # If we did, store them and pass it to the class
   var1, var2, var3 = var_list[0].get(), var_list[1].get(), var_list[2].get()
   try:
      action = ParmedActions.change(amber_prmtop, var1, var2, var3)
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
   messages.add_line(str(action))
   action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printinfo(amber_prmtop, messages):
   """ Prints all of the info in a given FLAG """
   # Unfortunately, since printinfo was one of the earlier functions added, it
   # was implemented by literally writing everything to sys.stdout. To minimize
   # changes needed in ParmedActions, we create a string-like class here that
   # has a "write" method, and change printinfo's outfile to point to this
   # class instead
   class writable_str(list):
      def __init__(self, mystr): list.__init__(self, mystr)
      def write(self, addition): self.extend(list(addition))
      def __str__(self): return ''.join(self)

   # Set up the window
   # variables we need for printInfo
   widget_list = [('Entry', '%FLAG you want info from')]
   # Variable list -- we need a single string
   var_list = [StringVar()]
   # description
   description = ' '.join(ParmedActions.printinfo.__doc__.split())
   cmd_window = _guiwidgets.ActionWindow('printInfo', amber_prmtop,
                             widget_list, var_list, description)
   cmd_window.wait_window()
   # Make sure we didn't cancel (or just press OK with no input), or just leave
   var = var_list[0].get()
   if not var: return
   # Now that we did something, do it
   return_chars = writable_str('')
   ParmedActions.printinfo.outfile = return_chars
   action = ParmedActions.printinfo(amber_prmtop, var)
   if not action.found:
      showerror('Not Found!', '%%FLAG %s not found!' % var.upper())
      return

   window = Toplevel()
   window.resizable(width=False, height=False)
   window.title('POINTER list in %s' % amber_prmtop)
   text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                        pady=5, width=85, height=20)
   text.pack()
   text.grid(row=0, column=0)
   text.text.insert(END, str(return_chars))
   window.wait_window()
      
   messages.add_line(str(action))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def addljtype(amber_prmtop, messages):
   """ Turns given mask into a new LJ atom type """
   # We need a mask, new radius, new epsilon, and for chamber topologies,
   # a new radius-1-4 and epsilon-1-4. Don't add the latter ones until we
   # know if we have a chamber prmtop or not.
   widget_list = [('MaskEntry', 'Atoms to make new LJ Type'),
                  ('Entry', 'New Radius (default old radius)'),
                  ('Entry', 'New Depth (default old depth)')]
   # We need 5 string variables, then get the description. 
   var_list = [StringVar(), StringVar(), StringVar(), StringVar(), StringVar()]
   description=('Turns given mask into a new LJ atom type. Uses the radius\n' +
                'and well depth from the first atom type in <mask> if none\n' +
                'are provided.')
   if amber_prmtop.chamber:
      widget_list += [('Entry', 'New Radius for 1-4 Terms'),
                      ('Entry', 'New Depth for 1-4 Terms')]
   # Create the window, open it, then wait for it to close
   cmd_window = _guiwidgets.ActionWindow('addLJType', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   vars_found = True in [bool(v.get()) for v in var_list]
   if not vars_found: return
   # addljtype expects any _non_specified variables to be None
   var_list = [v.get() for v in var_list]
   for i, v in enumerate(var_list):
      if not v: var_list[i] = None
   # If we did, store them and pass it to the class
   var1, var2, var3, var4, var5 = var_list[0], var_list[1], var_list[2], \
                                  var_list[3], var_list[4]
   try:
      action = ParmedActions.addljtype(amber_prmtop,var1,var2,var3,var4,var5)
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
   messages.add_line(str(action))
   action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def setmolecules(amber_prmtop, messages):
   """ Sets the molecules. Asks users if they want ions to be solute or not """
   title = 'Set Molecularity'
   desc = ('This will determine the molecular topology (ATOMS_PER_MOLECULE\n' +
           'and SOLVENT_POINTERS). Do you want the ions to be considered\n' +
           'part of the solute?')
   solute_ions = StringVar()
   namelist = ['Yes', 'No']
   cmd_window = _guiwidgets.RadioButtonWindow(
                  amber_prmtop, title, desc, solute_ions, namelist)
   # Wait until the window is destroyed, then get the variable and pass it
   # over to the class
   cmd_window.wait_window()
   sel = str(solute_ions.get())
   if not sel: return
   if sel == 'Yes': sel = 'True'
   else: sel = 'False'
   if sel:
      action = ParmedActions.setmolecules(amber_prmtop, 'solute_ions='+sel)
      action.execute()
      messages.add_line(str(action))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printdetails(amber_prmtop, messages):
   """ Prints details about a given Amber mask """
   title = 'Print Details'
   mask = StringVar()
   cmd_window = Toplevel()
   cmd_window.title(title)
   mask_entry = _guiwidgets.MaskEntry(cmd_window, amber_prmtop,
                     'Input an Amber Mask', mask, cmd_window)
   mask_entry.config(pady=10)
   mask_entry.pack()
   mask_entry.grid(row=0, column=0, sticky=N+E+S+W)
   button = Button(cmd_window, text='OK / Quit', command=cmd_window.destroy)
   button.pack()
   button.grid(row=1, column=0, sticky=N+E+S+W)
   cmd_window.wait_window()
   if not mask.get(): return
   # Print our mask
   window = Toplevel()
   window.resizable(width=False, height=False)
   window.title('Atom information for mask %s' % mask.get())
   text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                        pady=5, width=100, height=20)
   text.pack()
   text.grid(row=0, column=0)
   action = ParmedActions.printdetails(amber_prmtop, mask.get())
   text.text.insert(END, str(action))
   window.wait_window()
   messages.add_line('Printed Amber Mask details on [%s]' % mask.get())

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def netcharge(amber_prmtop, messages):
   """ Calculates the net charge and shows its value """
   title = 'Net Charge'
   mask = StringVar()
   cmd_window = Toplevel()
   mask_entry = _guiwidgets.MaskEntry(cmd_window, amber_prmtop,
                     'Mask from which to calculate charge', mask, cmd_window)
   mask_entry.config(pady=10)
   mask_entry.pack()
   mask_entry.grid(row=0, column=0, sticky=N+E+S+W)
   button = Button(cmd_window, text='OK / Quit', command=cmd_window.destroy)
   button.pack()
   button.grid(row=1, column=0, sticky=N+E+S+W)
   cmd_window.wait_window()
   if not mask.get(): return
   action = ParmedActions.netcharge(amber_prmtop, mask.get())
   chg = action.execute()
   showinfo('Net Charge', 'The net charge of [%s] is %.4f' % (mask.get(), chg))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def strip(amber_prmtop, messages):
   """ Strips a mask from the topology file """
   if amber_prmtop.chamber:
      showerror('Not Implemented', 'The strip command does not yet work with ' +
                'chamber topologies')
      return
   # We need a mask, new radius, new epsilon, and for chamber topologies,
   # a new radius-1-4 and epsilon-1-4. Don't add the latter ones until we
   # know if we have a chamber prmtop or not.
   widget_list = [('MaskEntry', 'Atoms to strip from topology')]
   # We need 5 string variables, then get the description. 
   var_list = [StringVar()]
   description=('Strips the selected atoms from the topology file. All\n' +
                'remaining atoms and parameters remain unchanged. Any\n' +
                'parameters associated with stripped atoms are removed.')
   # Create the window, open it, then wait for it to close
   cmd_window = _guiwidgets.ActionWindow('addLJType', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   var = var_list[0]
   if not var.get(): return
   try:
      action = ParmedActions.strip(amber_prmtop, var.get())
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
   messages.add_line(str(action))
   action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def addexclusions(amber_prmtop, messages):
   """ Adds atoms to other atoms' exclusion list """
   # We need 2 masks
   widget_list = [('MaskEntry', 'Atoms to add excluded atoms to'),
                  ('MaskEntry', 'Atoms to exclude from other mask')]
   # We have 2 mask variables
   var_list = [StringVar(), StringVar()]
   # Description
   description=('Allows you to add arbitrary excluded atoms to exclusion\n' +
       'lists. This omits all non-bonded interactions in the direct-space\n' +
       'calculation, but does not omit interactions from adjacent cells in\n' +
       'periodic simulations')
   cmd_window = _guiwidgets.ActionWindow('addExclusions', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   
   # Bail out if we didn't get any variables
   if not True in [bool(v.get()) for v in var_list]: return
   
   var1, var2 = var_list[0].get(), var_list[1].get()
   try:
      action = ParmedActions.addexclusions(amber_prmtop, var1, var2)
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
   messages.add_line(str(action))
   action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def changeprotstate(amber_prmtop, messages):
   # We need a mask and a state
   widget_list = [('MaskEntry', 'Residue to change protonation state'),
                  ('Entry', 'Protonation state to change to')]
   var_list = [StringVar(), StringVar()]
   description=('Changes the protonation state of a titratable residue that\n' +
         'can be treated with constant pH MD in Amber. The residue must be\n' +
         'defined in $AMBERHOME/AmberTools/src/etc/cpin_data.py')
   cmd_window = _guiwidgets.ActionWindow('changeProtState', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   
   # Bail out if we didn't get any variables
   if not True in [bool(v.get()) for v in var_list]: return
   
   var1, var2 = var_list[0].get(), var_list[1].get()
   try:
      action = ParmedActions.changeprotstate(amber_prmtop, var1, var2)
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
   messages.add_line(str(action))
   action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def scee(amber_prmtop, messages):
   # We need a value
   widget_list = [('Entry', '1-4 Electrostatic Scaling Factor')]
   var_list = [StringVar()]
   description = 'Adjust the scaling factor for 1-4 electrostatic interactions'

   cmd_window = _guiwidgets.ActionWindow('scee', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   
   var = var_list[0].get()

   # Bail out if we didn't get any variables
   if not var: return

   try:
      action = ParmedActions.scee(amber_prmtop, var)
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
   messages.add_line(str(action))
   action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def scnb(amber_prmtop, messages):
   # We need a value
   widget_list = [('Entry', '1-4 van der Waals Scaling Factor')]
   var_list = [StringVar()]
   description = 'Adjust the scaling factor for 1-4 van der Waals interactions'

   cmd_window = _guiwidgets.ActionWindow('scee', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   
   var = var_list[0].get()

   # Bail out if we didn't get any variables
   if not var: return

   try:
      action = ParmedActions.scnb(amber_prmtop, var)
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
   messages.add_line(str(action))
   action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printljtypes(amber_prmtop, messages):
   "Prints all of the atoms that have the same LJ type as atoms in a given mask"
   # We need a mask
   widget_list = [('MaskEntry', 'Atom Mask or LJ Type Index')]
   var_list = [StringVar()]
   description = ('Get a list of all atoms that share a common LJ type')
   cmd_window = _guiwidgets.ActionWindow('printLJTypes', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   
   var = var_list[0].get()

   # Bail out if we didn't get any variables
   if not var: return
   
   # Instantiate our action
   try:
      action = ParmedActions.printljtypes(amber_prmtop, var)
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return

   # Create our Info window
   window = Toplevel()
   window.resizable(width=False, height=False)
   window.title('LJ Type List')
   text = _guiwidgets.ExitingScrollText(window, None, spacing3=2, padx=5,
                                        pady=5, width=45, height=30)
   text.pack()
   text.grid(row=0, column=0)
   text.text.insert(END, str(action))
   window.wait_window()
      
   messages.add_line('Printed LJ types for [%s]' % var)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def changeljsingletype(amber_prmtop, messages):
   """ Changes radius/well depth of a single LJ type given by the mask """
   # We need a mask, radius, and well depth
   widget_list = [('MaskEntry', 'Mask to change LJ Type'),
                  ('Entry', 'New LJ Radius'), ('Entry', 'New LJ Depth')]
   var_list = [StringVar(), StringVar(), StringVar()]
   description = "Change a given atom type's LJ radius and well depth"
   # Create the window, open it, then wait for it to close
   cmd_window = _guiwidgets.ActionWindow('changeLJSingleType', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   vars_found = True in [bool(v.get()) for v in var_list]
   if not vars_found: return
   # addljtype expects any _non_specified variables to be None
   var_list = [v.get() for v in var_list]
   for i, v in enumerate(var_list):
      if not v: var_list[i] = None
   # If we did, store them and pass it to the class
   var1, var2, var3, = var_list[0], var_list[1], var_list[2]
   try:
      action = ParmedActions.changeljsingletype(amber_prmtop,var1,var2,var3)
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
   messages.add_line(str(action))
   action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def combinemolecules(amber_prmtop, messages):
   """ Combines 2 molecules into a single molecule """
   # We need a molecule #
   widget_list = [('Entry', 'Molecule Number')]
   var_list = [StringVar()]
   description = 'Combine the given molecule number with the next molecule'
   cmd_window = _guiwidgets.ActionWindow('combineMolecules', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   var = var_list[0].get()
   if not var: return
   # addljtype expects any _non_specified variables to be None
   try:
      action = ParmedActions.combinemolecules(amber_prmtop, var)
      messages.add_line(str(action))
      action.execute()
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def addcoarsegrain(amber_prmtop, messages):
   """ Adds coarse graining to topology file via Lula's algo """
   # This implementation doesn't exist anywhere yet, so disable it
   showerror('Warning', 'This functionality is not implemented in Amber yet!')
   return

   # We need a file
   fname = file_chooser('Coarse Grain Parameter', 
                        [('Coarse Grain Parameters', '*.cgparm'),
                         ('All Files', '*')])
   
   if not fname: return

   try:
      action = ParmedActions.addcoarsegrain(amber_prmtop, fname)
      messages.add_line(str(action))
      action.execute()
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
   
#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def definesolvent(amber_prmtop, messages):
   """ Allows you to define what you consider to be solvent molecules """
   # We need a molecule #
   widget_list = [('Entry', 'List of solvent residue names')]
   var_list = [StringVar()]
   description =('Tell ParmEd that these residue names should be considered\n' +
         'to be solvent residues. This should be a comma-delimited list')
   cmd_window = _guiwidgets.ActionWindow('defineSolvent', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   var = var_list[0].get()
   if not var: return
   # addljtype expects any _non_specified variables to be None
   try:
      action = ParmedActions.definesolvent(amber_prmtop, var)
      messages.add_line(str(action))
      action.execute()
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printbonds(amber_prmtop, messages):
   """ Prints bonds containing atoms from a given mask """
   # We need a mask, new radius, new epsilon, and for chamber topologies,
   # a new radius-1-4 and epsilon-1-4. Don't add the latter ones until we
   # know if we have a chamber prmtop or not.
   widget_list = [('MaskEntry', 'Atoms to analyze for bonds')]
   # We need 5 string variables, then get the description. 
   var_list = [StringVar()]
   description='Prints all bonds containing at least 1 atom in the given mask'
   # Create the window, open it, then wait for it to close
   cmd_window = _guiwidgets.ActionWindow('printBonds', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   var = var_list[0]
   if not var.get(): return
   try:
      action = ParmedActions.printbonds(amber_prmtop, var.get())
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
   messages.add_line('Printed BONDs for %s' % var.get())
   # Now make the text window
   window = Toplevel()
   window.resizable(width=False, height=False)
   window.title('BOND list in %s' % amber_prmtop)
   text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                        pady=5, width=65, height=20)
   text.pack()
   text.grid(row=0, column=0)
   text.text.insert(END, str(action))
   window.wait_window()
#  window.grab_release()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printangles(amber_prmtop, messages):
   """ Prints angles containing atoms from a given mask """
   # We need a mask, new radius, new epsilon, and for chamber topologies,
   # a new radius-1-4 and epsilon-1-4. Don't add the latter ones until we
   # know if we have a chamber prmtop or not.
   widget_list = [('MaskEntry', 'Atoms to analyze for angles')]
   # We need 5 string variables, then get the description. 
   var_list = [StringVar()]
   description='Prints all angles containing at least 1 atom in the given mask'
   # Create the window, open it, then wait for it to close
   cmd_window = _guiwidgets.ActionWindow('printAngles', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   var = var_list[0]
   if not var.get(): return
   try:
      action = ParmedActions.printangles(amber_prmtop, var.get())
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
   messages.add_line('Printed ANGLEs for %s' % var.get())
   # Now make the text window
   window = Toplevel()
   window.resizable(width=False, height=False)
   window.title('ANGLE list in %s' % amber_prmtop)
   text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                        pady=5, width=90, height=20)
   text.pack()
   text.grid(row=0, column=0)
   text.text.insert(END, str(action))
   window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printdihedrals(amber_prmtop, messages):
   """ Prints dihedrals containing atoms from a given mask """
   # We need a mask, new radius, new epsilon, and for chamber topologies,
   # a new radius-1-4 and epsilon-1-4. Don't add the latter ones until we
   # know if we have a chamber prmtop or not.
   widget_list = [('MaskEntry', 'Atoms to analyze for dihedrals')]
   # We need 5 string variables, then get the description. 
   var_list = [StringVar()]
   description=\
      'Prints all dihedrals containing at least 1 atom in the given mask'
   # Create the window, open it, then wait for it to close
   cmd_window = _guiwidgets.ActionWindow('printDihedrals', amber_prmtop,
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   var = var_list[0]
   if not var.get(): return
   try:
      action = ParmedActions.printdihedrals(amber_prmtop, var.get())
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
   messages.add_line('Printed DIHEDRALs for %s' % var.get())
   # Now make the text window
   window = Toplevel()
   window.resizable(width=False, height=False)
   window.title('DIHEDRAL list in %s' % amber_prmtop)
   text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                        pady=5, width=140, height=20)
   text.pack()
   text.grid(row=0, column=0)
   text.text.insert(END, str(action))
   window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def setbond(amber_prmtop, messages):
   """ Sets (adds or changes) a bond in the topology file """
   # We need 2 masks, a force constant, and an equilibrium distance
   widget_list = [('MaskEntry', 'First atom in bond'),
                  ('MaskEntry', 'Second atom in bond'),
                  ('Entry', 'Force constant (kcal/mol Ang**2)'),
                  ('Entry', 'Equilibrium Distance (Ang)')]
   # We need 4 variables
   var_list = [StringVar(), StringVar(), StringVar(), StringVar()]
   description = ('Sets a bond in the topology file with the given Force ' +
                  'constant in kcal/mol/Ang**2\nand the given equilibrium ' +
                  'distance in Angstroms. Both masks must specify only a\n' +
                  'single atom. If the bond exists, it will be replaced. If ' +
                  'it doesn\'t, it will be added.')
   # Create the window, open it, then wait for it to close
   cmd_window = _guiwidgets.ActionWindow('setBond', amber_prmtop, 
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   vars_found = True in [bool(v.get()) for v in var_list]
   if not vars_found: return
   # If we did, pass them through
   var_list = [v.get() for v in var_list]
   var1, var2, var3, var4 = var_list[0], var_list[1], var_list[2], var_list[3]
   try:
      action = ParmedActions.setbond(amber_prmtop, var1, var2, var3, var4)
      messages.add_line(str(action))
      action.execute()
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def setangle(amber_prmtop, messages):
   """ Sets (adds or changes) an angle in the topology file """
   # We need 3 masks, a force constant, and an equilibrium angle
   widget_list = [('MaskEntry', 'First atom in angle'),
                  ('MaskEntry', 'Second (middle) atom in angle'),
                  ('MaskEntry', 'Third atom in angle'),
                  ('Entry', 'Force constant (kcal/mol rad**2)'),
                  ('Entry', 'Equilibrium Angle (Degrees)')]
   # We need 5 variables
   var_list = [StringVar(), StringVar(), StringVar(), StringVar(), StringVar()]
   description = ('Sets an angle in the topology file with the given Force ' +
                  'constant in kcal/mol/rad**2\nand the given equilibrium ' +
                  'angle in Degrees. All three masks must specify only a\n' +
                  'single atom. If the angle exists, it will be replaced. If ' +
                  'it doesn\'t, it will be added.')
   # Create the window, open it, then wait for it to close
   cmd_window = _guiwidgets.ActionWindow('setAngle', amber_prmtop, 
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   vars_found = True in [bool(v.get()) for v in var_list]
   if not vars_found: return
   # If we did, pass them through
   var_list = [v.get() for v in var_list]
   var1,var2,var3,var4,var5 = var_list[0],var_list[1],var_list[2],var_list[3],\
                              var_list[4]
   try:
      action = ParmedActions.setangle(amber_prmtop, var1,var2,var3,var4,var5)
      messages.add_line(str(action))
      action.execute()
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def adddihedral(amber_prmtop, messages):
   """ Adds a dihedral (improper, multiterm, or normal) to the prmtop """
   # We need 4 masks, phi_k, periodicity, phase, scee/scnb, dihedral type
   widget_list = [('MaskEntry', 'First (end) atom in dihedral'),
                  ('MaskEntry', 'Second (middle) atom in dihedral'),
                  ('MaskEntry', 'Third (middle) atom in dihedral'),
                  ('MaskEntry', 'Fourth (end) atom in dihedral'),
                  ('Entry', 'Phi Force constant (kcal/mol)'),
                  ('Entry', 'Periodicity'),
                  ('Entry', 'Phase (Degrees)'),
                  ('Entry', 'EEL scaling factor'),
                  ('Entry', 'VDW scaling factor'),
                  ('Spinbox', 'Dihedral type', 'normal','multiterm','improper')]
   # We need 10 variables
   var_list = [StringVar() for i in range(10)]
   description = ('Adds a dihedral in the topology file with the given Phi ' +
                  'Force constant in kcal/mol the\ngiven phase in Degrees ' +
                  'and the given periodicity. All masks must specify only \n' +
                  'a single atom. The default Amber values for SCEE/SCNB are ' +
                  '1.2 and 2.0, respectively.\nSee the format documentation ' +
                  'for details about normal, multiterm, and improper dihedrals')
   # Create the window, open it, then wait for it to close
   cmd_window = _guiwidgets.ActionWindow('setDihedral', amber_prmtop, 
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   vars_found = True in [bool(v.get()) for v in var_list]
   if not vars_found: return
   # If we did, pass them through
   var_list = [v.get() for v in var_list]
   var1, var2, var3, var4 = var_list[0], var_list[1], var_list[2], var_list[3]
   var5, var6, var7, var8 = var_list[4], var_list[5], var_list[6], var_list[7]
   var9, var10 = var_list[8], var_list[9]
   # Fill scee/scnb in with default values
   if not var8: var8 = '1.2' # scee
   if not var9: var9 = '2.0' # scnb
   try:
      action = ParmedActions.adddihedral(amber_prmtop, var1, var2, var3, var4,
                                         var5, var6, var7, var8, var9, var10)
      messages.add_line(str(action))
      action.execute()
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def addatomicnumber(amber_prmtop, messages):
   """ Asks the user if they want to add ATOMIC_NUMBER to the prmtop """
   response = askyesno('addAtomicNumber',
                    'Do you want to add the ATOMIC_NUMBER section to %s?' % 
                    amber_prmtop)
   if response:
      action = ParmedActions.addatomicnumber(amber_prmtop)
      action.execute()
      messages.add_line(str(action))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def deletedihedral(amber_prmtop, messages):
   """ Deletes a dihedral between 4 given atoms """
   # We need 4 masks, phi_k, periodicity, phase, scee/scnb, dihedral type
   widget_list = [('MaskEntry', 'First (end) atom in dihedral'),
                  ('MaskEntry', 'Second (middle) atom in dihedral'),
                  ('MaskEntry', 'Third (middle) atom in dihedral'),
                  ('MaskEntry', 'Fourth (end) atom in dihedral')]
   # We need 10 variables
   var_list = [StringVar() for i in range(4)]
   description = ('Deletes dihedrals between the atoms specified in mask1, ' +
                  'mask2, mask3, and mask4.\nIt will try to match dihedrals ' +
                  'only in the order given or reverse order. Dihedrals are\n' +
                  'specified by atom N in mask1, mask2, mask3, and mask4,' +
                  'where N is the selected\natom in that mask.')
   # Create the window, open it, then wait for it to close
   cmd_window = _guiwidgets.ActionWindow('deleteDihedral', amber_prmtop, 
                     widget_list, var_list, description)
   cmd_window.wait_window()
   # See if we got any variables back
   vars_found = True in [bool(v.get()) for v in var_list]
   if not vars_found: return
   # If we did, pass them through
   var_list = [v.get() for v in var_list]
   var1, var2, var3, var4 = var_list[0], var_list[1], var_list[2], var_list[3]
   try:
      action = ParmedActions.deletedihedral(amber_prmtop, var1,var2,var3,var4)
      messages.add_line(str(action))
      action.execute()
   except Exception, err:
      showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
      return
