"""
The GUI components of xparmed.py
"""
import tkFileDialog, tkMessageBox
from Tkinter import *
from ParmedTools.gui.guiactions import gui_action_dispatcher
from ParmedTools.gui._guiwidgets import MessageWindow

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def action_list(root, amber_prmtop):
   """ Generates the action list on the given window """
   import math
   from ParmedTools import ParmedActions

   # Make a dict with each action pairing to their description so we can easily
   # display them
   actions = {}
   for item in ParmedActions.usages.keys():
      action_name = ParmedActions.usages[item].split()[0]
      # Skip some actions, because they only make sense in the context of
      # the command-line parmed. SetOverwrite is checked on saving files, and
      # is always set to True for xparmed.py, so skip that action here
      if item in ('parmout', 'go', 'quit', 'help', 'setoverwrite'): continue
      # Save outparm for later
      elif item == 'outparm':
         outparm_action = _Description(action_name,
                              getattr(ParmedActions, item).__doc__)
         continue
      actions[item] = _Description(action_name, 
                                   getattr(ParmedActions, item).__doc__)

   # We will pack each item next to a button with a short description. This
   # description button can be pressed to get a full description. We want 10
   # items per set of 2 columns, but a header for each one.
   ITEMS_PER_COLUMN = int(math.ceil((len(actions)) / 2.0))

   sorted_actions = actions.keys()
   sorted_actions.sort()
   root_frames = []
   current_frame = None
   column_number = 0
   message_frame = Frame(root, padx=5, pady=10)
   messages = MessageWindow(message_frame)
   for i, action in enumerate(sorted_actions):
      if i % ITEMS_PER_COLUMN == 0:
         if current_frame:
            current_frame.pack()
            current_frame.grid(column=column_number, row=0, sticky=N+S+E+W)
            column_number += 1
         current_frame = Frame(root, padx=5, pady=10)
         root_frames.append(current_frame)
         my_label = Label(current_frame, text='Action Name')
         my_label.pack()
         my_label.grid(column=0, row=0, sticky=N+S+E+W)
         my_label = Label(current_frame, text='Description')
         my_label.pack()
         my_label.grid(column=1, row=0, sticky=N+S+E+W)

      my_button = ActionButton(current_frame, amber_prmtop, messages, 
                               actions[action])
      my_button.pack()
      my_button.grid(column=0, row=i%ITEMS_PER_COLUMN+1, sticky=N+S+E+W)
      my_button = HelpButton(current_frame, actions[action])
      my_button.pack()
      my_button.grid(column=1, row=i%ITEMS_PER_COLUMN+1, sticky=N+S+E+W)
      i += 1
   
   # Print out the last column
   current_frame.pack()
   current_frame.grid(column=column_number, row=0, sticky=N+S+E+W)
   column_number += 1
   message_frame.pack()
   message_frame.grid(column=column_number, row=0, sticky=N+S+E+W)
   column_number += 1
   # Now print out the exit button
   end1 = int(math.ceil(column_number / 2.0))
   span2 = column_number - end1
   my_button = OutparmButton(root, amber_prmtop, messages, outparm_action)
   my_button.pack()
   my_button.grid(column=0, row=1, sticky=N+S+E+W, columnspan=end1)
   my_button = Button(root, text='Exit xParmEd', command=root.quit)
   my_button.pack()
   my_button.grid(column=end1, row=1, sticky=N+S+E+W, columnspan=column_number)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class _Description(object):
   """ This formats a description so it can be displayed in different places """
   #====================================

   def __init__(self, action_name, desc_string):
      self.desc_string = desc_string
      self.action_name = action_name

   #====================================

   def button_format(self):
      """ 
      Returns a string formatted to fit inside a button. 2 lines of 20
      characters each. These can be changed
      """
      MAX_LENGTH = 20
      MAX_LINES = 1
      words = self.desc_string.split()
      desc = ''
      wordnum = 0
      done = False
      for i in range(MAX_LINES):
         line = ''
         # For the last line, allow room for an elipsis (...)
         if i == MAX_LINES - 1: line_len = MAX_LENGTH - 3
         else: line_len = MAX_LENGTH
         while len(line) < line_len:
            wordnum += 1
            if len(line) + len(words[wordnum-1]) > line_len:
               wordnum -= 1
               break
            elif len(line) + len(words[wordnum-1]) == line_len:
               line += words[wordnum-1]
               if wordnum == len(words): done = True
               break
            else: # still less than line_len
               line += words[wordnum-1] + ' '
               if wordnum == len(words):
                  done = True
                  break
         # end while
         desc += '\n' + line
         if done: break
      # end for

      # If we didn't hit the end of the description, we need to add an ellipsis
      if not done: desc += '...'
      # Now strip off the leading/trailing whitespace and return the result
      desc = desc.strip()
      return desc

   #====================================

   def info_box(self):
      """ Display an info button """
      tkMessageBox.showinfo(title=self.action_name, 
               message=str(' '.join(self.desc_string.split())).strip())

   #====================================

def test_description():
   """ A method to test the implementation of the _Description class """
   my_desc1 = _Description('changeLJPair',
                           'Changes a particular Lennard Jones pair ' +
                           'based on a given (pre-combined) epsilon/Rmin')
   my_desc2 = _Description('scee', 'Change the scee scaling factor')

   check1 = my_desc1.button_format()
   check2 = my_desc2.button_format()
   correct1 = 'Changes a particular\nLennard Jones ...'
   correct2 = 'Change the scee \nscaling factor'

   print 'Testing button_format:\nmy_desc1:\n'
   print check1
   print '\nTesting button_format:\nmy_desc2:\n'
   print check2

   print ''

   if check1 == correct1: print 'TEST 1 PASSED'
   else: print 'TEST 1 FAILED'
   if check2 == correct2: print 'TEST 2 PASSED'
   else: print 'TEST 2 FAILED'

   print '\nTesting info boxes. They should display the full description in an'
   print 'info box. Check that they pass.'

   check1 = my_desc1.info_box()
   check2 = my_desc2.info_box()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class HelpButton(Button):
   """ A button to display a help """
   def __init__(self, root, action_desc):
      self.action_desc = action_desc
      Button.__init__(self, root, text=action_desc.button_format(),
                      command=self.execute)
   def execute(self):
      self.action_desc.info_box()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class ActionButton(Button):
   " A button to dispatch the command to a particular function in guiactions "
   def __init__(self, root, amber_prmtop, messages, action_desc):
      self.amber_prmtop = amber_prmtop
      self.action_desc = action_desc
      self.messages = messages
      Button.__init__(self, root, text=action_desc.action_name,
                      command=self.execute)
   def execute(self):
      """ Sends the action name to the GUI dispatcher """
      gui_action_dispatcher(self.amber_prmtop, self.action_desc.action_name,
                            self.messages)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class OutparmButton(ActionButton):
   " An ActionButton specifically for outparm "
   def __init__(self, root, amber_prmtop, messages, action_desc):
      self.amber_prmtop = amber_prmtop
      self.action_desc = action_desc
      self.messages = messages
      Button.__init__(self, root, text='Write a New Topology File Now ' +
                      '(can be used any number of times)', command=self.execute)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

if __name__ == '__main__':
   # If this is the main script, run all of the tests
   test_description()
