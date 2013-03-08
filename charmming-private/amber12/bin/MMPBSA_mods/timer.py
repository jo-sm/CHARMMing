"""
This is a module that contains the timer facility for MMPBSA.py. It's 
useful for profiling the performance of the script.                   

Last updated: 04/17/2011                                     

                           GPL LICENSE INFO                             

Copyright (C) 2009 - 2011 Dwight McGee, Billy Miller III, and Jason Swails

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""

from time import time
from MMPBSA_mods.exceptions import TimerException

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Timer(object):
   """ Timer class. It adds new timers and keeps track of how much time has been
       spent. """

   # Do we force the issue? i.e. : If we want to start a timer that doesn't 
   # exist, add one and start it. Allow turning a non-off timer off, or 
   # turning one on that is already on
   force = False

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __init__(self):
      """ Declares the global timer """
      self.timers = {'global' : -time()}
      self.descriptions = {'global' : 'Total time taken:'}
      self.active_timers = ['global']
      self.timer_names = ['global']
      self.units = 'sec.'

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def AddTimer(self, timer_name, description):
      """ Add a new timer """
      # Check to make sure it doesn't already exist. Then add it if it doesn't
      if timer_name in self.timer_names:
         if not self.force:
            raise(TimerException("%s is already a timer!" % timer_name))
         else:
            return 0

      self.timers[timer_name] = 0
      self.timer_names.append(timer_name)
      self.descriptions[timer_name] = description

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def StartTimer(self, timer_name):
      """ Start the specified timer """
      # Make sure the timer is not on already
      if timer_name in self.active_timers and not self.force:
         raise(TimerException("%s is already active!" % timer_name))
      
      # Check to make sure we've added the timer or not. If not, then add it
      if not timer_name in self.timer_names:
         if not self.force:
            raise TimerException("Failed to start non-existent timer %s" % 
                                 timer_name)
         else:
            self.AddTimer(timer_name, '')

      self.timers[timer_name] -= time()

      # This timer is now on
      self.active_timers.append(timer_name)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def StopTimer(self, timer_name):
      """ End the specified timer """
      # First make sure that it is a timer to begin with, or raise an exception
      if not timer_name in self.timer_names:
         raise(TimerException("%s is not a timer!" % timer_name))

      # Next make sure that the timer is on in the first place
      if not timer_name in self.active_timers:
         if not self.force:
            raise(TimerException("%s timer is not on!" % timer_name))
         else:
            return 0

      # Now if it's on, end it and remove it from the list of active timers
      self.timers[timer_name] += time()
      self.active_timers.pop(self.active_timers.index(timer_name))

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def EndAll(self):
      """ End all of the timers """
      while len(self.active_timers) > 0:
         self.StopTimer(self.active_timers[0])

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Done(self):
      """ Tell Timer we are done, and we can convert units to best 
          human-readable option 
      """

      # Make sure all timers are ended
      self.EndAll()
      
      tfactor = 1
      # Now test the magnitude of the global timer so we can decide what the 
      # reported units should be
      if self.timers['global'] > 60 * 60 * 48:
         self.units = 'days'
         tfactor = 60 * 60 * 24
      elif self.timers['global'] > 60 * 60:
         self.units = 'hr.'
         tfactor = 60 * 60
      elif self.timers['global'] > 60:
         self.units = 'min.'
         tfactor = 60

      for timer in self.timer_names:
         self.timers[timer] /= tfactor

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Print(self, timer):
      """ Prints the value of the timer """
      return "%-40s %8.3f %s" % (self.descriptions[timer], self.timers[timer], self.units)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
