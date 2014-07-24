"""
CHARMMing Scheduler Daemon
Joshua Smock

A dispatcher for various schedulers supported for job queueing in CHARMMing. 
Listens on port 9995 as a ZMQ PULL socket
"""

import os
import sys
import socket
import threading
import time
import zmq

class daemon():
  # Type is the kind of supported scheduler (pbs, simple) and the port it resides on
  def __init__(self, type, address, port):
    module = __import__(type)
    queue_ = getattr(module, 'Queue')
    self.queue = queue_(address, port)

  def run(self):
    # We query jobs every 0.5 seconds
    # If a job status changes, we update the database and notify the notification system to send out notifications
    while True:
      # Listen on the zeromq socket for a job submission
      print "Querying jobs"
      time.sleep(0.5)

if __name__ == "__main__":
  if len(sys.argv) < 3:
    sys.exit('Error: Scheduler daemon requires a scheduler class (e.g. pbs, simple), an address, and a port.')
  _daemon = daemon(sys.argv[0], sys.argv[1], sys.argv[2])
  _daemon.run()
