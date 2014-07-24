# Very simple daemon that acts as a scheduler and takes care of sending jobs
# serially as they come in. Only supports running one program at a time

class Queue:
  def __init__(self, address, port):
    self.address = address
    self.port = port
