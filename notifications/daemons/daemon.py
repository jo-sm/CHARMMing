"""
CHARMMing Websockets Notifications Daemon
Joshua Smock

Daemon that handles notifications to/from ZeroMQ and the browser. 
Is only a client and listens on port 55140, publishes on port
55150. Need to extend to support multiple services, as this only
serves a single service (the scheduler)
"""
import sockjs.tornado
import zmq
import logging
from zmq.eventloop import ioloop
from zmq.eventloop.zmqstream import ZMQStream
ioloop.install()

# Handles the websocket connections
class NotificationsWebsocket(sockjs.tornado.SockJSConnection):
  clients = set()

  context = zmq.Context()

  # Set up the publisher socket
  publisher_socket = context.socket(zmq.PUB)
  publisher_socket.bind('tcp://127.0.0.1:55150')

  # Set up the subscriber socket
  subscriber_socket = context.socket(zmq.SUB)
  subscriber_socket.connect('tcp://127.0.0.1:55140')
  subscriber_socket.setsockopt(zmq.SUBSCRIBE, '') # No channel prefix

  def on_open(self, request):
    print('Client connected: ', request.cookies)

  def on_message(self, message):
    print('Message: ', message)
    self.publish_stream.send_unicode(message)

  def on_close(self):
    print('Client disconnected')
    self.clients.remove(self)
