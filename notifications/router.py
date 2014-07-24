from notifications.daemons.daemon import NotificationsWebsocket
import sockjs.tornado
from sockjs.tornado.basehandler import BaseHandler

class HelloMoto(BaseHandler):
  def get(self):
    self.set_header('Content-Type', 'text/plain; charset=UTF-8')
    self.write("Hello moto\n")

def routes():
  url = '/ns'
  routes = sockjs.tornado.SockJSRouter(NotificationsWebsocket, url)._transport_urls
  for i, route in enumerate(routes):
    if route[0] == url + '/?':
      routes[i] = (routes[i][0], HelloMoto, routes[i][2])
  return routes
