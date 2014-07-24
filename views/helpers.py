import string
import random
from charmming.models import User

def id_generator(size=8, chars=string.ascii_uppercase + string.digits):
  return ''.join(random.choice(chars) for _ in range(size))

"""
There are three problems with the following method:
1) It does the calling on objects, i.e. I call User.most_recent_sessions(n) to get the n most recent sessions
2) This performs raw SQL
3) It is a function rather than a method on the User class

Problem 1) is not a big issue, but will need to be reevaluated at a later time to determine if this is the best
way to accomplish this. Problem 2) cannot be dealt with in Django's ORM, as far as I know, but because
it does not perform any updates and is ANSI SQL compliant, it shouldn't pose any issues between databases.
Problem 3) will be solved in a later update.

See http://www.dabapps.com/blog/higher-level-query-api-django-orm/
"""
def most_recent_sessions(n):
  users = User.objects.raw("select count(user.id) as count, user.*, session.last_accessed from ( select user_id, max(last_accessed) as last_accessed from charmming_session group by user_id) session join charmming_user user on user.id = session.user_id order by session.last_accessed desc limit %s;", [n])
  return list(users)
