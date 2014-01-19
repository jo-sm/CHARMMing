#
#                            PUBLIC DOMAIN NOTICE
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the authors' official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software is freely available
#  to the public for use.  There is no restriction on its use or
#  reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, NIH and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. NIH, NHLBI, and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any
#  particular purpose.

import pychm, re
import charmming_config
import os, json, sha, tempfile, cPickle
from httplib import HttpConnection
from api.models import APIUser, APIJob, APIOptimization, APIEnergy
from django.http import HttpResponseRedirect, HttpResponse

def validateUser(postdata):

    if not postdata.has_key('data'):
        return -1, None

    ddct = json.loads(postdata['data'])
    if not ddct.has_key('user_id'):
        return -2, None
    try:
        client = APIUser.objects.get(callerName=ddct['user_id'])
    except:
        return -2, None

    if not ddct.has_key('auth'):
        return -2, None

    # verify the authentication hash sent by the client against their API key
    # allow for some skew between client and server clocks and transmit time
    okc = False
    now = int(time.time())
    for i in range(now-10,now+10): 
        hstr = client.callerKey + str(i)
        hash = sha.new(hstr)
        if hash.hexdgest() == ddct['auth']:
            okc = True
            break
    if not okc:
        return -2, None

    return 0, ddct
        

def setUPDirectory(client,req):
    dirname = tempfile.mkdtemp(prefix='%s/%s/job' % (charmming_config['user_home'],client['user_id']))

    if req.has_key('start_struct'):
        fp = open('%s/start.pdb' % dirname, 'w')
        fp.write(req['start_struct'])
        fp.close()
    elif req.has_key('start_struct_url'):
        m = re.search('^http://([^/]+)(.*?)$',req['start_struct_url'])
        if not m:
            return None

        conn = HTTPConnection(m.group(0))
        conn.request("GET", m.group(1))
        resp = conn.getresponse()
        fp = open('%s/start.pdb' % dirname, 'w')
        if resp.status != 200:
            return None
        pdb_file = resp.read()

        fp.write(pdb_file)
        fp.close()
    else:
        return None

    try:
        mol = pychm.io.pdb.get_molFromPDB('%s/start.pdb' % self.directory)
        for segment in mol.iter_segs():
            pass
    except:
        return None

    fp = open('%s/pickle.dat' % dirname,'w')
    cPickle.dump(mol,fp)
    fp.close()

def energyCall(request):
    response = {}

    err, rcall = validateUser(request.POST)
    if err < 0:
        response['errcode'] = err
        s = json.dumps(response)
        return HttpResponse(s)
    
    client = APIUser.objects.get(callerName=rcall['user_id'])
    req = json.loads(rcall)
    dirname = setUPDirectory(client,req)
    if not dirname:
        response['errcode'] = -100
        s = json.dumps(response)
        return HttpResponse(s)

    myjob = APIEnergy()
    myjob.jobType = 1
    myjob.directory = dirname
    myjob.save()
    response = myjob.run(req)

    s = json.dumps(response)
    return HttpResponse(s)

def minimizeCall(request):
    response = {}

    err, rcall = validateUser(request.POST)
    if err < 0:
        response['errcode'] = err
        s = json.dumps(response)
        return HttpResponse(s)

    client = APIUser.objects.get(callerName=rcall['user_id'])
    req = json.loads(rcall)

    s = json.dumps(response)
    return HttpResponse(s)
