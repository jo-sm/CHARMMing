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
import os, json, sha, tempfile, cPickle, time
import traceback, socket
import output
import openbabel
from httplib import HTTPConnection
from api.models import APIUser, APIJob, APIOptimization, APIEnergy
from django.http import HttpResponseRedirect, HttpResponse
from django.template.loader import get_template
from django.template import Context

def validateUser(data):

    logfp = open('/tmp/validateuser.log','w')
    logfp.write('rewind the octopus!\n')
    logfp.flush()

    logfp.write('unserialized OK\n')
    if not data.has_key('user_id'):
        logfp.write('data does not have user_id key\n')
        logfp.close()
        return -2, None
    try:
        client = APIUser.objects.get(callerName=data['user_id'])
    except:
        logfp.write('could not look up user %s\n' % data['user_id'])
        logfp.close()
        return -2, None

    if not data.has_key('auth'):
        logfp.write('data has no auth key\n')
        logfp.close()
        return -2, None

    # verify the authentication hash sent by the client against their API key
    # allow for some skew between client and server clocks and transmit time
    okc = False
    now = int(time.time())
    for i in range(now-10,now+10): 
        hstr = client.callerKey + str(i)
        hash = sha.new(hstr)
        logfp.write('try has string: %s hash = %s\n' % (hstr,hash.hexdigest()))
        if hash.hexdigest() == data['auth']:
            okc = True
            break
    if not okc:
        logfp.write('bad times in hash land.\n')
        logfp.close()
        return -2, None

    logfp.write('all square, baby!\n')
    logfp.close()
    return 0, data
        
def badRes(seg,dname):

    logfp = open('/tmp/api_badres.txt','w')
    
    streamList = []
    badResList = []
    for res in seg.iter_res():
        logfp.write('dealing with res %s\n' % res.resName)
        if res.resName in badResList:
            continue

        filename_noh  = '%s/%s-%s.pdb' % (dname,seg.segid,res.resName)
        filename_mol2 = '%s/%s-%s.mol2' % (dname,seg.segid,res.resName)

        res.write(filename_noh,outformat='pdborg')
        logfp.write('wrote out pdb\n')

        obconv = openbabel.OBConversion()
        mol = openbabel.OBMol()
        obconv.SetInAndOutFormats("pdb","mol2")
        obconv.ReadFile(mol, filename_noh.encode("utf-8"))
        mol.AddHydrogens()
        mol.SetTitle(res.resName)
        obconv.WriteFile(mol, filename_mol2.encode("utf-8"))
        logfp.write('converted it to MOL2\n')

        header = '# USER_IP 165.112.184.30 USER_LOGIN apiuser\n\n'
        fp = open(filename_mol2,'r')
        payload = fp.read()
        fp.close()
        content = header + payload

        recvbuff = ''
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.connect((charmming_config.cgenff_host, charmming_config.cgenff_port))
            s.sendall(content)
            s.shutdown(socket.SHUT_WR)
            while True:
                data = s.recv(1024)
                if not data:
                    break
                recvbuff += data

        except Exception, e:
            return None

        logfp.write('got response from paramChem\n')
        fp = open('%s/stream-%s-%s.str' % (dname,seg.segid,res.resName), 'w')
        instr = False
        for line in recvbuff.split('\n'):
            if not instr:
                if line.startswith('*'):
                    instr = True
                else:
                    continue
            fp.write(line + '\n')
        fp.close()
        streamList.append('%s/stream-%s-%s.str' % (dname,seg.segid,res.resName))

    logfp.write('All done!\n')
    logfp.close()
        
    return streamList

def setUPDirectory(client,req,job):
    dirname = tempfile.mkdtemp(prefix='%s/%s/job' % (charmming_config.user_home,client.callerName))
    os.chmod(dirname, 0775)

    if req.has_key('start_struct'):
        fp = open('%s/start.pdb' % dirname, 'w')
        fp.write(req['start_struct'])
        fp.close()
    elif req.has_key('start_struct_url'):
        m = re.search('^http://([^/]+)(.*?)$',req['start_struct_url'])
        if not m:
            return None

        conn = HTTPConnection(m.group(1))
        conn.request("GET", m.group(2))
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
        pdb = pychm.io.pdb.PDBFile('%s/start.pdb' % dirname)
        mol = pdb[0]
        mol.parse()
    except:
        logfp = open('/tmp/setupdir2.log','w')
        logfp.write('It is an abomination!\n')
        logfp.flush()
        traceback.print_exc(file=logfp)
        logfp.close()

        return None

    segList = []
    rtfList = set()
    prmList = set()
    strList = set()
    os.chdir(dirname)
    mainlog = open('/tmp/setupdir.log','w')
    for seg in mol.iter_seg():
        segList.append(seg.segid)

        seg.write('%s/segment-%s.pdb' % (dirname,seg.segid), outformat='charmm')

        # build it
        myRtfs = []
        myPrms = []
        myStrs = []
        if seg.segType == 'good':
            myRtfs.append('%s/toppar/%s' % (charmming_config.data_home,charmming_config.default_pro_top))
            myPrms.append('%s/toppar/%s' % (charmming_config.data_home,charmming_config.default_pro_prm))
            myStrs.append('%s/toppar/toppar_water_ions.str' % charmming_config.data_home)
            genopts = 'first none last none noangle nodihedral'
            doic = False
        elif seg.segType == 'pro':
            myRtfs.append('%s/toppar/%s' % (charmming_config.data_home,charmming_config.default_pro_top))
            myPrms.append('%s/toppar/%s' % (charmming_config.data_home,charmming_config.default_pro_prm))
            genopts = ''
            doic = True
        elif seg.segType == 'bad':
            badtp = badRes(seg,dirname)
            if not badtp:
                return None
            myRtfs.append('%s/toppar/top_all36_cgenff.rtf' % charmming_config.data_home)
            myPrms.append('%s/toppar/par_all36_cgenff.prm' % charmming_config.data_home)
            for item in badtp:
                myStrs.append(item)
            genopts = 'first none last none'
            doic = False
        else:
            # let's not deal with anything else yet...
            logfp = open('/tmp/setupdir3.log','w')
            logfp.write('Cannot deal with segtype %s!\n' % seg.segType)
            logfp.close()
            return None

        for x in myRtfs:
            rtfList.add(x)
        for x in myPrms:
            prmList.add(x)
        for x in myStrs:
            strList.add(x)

        template_dict = {'rtf_list': myRtfs, 'prm_list': myPrms, 'str_list': myStrs, \
                         'segname': seg.segid, 'genopts': genopts, 'doic': doic}
        t = get_template('%s/mytemplates/input_scripts/api_builds.inp' % charmming_config.charmming_root)
        charmm_inp = output.tidyInp(t.render(Context(template_dict)))

        inp_file = '%s/buildseg-%s.inp' % (dirname,seg.segid)
        out_file = '%s/buildseg-%s.out' % (dirname,seg.segid)
        fp = open(inp_file, 'w')
        fp.write(charmm_inp)
        fp.close()

        cmd = '%s < %s > %s' % (charmming_config.charmm_exe,inp_file,out_file)
        mainlog.write(cmd + '\n')
        os.system(cmd)

    # all segments built -- build final structure
    template_dict = {'rtf_list': rtfList, 'prm_list': prmList, 'str_list': strList, 'seg_list': segList}
    t = get_template('%s/mytemplates/input_scripts/api_append.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))

    job.rtfList = ''
    job.prmList = ''
    job.strList = ''
    job.segList = ''
    for rtf in rtfList:
        job.rtfList += '%s ' % rtf
    for prm in prmList:
        job.prmList += '%s ' % prm
    for str in strList:
        job.strList += '%s ' % str
    for seg in segList:
        job.segList += '%s ' % seg

    inp_file = '%s/buildstruct.inp' % dirname
    out_file = '%s/buildstruct.out' % dirname

    fp = open(inp_file, 'w')
    fp.write(charmm_inp)
    fp.close()

    cmd = '%s < %s > %s' % (charmming_config.charmm_exe,inp_file,out_file)
    mainlog.write(cmd + '\n')
    os.system(cmd)

    fp = open('%s/pickle.dat' % dirname,'w')
    cPickle.dump(mol,fp)
    fp.close()

    return dirname

def energyCall(request):
    response = {}
    logfp = open("/tmp/energycall.log","w")
    logfp.write("Piss on the floor!\n")
    logfp.flush()

    if request.method == 'POST':
        if not request.POST.has_key('data'):
            logfp.write('no data\n')
            logfp.close()
            response['errcode'] = -10
            s  = json.dumps(response)
            return HttpResponse(s, content_type='text/plain')

        ddct = json.loads(request.POST['data'])
        err, req = validateUser(ddct)
    elif request.method == 'GET':
        err, req = validateUser(request.GET)
    else:
        return HttpResponse('bad request type')

    if err < 0:
        logfp.write("validateUser failed...\n")

        response['errcode'] = err
        s = json.dumps(response)
        
        logfp.write('sendback: %s\n' %s)
        logfp.close()
        return HttpResponse(s, content_type='text/plain')
    
    logfp.write("Past validation...\n")
    logfp.flush()
    client = APIUser.objects.get(callerName=req['user_id'])

    myjob = APIEnergy()
    dirname = setUPDirectory(client,req,myjob)
    if not dirname:
        logfp.write("setUpDirectory failed...\n")

        response['errcode'] = -100
        s = json.dumps(response)
        logfp.write('send back %s\n' % s)
        logfp.close()
        return HttpResponse(s, content_type='text/plain')

    logfp.write("Serving the shizzas up from %s\n" % dirname)

    myjob.jobType = 1
    myjob.directory = dirname
    myjob.user = client
    if req.has_key('implicit_solvent'):
        myjob.implicitSolvent = 'gbmv'
    else:
        myjob.implicitSolvent = 'none'
    myjob.save()
    response = myjob.run(req)

    # test, I guess nothing went wrong here...
    s = json.dumps(response)
    logfp.write('sendback: %s\n' %s)
    logfp.close()
    return HttpResponse(s, content_type='text/plain')

def minimizeCall(request):
    response = {}

    if request.method == 'POST':
        if not request.POST.has_key('data'):
            logfp.write('no data\n')
            logfp.close()
            response['errcode'] = -10
            s  = json.dumps(response)
            return HttpResponse(s, content_type='text/plain')

        ddct = json.loads(request.POST['data'])
        err, req = validateUser(ddct)
    elif request.method == 'GET':
        err, req = validateUser(request.GET)
    else:
        return HttpResponse('bad request type')
    if err < 0:
        response['errcode'] = err
        s = json.dumps(response)
        return HttpResponse(s, content_type='text/plain')

    client = APIUser.objects.get(callerName=req['user_id'])
    myjob = APIOptimization()
    dirname = setUPDirectory(client,req,myjob)
    if not dirname:
        response['errcode'] = -100
        s = json.dumps(response)
        return HttpResponse(s, content_type='text/plain')

    myjob.jobType = 2
    myjob.directory = dirname
    myjob.user = client
    if req.has_key('implicit_solvent'):
        myjob.implicitSolvent = 'gbmv'
    else:
        myjob.implicitSolvent = 'none'
    myjob.save()
    response = myjob.run(req)

    s = json.dumps(response)
    return HttpResponse(s, content_type='text/plain')
