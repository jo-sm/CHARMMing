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
import charmming_config
import pychm, output
import os
from django.db import models
from django.template.loader import get_template
from django.template import Context
from minimization.models import minimizeTask
from structure.models import energyTask,WorkingFile

class APIUser(models.Model):
    callerName = models.CharField(max_length=40,null=False,default=None)
    callerKey = models.CharField(max_length=40,null=False,default=None)

class APIJob(models.Model):
    jobType = models.PositiveIntegerField(default=0)
    user = models.ForeignKey(APIUser)
    timestamp = models.DateField(auto_now_add=True)
    directory = models.CharField(max_length=30,null=False)
    rtfList = models.CharField(max_length=1000,null=True)
    prmList = models.CharField(max_length=1000,null=True)
    strList = models.CharField(max_length=1000,null=True)
    segList = models.CharField(max_length=1000,null=True)

    def parse_energy(self,fp):
        foundEne = False
        for line in fp:
            line = line.strip()

            if line.startswith('ENER>'):
                foundEne = True
                try:
                    eneval = float(line.split()[2])
                except:
                    return -9999999.99

        if not foundEne:
            return -9999999.99
        return eneval

class APIEnergy(APIJob):
    implicitSolvent = models.CharField(max_length=8,null=False,default=None)
    EnergyValue = models.FloatField(null=True)
    task = models.ForeignKey(energyTask,null=True)

    def run(self,req):
        response = {}

        if self.implicitSolvent == 'gbmv':
            template_dict = {'rtf_list': self.rtfList.split(), 'prm_list': self.prmList.split(), \
                             'str_list': self.strList.split(), 'gbmv': True}
            t = get_template('%s/mytemplates/input_scripts/api_energy.inp' % charmming_config.charmming_root)
            charmm_inp = output.tidyInp(t.render(Context(template_dict)))

            os.chdir(self.directory)
            inp_file = '%s/energy.inp' % self.directory
            out_file = '%s/energy.out' % self.directory

            fp = open(inp_file,'w')
            fp.write(charmm_inp)
            fp.close()

            cmd = '%s < %s > %s' % (charmming_config.charmm_exe,inp_file,out_file)
            os.system(cmd)

            eneout   =  out_file
        else:
             eneout = '%s/buildstruct.out' % self.directory

        # this is easy, just parse the energy from the file...
        try:
            fp = open(eneout, 'r')
        except:
            response['errcode'] = -999
            return response

        if not fp:
            response['errcode'] = -999
            return response

        ene = self.parse_energy(fp)
        fp.close()

        if ene < -9000000:
            response['errcode'] = -980
        else:
            response['errcode'] = 0
            response['eneval']  = ene

        return response

class APIOptimization(APIJob):
    implicitSolvent = models.CharField(max_length=8,null=False,default=None)
    nSDStep = models.PositiveIntegerField(default=0)
    nABNRStep = models.PositiveIntegerField(default=0)
    FinalEnergyValue = models.FloatField(null=True)
    task = models.ForeignKey(minimizeTask,null=True)

    def run(self,req):
        response = {}

        # NB: for now we just hard-code in number of steps
        template_dict = {'rtf_list': self.rtfList.split(), 'prm_list': self.prmList.split(), \
                         'str_list': self.strList.split(), 'gbmv': (self.implicitSolvent == 'gbmv'), \
                         'sdsteps': 100, 'abnrsteps': 500}
        t = get_template('%s/mytemplates/input_scripts/api_minimization.inp' % charmming_config.charmming_root)
        charmm_inp = output.tidyInp(t.render(Context(template_dict)))

        inp_file = '%s/minimize.inp' % self.directory
        out_file = '%s/minimize.out' % self.directory

        fp = open(inp_file,'w')
        fp.write(charmm_inp)
        fp.close()

        os.chdir(self.directory)
        cmd = '%s < %s > %s' % (charmming_config.charmm_exe,inp_file,out_file)
        os.system(cmd)

        try:
            fp = open(out_file, 'r')
        except:
            response['errcode'] = -999
            return response

        if not fp:
            response['errcode'] = -999
            return response

        ene = self.parse_energy(fp)
        fp.close()

        if ene < -9000000:
            response['errcode'] = -980
        else:
            response['errcode'] = 0  
            response['eneval']  = ene

        return response
