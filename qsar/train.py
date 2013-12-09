# -*- coding: utf-8 -*-
import random
import cPickle
from tempfile import mkstemp
import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.ML.DecTree.BuildSigTree import SigTreeBuilder
from rdkit.ML.Composite.Composite import Composite
from rdkit.ML.DecTree import CrossValidate
from rdkit.ML import ScreenComposite
from sklearn.ensemble import RandomForestRegressor
from sklearn import tree,metrics

def get_sd_properties(filename):
  path = filename
  props = []
  for x in Chem.SDMolSupplier(str(path)):
    if x is not None:
      props = list(x.GetPropNames())
      break
  if len(props) != 0:
      choices = [(p,p) for p in props]
  else:
      choices = []
  return choices

def is_valid_sd(filename):
    path = filename
    for x in Chem.SDMolSupplier(str(path)):
        if x is not None:
            return True
    return False
          
def active_inactive(ms,activity_property):
    types_count = dict()
    types = []
    for m in ms:
      value = m.GetProp(str(activity_property))
      if value in types_count:
        types_count[value] += 1
      else:
        types_count[value] = 1
    if len(types_count) != 2:
      return "","",True
    types = types_count.keys()
    if types_count[types[0]] > types_count[types[1]]:
      active = types[1]
      inactive = types[0]
    else:
      active = types[0]
      inactive = types[1]
    return active,inactive,False
    
def check_regression_property(ms,activity_property):
    types_count = dict()
    for m in ms:
        value = m.GetProp(str(activity_property))
        try:
            v = float(value)
        except ValueError:
            return False
        types_count[v] = 1
    if len(types_count) < 10:
        return False
    return True
                                                        
#YP
def check_activity_property(ms, activity_property): 
    types_count = dict()
    types = []
    for m in ms:
      value = m.GetProp(str(activity_property))
      if value in types_count:
        types_count[value] += 1
      else:
        types_count[value] = 1
    if len(types_count) != 2:
      return True
    
    return False

def train_model(pts):
    cmp = Composite()
    cmp.Grow(pts,attrs=[1],nPossibleVals=[2],nTries=500,randomDescriptors=3,
             buildDriver=CrossValidate.CrossValidationDriver,
                      treeBuilder=SigTreeBuilder,needsQuantization=False,
                      maxDepth=2,silent=True)
    return cmp

def train_model_regression(fps,acts):
    nclf = RandomForestRegressor(n_estimators=500, max_depth=5,random_state=0,n_jobs=4)
    nclf = nclf.fit(fps,acts)
    return nclf

def auc(pts,cmp):
    positive=0;
    decoy=0
    pa_act = []
    for p in pts:
      [pred,prob] = cmp.ClassifyExample(p)
      pa=prob
      if pred==0:
        pa=1-prob
      act = p[2]
      if act == 1:   
        positive+=1
      else:
        decoy+=1
      pa_act.append({"pa":pa,"act":act});
      
    sorted_pa_act = sorted(pa_act, key = lambda p: p['pa'])
    auc = 0.0
    left=decoy
    for p in sorted_pa_act:
      if p['act'] == 1:
        auc += (1-1.*left/decoy)
      else:
        left-=1
                              
    auc/=positive
    return auc

def r2(fps,acts,nclf):
    preds=nclf.predict(fps)
    r = metrics.r2_score(acts,preds)
    return r
    
def y_randomization(pts):
    acts = [x[2] for x in pts]
    random.shuffle(acts)
    rand_pts = []
    for i,p in enumerate(pts):
        rand_pts.append([p[0],p[1],acts[i]])
    rand_cmp = train_model(rand_pts)
    rand_auc = auc(rand_pts,rand_cmp)
    return rand_auc

def y_randomization_r2(fps,acts):
    rand_acts = acts[:]
    random.shuffle(rand_acts)
    rand_cmp = train_model_regression(fps,rand_acts)
    rand_r2 = r2(fps,rand_acts,rand_cmp)
    return rand_r2
                                    
                                    
def cross_validation(pts):
    random.shuffle(pts)
    avg_auc = 0
    for i in range(5):
        test_set = []
        train_set = []
        for j,p in enumerate(pts):
            if ( j % 5 == i):
                test_set.append(p)
            else:
                train_set.append(p)
        train_cmp = train_model(train_set)
        test_auc = auc(test_set,train_cmp)
        avg_auc += test_auc
    avg_auc /= 5
    return avg_auc

def cross_validation_r2(fps,acts):
    pts = zip(fps,acts)
    random.shuffle(pts)
    cross_fps,cross_acts = zip(*pts)   
    avg_r2 = 0
    for i in range(5):
        test_fp = [] 
        train_fp = []
        test_act = []
        train_act = []
        for j,p in enumerate(cross_fps):
            if ( j % 5 == i):
                test_fp.append(p)
            else:
                 train_fp.append(p)
        for j,a in enumerate(cross_acts):
            if ( j % 5 == i):
                test_act.append(a)
            else:
                train_act.append(a)
        train_cmp = train_model_regression(train_fp,train_act)
        test_r2 = r2(test_fp,test_act,train_cmp)
        avg_r2 += test_r2
    avg_r2 /= 5
    return avg_r2
                                                                                                                            

def calculate_threshold(pts,cmp):
    positive=0;
    decoy=0
    pa_act = []
    for p in pts:
      [pred,prob] = cmp.ClassifyExample(p)
      pa=prob
      if pred==0:
        pa=1-prob
      act = p[2]
      if act == 1:   
        positive+=1
      else:
        decoy+=1
      pa_act.append({"pa":pa,"act":act});

    recall = 1.
    precision = 0.
    count = 0
    recall,precision,threshold = recall_precision(0.,1.,pa_act,recall,precision,count)
    return recall,precision,threshold

def recall_precision(low,top,pa_act,old_recall,old_precision,count):
    count += 1
    threshold = 0.5*(low+top)
    actives=0
    predicted=0
    TP=0
    TN=0
    for p in pa_act:
        act=0
        if p['act']==1:
            actives+=1
            act=1
        pa=p['pa']
        if pa>=threshold:
            predicted+=1
            if act==1:
                TP+=1
        else:
            if act==0:
                TN+=1

    recall=1.*TP/actives
    if predicted != 0:
        precision=1.*TP/predicted
    else:
        precision=0
    if recall == old_recall and precision == old_precision and recall != 0 and precision != 0:
        return recall,precision,threshold
    if recall == precision and recall != 0:
        return recall,precision,threshold
    if count > sys.getrecursionlimit() / 2:
            return recall,precision,threshold
    if recall > precision:
        return recall_precision(threshold,top,pa_act,recall,precision,count)
    else:
        return recall_precision(low,threshold,pa_act,recall,precision,count)
        

def load_model(name):
    cmp = cPickle.load(file(name,'rb'))
    return cmp

def run_prediction(cmp,name,threshold,active,inactive,activity_property,out):
    actives=0
    predicted=0
    TP=0
    i=0
    w = Chem.SDWriter(out)
    for m in  Chem.SDMolSupplier(name) :
        if m is not None:
            act=0
            fp=AllChem.GetMorganFingerprintAsBitVect(m,2,2048)
            [pred,prob] = cmp.ClassifyExample([i,fp,act])
            pa=prob
            if pred==0:
                pa=1-prob
            m.SetProp('PREDICTED_SCORE',str(pa))
            if pa > threshold:
                m.SetProp('PREDICTED_'+str(activity_property),str(active))
            else:
                m.SetProp('PREDICTED_'+str(activity_property),str(inactive))
            w.write(m)
            i+=1
           
            if m.HasProp(str(activity_property)):
                real_act = m.GetProp(str(activity_property))
                if real_act == str(active):
                    actives += 1
                    if pa > threshold:
                        TP += 1
                if pa > threshold:
                    predicted += 1
    if actives != 0:
        recall=1.*TP/actives
    else:
        recall=0.
    if predicted != 0:
        precision=1.*TP/predicted
    else:
        precision=0.
    return recall, precision

def run_prediction_regression(nclf,name,activity_property,out):
    w = Chem.SDWriter(out)
    preds = []
    real_acts = []
    count = 0
    for m in  Chem.SDMolSupplier(name) :
        if m is not None:
            fp=AllChem.GetMorganFingerprintAsBitVect(m,2,2048)                        
            fps = [fp]
            pred=nclf.predict(fps)
            m.SetProp('PREDICTED_'+str(activity_property),str(pred[0]))
            w.write(m)
            if m.HasProp(str(activity_property)):
                real_act = m.GetProp(str(activity_property))
                preds.append(pred[0])
                real_acts.append(float(real_act))
                count += 1
    r = 0.
    if count > 10:
        r = metrics.r2_score(real_acts,preds)
    return r
    
