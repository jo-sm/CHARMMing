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
from sklearn import tree,metrics
from sklearn import svm,linear_model,ensemble,tree
#from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB

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

def train_model(pts,type):
    if type != "rf":
        return train_model_other(pts,type)
    return train_model_rf(pts)
        
def train_model_rf(pts):
    cmp = Composite()
    cmp.Grow(pts,attrs=[1],nPossibleVals=[2],nTries=500,randomDescriptors=3,
             buildDriver=CrossValidate.CrossValidationDriver,
                      treeBuilder=SigTreeBuilder,needsQuantization=False,
                      maxDepth=2,silent=True)
    return cmp

def train_model_other(pts,type):
    fps = [x[1] for x in pts]
    act = [x[2] for x in pts]
    clf = eval(type+"()")
    clf.fit(fps, act)  
    return clf
    
def classify(p,cmp,type):
    if type != "rf":
        fp = p[1]
        a = cmp.predict(fp);
        pred = a[0]
        if "decision_function" in dir(cmp):
            b = cmp.decision_function(fp);
        else:
            b = cmp.predict_proba(fp);
        try:
            prob = b[0][0]
        except IndexError:
            try:
                prob = b[0]
            except IndexError:
                prob = b
        if "decision_function" in dir(cmp):
            if pred == 0:
                prob = 1-prob
        else:
            if pred == 1:
                prob = 1-prob
        return pred,prob
    return cmp.ClassifyExample(p)
    
def train_model_regression(fps,acts,type):
    nclf = eval(type+"()")
#    RandomForestRegressor(n_estimators=500, max_depth=5,random_state=0,n_jobs=4)
    nclf = nclf.fit(fps,acts)
    return nclf
    
def auc(pts,cmp,type):
    positive=0;
    decoy=0
    pa_act = []
    for p in pts:
      [pred,prob] = classify(p,cmp,type) #cmp.ClassifyExample(p)
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
    
def y_randomization(pts,type):
    acts = [x[2] for x in pts]
    random.shuffle(acts)
    rand_pts = []
    for i,p in enumerate(pts):
        rand_pts.append([p[0],p[1],acts[i]])
    rand_cmp = train_model(rand_pts,type)
    rand_auc = auc(rand_pts,rand_cmp,type)
    return rand_auc

def y_randomization_r2(fps,acts,type):
    rand_acts = acts[:]
    random.shuffle(rand_acts)
    rand_cmp = train_model_regression(fps,rand_acts,type)
    rand_r2 = r2(fps,rand_acts,rand_cmp)
    return rand_r2
                                    
                                    
def cross_validation(pts,type):
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
        train_cmp = train_model(train_set,type)
        test_auc = auc(test_set,train_cmp,type)
        avg_auc += test_auc
    avg_auc /= 5
    return avg_auc

def cross_validation_r2(fps,acts,type):
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
        train_cmp = train_model_regression(train_fp,train_act,type)
        test_r2 = r2(test_fp,test_act,train_cmp)
        avg_r2 += test_r2
    avg_r2 /= 5
    return avg_r2
                                                                                                                            

def calculate_threshold(pts,cmp,type):
    positive=0;
    decoy=0
    pa_act = []
    low = sys.float_info.max
    top = -(sys.float_info.max/2)
    for p in pts:
      [pred,prob] = classify(p,cmp,type) #cmp.ClassifyExample(p)
      pa=prob
      if prob > top:
          top = prob
      if prob < low:
          low = prob
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
    recall,precision,threshold = recall_precision(low,top,pa_act,recall,precision,count)
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

def run_prediction(cmp,name,threshold,active,inactive,activity_property,out,type):
    actives=0
    predicted=0
    TP=0
    i=0
    w = Chem.SDWriter(out)
    for m in  Chem.SDMolSupplier(name) :
        if m is not None:
            act=0
            fp=AllChem.GetMorganFingerprintAsBitVect(m,2,2048)
            [pred,prob] = classify([i,fp,act],cmp,type) #cmp.ClassifyExample([i,fp,act])
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

name_to_type_cat = dict()
name_to_type_cat["Random Forest (SAR) Categorization"] = ["rf",True]
name_to_type_cat["SVM (SAR) Categorization"] = ["svm.SVC",True]
name_to_type_cat["Random Forest (QSAR) Regression"] = ["ensemble.RandomForestRegressor",False]
name_to_type_cat["SVM (QSAR) Regression"] = ["svm.SVR",False]
name_to_type_cat["Logit (SAR) Categorization"] = ["linear_model.LogisticRegression",True]
name_to_type_cat["SGD (SAR) Categorization"] = ["linear_model.SGDClassifier",True]
#name_to_type_cat["Nearest Neighbors (SAR) Categorization"] = ["KNeighborsClassifier",True]
name_to_type_cat["Naive Bayes (SAR) Categorization"] = ["GaussianNB", True]
#name_to_type_cat["AdaBoost (SAR) Categorization"] = ["ensemble.AdaBoostClassifier", True]
name_to_type_cat["Gradient Boosting (SAR) Categorization"] = ["ensemble.GradientBoostingClassifier", True]
name_to_type_cat["Decision Tree (SAR) Categorization"] = ["tree.DecisionTreeClassifier",True]
name_to_type_cat["Gradient Boosting (QSAR) Regression"] = ["ensemble.GradientBoostingRegressor",False]
name_to_type_cat["Decision Tree (QSAR) Regression"] = ["tree.DecisionTreeRegressor",False]
#name_to_type_cat["Least Squares Linear (QSAR) Regression"] = ["linear_model.LinearRegression",False]
name_to_type_cat["Ridge (QSAR) Regression"] = ["linear_model.Ridge",False]
name_to_type_cat["Lasso (QSAR) Regression"] = ["linear_model.Lasso",False]
name_to_type_cat["Elastic Net (QSAR) Regression"] = ["linear_model.ElasticNet",False]
name_to_type_cat["SGD (QSAR) Regression"] = ["linear_model.SGDRegressor",False]

def categorization_regression(name):
    type,category = name_to_type_cat[name]
    return type,category
