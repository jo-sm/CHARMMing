# -*- coding: utf-8 -*-
import cPickle
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
from train import train_model,auc,y_randomization,cross_validation,calculate_threshold
#from qsar.models import qsar_models
import common

input_sdf = sys.argv[1]
output_model = sys.argv[2]
activity_property = sys.argv[3]
active = sys.argv[4]
output_txt = sys.argv[5]
model_id = sys.argv[6]
user_id = sys.argv[7]
type = sys.argv[8]


pts = []
i = 0
for m in Chem.SDMolSupplier(str(input_sdf)):
	if m is not None:
		fp = AllChem.GetMorganFingerprintAsBitVect(m,2,2048)
        	if m.GetProp(str(activity_property)) == active:
        		act=1
        	else:
        		act=0
                i += 1
        	pts.append([i,fp,act])

cmp = train_model(pts,type)
cPickle.dump(cmp,file(output_model,'wb+'));
self_auc = auc(pts,cmp,type)
self_auc = round(self_auc,3)
recall,precision,threshold = calculate_threshold(pts,cmp,type)
recall = round(recall,3)
precision = round(precision,3)
auc_rand = y_randomization(pts,type)
auc_rand = round(auc_rand,3)
cross_auc = cross_validation(pts,type)
cross_auc = round(cross_auc,3)
f = open(output_txt,"w")
f.write(str(self_auc))
f.write("\n")
f.write(str(recall))
f.write("\n")
f.write(str(precision))
f.write("\n")
f.write(str(threshold))
f.write("\n")
#qsar_model=qsar_models.objects.get(id=model_id)
common.AssignObjectAttribute(user_id, model_id, "qsar_qsar_models", "Recommended threshold", str(threshold))
common.AssignObjectAttribute(user_id, model_id, "qsar_qsar_models", "Self AUC", str(self_auc))
common.AssignObjectAttribute(user_id, model_id, "qsar_qsar_models", "Y-randomization", str(auc_rand))
common.AssignObjectAttribute(user_id, model_id, "qsar_qsar_models", "5-fold cross-validation", str(cross_auc))
common.AssignObjectAttribute(user_id, model_id, "qsar_qsar_models", "Recall", str(recall))
common.AssignObjectAttribute(user_id, model_id, "qsar_qsar_models", "Precision", str(precision))

f.write(str(auc_rand))
f.write("\n")
f.write(str(cross_auc))
f.write("\n")
f.close()
