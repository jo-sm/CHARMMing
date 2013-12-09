# -*- coding: utf-8 -*-
import cPickle
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
from train import train_model_regression,r2,y_randomization_r2,cross_validation_r2
import common

input_sdf = sys.argv[1]
output_model = sys.argv[2]
activity_property = sys.argv[3]
output_txt = sys.argv[4]
model_id = sys.argv[5]
user_id = sys.argv[6]

fps = []
acts = []
for m in Chem.SDMolSupplier(str(input_sdf)):
	if m is not None:
		fp = AllChem.GetMorganFingerprintAsBitVect(m,2,2048)
        	act = float(m.GetProp(str(activity_property)))
        	fps.append(fp)
        	acts.append(act)
        	
cmp = train_model_regression(fps,acts)
cPickle.dump(cmp,file(output_model,'wb+'));
self_r2 = r2(fps,acts,cmp)
self_r2 = round(self_r2,3)
r2_rand = y_randomization_r2(fps,acts)
r2_rand = round(r2_rand,3)
cross_r2 = cross_validation_r2(fps,acts)
cross_r2 = round(cross_r2,3)
f = open(output_txt,"w")
f.write(str(self_r2))
f.write("\n")
common.AssignObjectAttribute(user_id, model_id, "qsar_qsar_models", "Self R2", str(self_r2))
common.AssignObjectAttribute(user_id, model_id, "qsar_qsar_models", "Y-randomization", str(r2_rand))
common.AssignObjectAttribute(user_id, model_id, "qsar_qsar_models", "5-fold cross-validation", str(cross_r2))

f.write(str(r2_rand))
f.write("\n")
f.write(str(cross_r2))
f.write("\n")
f.close()
