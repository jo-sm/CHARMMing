# -*- coding: utf-8 -*-
import sys
from train import load_model, run_prediction_regression

saved_model = sys.argv[1]
name = sys.argv[2]
activity_property = sys.argv[3]
out = sys.argv[4]
output_txt = sys.argv[5]

cmp = load_model(saved_model)
pred_r2 = run_prediction_regression(cmp,str(name),activity_property,str(out))

f = open(output_txt,"w")
f.write(str(pred_r2))
f.write("\n")
f.close()
