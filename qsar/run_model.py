# -*- coding: utf-8 -*-
import sys
from train import load_model, run_prediction

saved_model = sys.argv[1]
name = sys.argv[2]
threshold = sys.argv[3]
active = sys.argv[4]
inactive = sys.argv[5]
activity_property = sys.argv[6]
out = sys.argv[7]
output_txt = sys.argv[8]

cmp = load_model(saved_model)
recall,precision = run_prediction(cmp,str(name),float(threshold),active,inactive,activity_property,str(out))
recall = round(recall,3)
precision = round(precision,3)

f = open(output_txt,"w")
f.write(str(recall))
f.write("\n")
f.write(str(precision))
f.write("\n")
f.close()
