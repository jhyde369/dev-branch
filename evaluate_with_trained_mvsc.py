from laldetchar.idq import idq
import numpy as np
import ConfigParser
from laldetchar.idq import time_lock as tl
import plotROC as pR


def evaluate_all_combinations(patfile, trainedforestbase, trainedforesttend, config):
    original_usertag = config.get("general", "usertag") ### remember what the tag started as
    original_variables = config.get("mvsc", "z") ### remember the variables
    print (original_variables)
    for tstart in ["no_tstart", "yes_tstart"]: ### iterate over all possible combinations
        for tend in ["no_tend", "yes_tend"]:
            config.set("general", "usertag", "%s-%s-%s"%(original_usertag, tstart, tend) ) ### change the usertag as needed
            if(tstart == "no_tstart" and tend == "no_tend"):
                config.set("mvsc", "z", "%s,%s,%s"%(original_variables, "time_locked", "time_end"))
            elif(tstart == "no_tstart"):
                config.set("mvsc", "z", "%s,%s"%(original_variables, "time_locked"))                
            elif(tend == "no_tend"):
                config.set("mvsc", "z", "%s,%s"%(original_variables, "time_end"))
            else:
                config.set("mvsc","z", original_variables)
            trainedforest = trainedforestbase + "%s-%s/"%(tstart,tend) + trainedforestend
            ranked_file = trainedforest[:-3] + "dat"
            gps_start_time = -np.infty
            gps_end_time = np.infty
            dir = "."
            classifiersD, mla, ovl = idq.config_to_classifiersD( config )
            classD = classifiersD['mvsc']
            miniconfig = classD['config']
            print(ranked_file)
            print(trainedforest)
            print(config.get("mvsc", "z"))
            idq.forest_evaluate(patfile, trainedforest, ranked_file, miniconfig, gps_start_time, gps_end_time, dir)
    config.set("general", "usertag", original_usertag)
    config.set("mvsc", "z", original_variables)


config = ConfigParser.SafeConfigParser()
patfile = "/home/jessica.hyde/test/dev-branch/now//train/1117329857_1118336935//L1_mla-1117329857-1007078-B.pat"
config.read("/home/jessica.hyde/lalsuites/src/laldetchar/python/laldetchar/idq/idq.ini")
trainedforestbase = "/home/jessica.hyde/test/dev-branch/now//train/1117329857_1118336935//mvsc/"
trainedforestend = "L1_mla-1117329857-1007078-A.spr"

#evaluate_all_combinations(patfile, trainedforestbase, trainedforestend, config)
#pR.plot_ROC_combinations("ER7_ROC_Curves-A.png", 1117329857, 1118336935)
gpsstart = 1117329857 
gpsstop = 1118336935 
### get start and end seg times formatted as list of lists
tl.get_segments(gpsstart, gpsstop)
scisegs = tl.parse_segments("segment_times.txt")

newsegs = []
for seg in scisegs:
    if (seg[1]-seg[0] >= 3600):
        newsegs.append(seg)
scisegs=newsegs

for seg in scisegs:
    patf = idq.pat("/home/jessica.hyde/test/dev-branch/now/train/%s_%s/" % (gpsstart, gpsstop), "L1", "LowThresh", seg[0], seg[1]-seg[0])
    patfile = patf
    trainedforestbase = "/home/jessica.hyde/test/dev-branch/now/train/%s_%s/mvsc/" % (gpsstart, gpsstop)
    trainedforestend = "L1_mla-%s-%s-A.spr" % (int(seg[0]), int(seg[1]-seg[0]))
    evaluate_all_combinations(patfile, trainedforestbase, trainedforestend, config)
    pR.plot_ROC_combinations("LowThreshROCSegment%s-%s-A.png" % (seg[0], seg[1]), gpsstart, gpsstop, seg[0], seg[1])
#    trainedforestend = "L1_mla-%s-%s-B.spr" % (int(seg[0]), int(seg[1]-seg[0]))
#    evaluate_all_combinations(patfile, trainedforestbase, trainedforestend, config)
#    pR.plot_ROC_combinations("LowThreshROCSegment%s-%s-B.png" % (seg[0], seg[1]), gpsstart, gpsstop, seg[0], seg[1])
    

"""for seg in scisegs:
    patf = idq.pat("/home/jessica.hyde/test/dev-branch/now/train/%s_%s/" % (gpsstart, gpsstop), "L1", "", seg[0], seg[1]-seg[0])
    patfile = patf[:-4] + "-B.pat"
    trainedforestbase = "/home/jessica.hyde/test/dev-branch/now/train/%s_%s/mvsc/" % (gpsstart, gpsstop)
    trainedforestend = "L1_mla-%s-%s-A.spr" % (int(seg[0]), int(seg[1]-seg[0]))
    evaluate_all_combinations(patfile, trainedforestbase, trainedforestend, config)
    pR.plot_ROC_combinations("ROCSegment%s-%s.png" % (seg[0], seg[1]), gpsstart, gpsstop, seg[0], seg[1])"""
