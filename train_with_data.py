from laldetchar.idq import idq
import numpy as np
import os
import ConfigParser
config = ConfigParser.SafeConfigParser()
config.read("/home/jessica.hyde/lalsuites/src/laldetchar/python/laldetchar/idq/idq.ini")
pat = "/home/jessica.hyde/test/dev-branch/now//train/1117329856_1117331656/L1_mlaA-1117329856-207080.pat" 
traindir = config.get('general', 'traindir')
gpsstart = 1117329856 
gpsstop = 1117331656
cwd = os.getcwd()
classifiersD, mla, ovl = idq.config_to_classifiersD( config )

classifiers = sorted(classifiersD.keys())
classD = classifiersD['mvsc']
miniconfig = classD['config']
usertag = "A-both"
output_dir = "/home/jessica.hyde/test/dev-branch/now//train/1117329856_1117331656/"
dag_classifiers = ["mvsc"]
for classifier in dag_classifiers: ### these are trained with a dag
    classD = classifiersD[classifier]
    train_cache = dict( (classifier, idq.Cachefile(idq.cache(traindir, classifier, tag='_train%s'%usertag))))
    train_dir = "%s/%s/"%(output_dir, classifier)
    flavor = classD['flavor']
    if not os.path.exists(train_dir):
        os.makedirs(train_dir)
    else:
        os.system("rm %s/*.dag*"%train_dir) ### remove existing condor files!
    (submit_dag_exit_status, dag_file) = idq.dag_train(flavor, pat,  train_cache[classifier], miniconfig, train_dir, cwd)

