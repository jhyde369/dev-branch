import time
import sys
import os
import traceback
import logging
from random import shuffle
import numpy
import copy

import ConfigParser
from optparse import OptionParser

import subprocess as sp
import multiprocessing as mp

from laldetchar.idq import idq
#from laldetchar.idq import reed as idq
from laldetchar.idq import event
from laldetchar.idq import auxmvc_utils
from laldetchar.idq import time_lock as tl
from laldetchar.idq import pat_shuffle

from pylal import frutils

from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables
from glue.ligolw import table

from laldetchar import git_version

########################################

__author__ = \
    'Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>), Jessica Hyde (<jessica.hyde@ligo.org>)'
__version__ = git_version.id
__date__ = git_version.date

description = \
    """This program generates .pat files with enriched features such as the time since the beginning of the segment. """

#===================================================================================================

parser = OptionParser(version='Name: %%prog\n%s'
                      % git_version.verbose_msg, usage='%prog [options]'
                      , description=description)
parser.add_option('-c', '--config', default='idq.ini', type='string', help='configuration file')

parser.add_option('-s', '--gpsstart', dest="gpsstart", default=False, type='int', help='a GPS start time for the analysis. If default, gpsstart is calculated from the current time.')
parser.add_option('-e', '--gpsstop', dest="gpsstop", default=False, type='int', help='a GPS stop time for the analysis. If default, gpsstop is calculated from the current time.')
parser.add_option('-b', '--lookback', default='0', type='string', help="Number of seconds to look back and get data for training. Default is zero.\
    Can be either positive integer or 'infinity'. In the latter case, the lookback will be incremented at every stride and all data after --gps-start will be used in every training.")

parser.add_option('-l', '--log-file', default='idq_train.log', type='string', help='log file')

parser.add_option("", "--ignore-science-segments", default=False, action="store_true", help="ignores science segments and trains over all available data within range")

parser.add_option("", "--include-time-locked", default=False, action="store_true", help="includes the time since the beginning of the segment for each trigger")

parser.add_option("", "--include-time-until-end", default=False, action="store_true", help="includes the time until the end of the segment for each trigger")

parser.add_option('-f','--force',default=False, action='store_true', help="forces classifiers to be trained, \
    and if errors are recovered it passes them along. Intended for use when we must have a trained classifier or else we should raise an error. \
    This may produce an infinite loop if the condor dags get lost (for mla classifiers). Use with caution.")
parser.add_option('', '--min-seg-length', dest = "min_seg_length", default = None, type='int', help="the minimum segment length to be considered")
(opts, args) = parser.parse_args()

if opts.lookback != "infinity":
    opts.lookback = int(opts.lookback)
print(opts.include_time_locked)
print(opts.include_time_until_end)
##################################################

#===================================================================================================
### setup logger to record processes
logger = idq.setup_logger('idq_logger', opts.log_file, sys.stdout, format='%(asctime)s %(message)s')

sys.stdout = idq.LogFile(logger)
sys.stderr = idq.LogFile(logger)

#================================================================================================
### read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(opts.config)

ifo = config.get('general', 'ifo')

usertag = config.get('general', 'usertag')
if usertag:
    usertag = "_%s"%usertag

#================================================================================================

auxmvc_coinc_window = config.getfloat('realtime', 'padding')
auxmc_gw_signif_thr = config.getfloat('general', 'gw_kwsignif_thr')

auxmvc_selected_channels = config.get('general','selected-channels')
auxmvc_unsafe_channels = config.get('general','unsafe-channels')

#========================
# realtime
#========================
realtimedir = config.get('general', 'realtimedir')

padding = config.getfloat('realtime', 'padding')

clean_rate = config.getfloat('realtime', 'clean_rate')
clean_window = config.getfloat('realtime', 'clean_window')
clean_threshold = config.getfloat('realtime', 'clean_threshold')

#========================
# train
#========================
traindir = config.get('general', 'traindir')
"""if ovl: ### need snglchandir
   snglchndir = config.get('general', 'snglchndir') """

stride = config.getint('train', 'stride')
delay = config.getint('train', 'delay')

#train_script = config.get('condor', 'train')


#build_auxmvc_vectors = mla and (not os.path.exists(realtimedir)) ### if realtimedir does not exist, we cannot rely on patfiles from the realtime job
                                                                 ### we need to build our own auxmvc_vectors

max_gch_samples = config.getint("train", "max-glitch-samples")
max_cln_samples = config.getint("train", "max-clean-samples")

#========================
# data discovery
#========================
if not opts.ignore_science_segments:
    ### load settings for accessing dmt segment files
#    dmt_segments_location = config.get('get_science_segments', 'xmlurl')
    dq_name = config.get('get_science_segments', 'include')
#    dq_name = config.get('get_science_segments', 'include').split(':')[1]
    segdb_url = config.get('get_science_segments', 'segdb')
else:
    dq_name = None ### needs to be defined

### define the gravitational wave channel/what is a glitch
gwchannel = config.get('general', 'gwchannel')
gwthreshold = config.getfloat('general', 'gw_kwsignif_thr')

### kleineWelle config
GWkwconfigpath = config.get('data_discovery', 'GWkwconfig')
GWkwconfig = idq.loadkwconfig(GWkwconfigpath)
GWkwbasename = GWkwconfig['basename']
GWgdsdir = config.get('data_discovery', 'GWgdsdir')
GWkwstride = int(float(GWkwconfig['stride']))

GWkwtrgdir = "%s/%s"%(GWgdsdir, GWkwbasename)

AUXkwconfigpath = config.get('data_discovery', 'AUXkwconfig')
AUXkwconfig = idq.loadkwconfig(AUXkwconfigpath)
AUXkwbasename = AUXkwconfig['basename']
AUXgdsdir = config.get('data_discovery', 'AUXgdsdir')
AUXkwstride = int(float(AUXkwconfig['stride']))

AUXkwtrgdir = "%s/%s"%(AUXgdsdir, AUXkwbasename)

identical_trgfile = (GWgdsdir == AUXgdsdir) and (GWkwbasename == AUXkwbasename) ### used to avoid re-loading the same file for both GW and AUX triggers

#==================================================
### current time and boundaries

t = int(idq.nowgps())

cwd = os.getcwd()
gpsstop = opts.gpsstop
if not gpsstop: ### stop time of this analysis
    logger.info('computing gpsstop from current time')
    gpsstop = t ### We do not require boundaries to be integer multiples of stride

gpsstart = opts.gpsstart
if not gpsstart:
    logger.info('computing gpsstart from gpsstop')
    gpsstart = gpsstop - stride

lookback = 0

### directory into which we write data
output_dir = "%s/%d_%d/"%(traindir, gpsstart, gpsstop)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)



##### Get Science Segments

logger.info('test - Begin: querying science segments')
"""if not opts.ignore_science_segments:
    try:
        ### this returns a string
        seg_xml_file = idq.segment_query(config, gpsstart - lookback , gpsstop, url=segdb_url)

        ### load xml document
        xmldoc = ligolw_utils.load_fileobj(seg_xml_file, contenthandler=ligolw.LIGOLWContentHandler)[0]

        ### science segments xml filename
        seg_file = idq.segxml(output_dir, "_%s"%dq_name, gpsstart - lookback , lookback+stride)

        logger.info('writing science segments to file : '+seg_file)
        ligolw_utils.write_filename(xmldoc, seg_file, gz=seg_file.endswith(".gz"))

    except Exception as e:
        traceback.print_exc()
        logger.info('ERROR: segment generation failed. Skipping this training period.')

    logger.info('finished segment generation')

    (scisegs, coveredseg) = idq.extract_dq_segments(seg_file, dq_name)
    print ("Scisegs: " + str(scisegs))
else:
    scisegs=None"""

### get segments and write to xml doc
gps_start_time = gpsstart
gps_end_time = gpsstop
print("start: " +str( gps_start_time) + ", end: " +str( gps_end_time))
proc = sp.Popen(["ligolw_segment_query_dqsegdb", "-t", "https://dqsegdb5.phy.syr.edu", "-q", "-a", "L1:DMT-ANALYSIS_READY:1", "-s", str(gps_start_time), "-e", str(gps_end_time), "-o", "sample_segments.xml", "-S"], stdout = sp.PIPE)
ans = proc.communicate()[0]
### get start and end times from xml doc, write to output doc
output_file = open("segment_times.txt", "w")
proc = sp.Popen(["ligolw_print", "-c" "start_time", "-c" "end_time", "-t" "segment", "sample_segments.xml"], stdout=output_file)
ans = proc.communicate()[0]
output_file.close()

### get start and end seg times formatted as list of lists
scisegs = tl.parse_segments("segment_times.txt")
print("Segs: " + str(scisegs))
if opts.ignore_science_segments:
    scisegs=None

### remove segments that are too short
if opts.min_seg_length:
    newsegs=[]
    for i in range(len(scisegs)):
        if (scisegs[i][1]-scisegs[i][0] >= opts.min_seg_length):
            newsegs.append(scisegs[i])
    scisegs = newsegs
    print("New Segs: " + str(scisegs))
#==================================================================

trigger_dict = idq.retrieve_kwtrigs(GWgdsdir, GWkwbasename, gpsstart, gpsstop-gpsstart, GWkwstride, sleep=0, ntrials=1, logger=logger, segments=scisegs) ### go find GW kwtrgs
if gwchannel not in trigger_dict:
    trigger_dict[gwchannel] = []
print("Number trigs before theshholding: " +  str(len(trigger_dict[gwchannel])))
### keep only relevant gchs
if len(trigger_dict[gwchannel]) > max_gch_samples:
    trigger_dict.resort() ### make sure they're in the correct order
    trigger_dict[gwchannel] = trigger_dict[gwchannel][-max_gch_samples:]
print("Number of trigs after max cap: " + str(len(trigger_dict[gwchannel])))


if not identical_trgfile: ### add AUX triggers
    logger.info('looking for additional AUX triggers')
    aux_trgdict = idq.retrieve_kwtrigs(AUXgdsdir, AUXkwbasename, gpsstart-lookback-padding, lookback+stride+padding, AUXkwstride, sleep=0, ntrials=1, logger=logger, segments=scisegs) ### find AUX kwtrgs
    if aux_trgdict == None:
        logger.warning('  no auxiliary triggers were found')
        ### we do not skip, although we might want to?
    else:
        triggger_dict.add( aux_trgdict )
        trigger_dict.resort()

trigger_dict.include([[gpsstart-lookback-padding, gpsstop+padding]])
trigger_dict.include([[gpsstart-lookback, gpsstop]], channels=[gwchannel])
print("Number of trigs after windowing: " + str(len(trigger_dict[gwchannel])))

### define cleans
logger.info('  constructing cleans')
dirtyseg = event.vetosegs(trigger_dict[gwchannel], clean_window, clean_threshold)

if not opts.ignore_science_segments:
    clean_gps = sorted(event.randomrate(clean_rate, scisegs)) ### generate random clean times as a poisson time series within scisegs
else:
    clean_gps = sorted(event.randomrate(clean_rate, [[gpsstart-lookback, gpsstart + stride]])) ### generate random clean times as a poisson time series within analysis range
clean_gps = [ l[0] for l in event.exclude( [[gps] for gps in clean_gps], dirtyseg, tcent=0)] ### keep only those gps times that are outside of dirtyseg
print("Number of cleans: " + str(len(clean_gps)))

## keep only the most relevant cleans
if len(clean_gps) > max_cln_samples:
    clean_gps = clean_gps[-max_cln_samples:]
print("Number of cleans after max cap: " + str(len(clean_gps)))

### keep only times that are within science time
if not opts.ignore_science_segments:
    logger.info('  filtering trigger_dict through scisegs')
    trigger_dict.include(scisegs) ### already loaded into memory above here
print("Number of trigs in scisegs: " + str(len(trigger_dict[gwchannel])))

"""gwtrigs = trigger_dict[gwchannel][:]
shuffle(gwtrigs)
print("Shuffled")
firsthalf = gwtrigs[:len(gwtrigs)/2]
secondhalf = gwtrigs[len(gwtrigs)/2:]
print("Starting copies")
first_trigger_dict = copy.deepcopy(trigger_dict)
print("First Copy")
second_trigger_dict = copy.deepcopy(trigger_dict)
print("Second Copy")
first_trigger_dict[gwchannel] = firsthalf
second_trigger_dict[gwchannel] = secondhalf"""

print("Number of triggers def into build_auxmvc_vectors(): " + str(len(trigger_dict[gwchannel])))
### build vectors, also writes them into pat
pat = idq.pat(output_dir, ifo, usertag, gpsstart-lookback, lookback+gpsstop-gpsstart)
logger.info('  writing %s'%pat)
idq.build_auxmvc_vectors(trigger_dict, gwchannel, auxmvc_coinc_window, auxmc_gw_signif_thr, pat, gps_start_time=gpsstart,
    gps_end_time=gpsstop,  channels=auxmvc_selected_channels, unsafe_channels=auxmvc_unsafe_channels, clean_times=clean_gps,
    clean_window=clean_window, filter_out_unclean=False, locked_segments = "segment_times.txt", include_time_locked = opts.include_time_locked, 
    include_time_until_end = opts.include_time_until_end)
print("Done with pat file building")
ptas_exit_status = 0 ### used to check for success

if ptas_exit_status != 0: ### check that process executed correctly
    logger.warning('WARNING: Preparing training auxmvc samples failed')
    if opts.force:
        raise StandardError, "auxmvc samples required for successful training"
    else:
        logger.warning('WARNING: skipping re-training the MLA classifiers')
else:
    ### figure out training set size
    ### load auxmvc vector samples
    
    logger.info('Done: preparing training auxmvc samples')

    (firstpat, secondpat) = pat_shuffle.split_pat_file(pat)

    for pat in [firstpat, secondpat]:
        auxmvc_samples = auxmvc_utils.ReadMVSCTriggers([pat], Classified=False)

        ### get numbers of samples
        N_clean = len(auxmvc_samples[numpy.nonzero(auxmvc_samples['i']==0)[0], :])
        N_glitch = len(auxmvc_samples[numpy.nonzero(auxmvc_samples['i']==1)[0], :])

        del auxmvc_samples

        #================================================================================
        dags = {}
        min_num_cln = int(config.get("mvsc", "min_num_cln"))
        min_num_gch = int(config.get("mvsc", "min_num_gch"))
        print(min_num_cln, min_num_gch, N_clean, N_glitch)
        classifier = "mvsc"
        flavor = "forest"
        #if opts.force or ((N_clean >= min_num_cln) and (N_glitch >= min_num_gch)):
        if opts.force or (N_clean >= min_num_cln and N_glitch >= min_num_gch):
            train_dir = "%s/%s/"%(output_dir, classifier)
            if not os.path.exists(train_dir):
                os.makedirs(train_dir)
            else:
                os.system("rm %s/*.dag*"%train_dir) ### remove existing condor files!

            original_usertag = config.get("general", "usertag") ### remember what the tag started as
            original_variables = config.get("mvsc", "z") ### remember the variables
            print (original_variables)
            for tstart in ["no_tstart", "yes_tstart"]: ### iterate over all possible combinations
                for tend in ["no_tend", "yes_tend"]:
                    config.set("general", "usertag", "%s-%s-%s"%(usertag+original_usertag, tstart, tend) ) ### change the usertag as needed
                    if(tstart == "no_tstart" and tend == "no_tend"):
                        config.set("mvsc", "z", "%s,%s,%s"%(original_variables, "time_locked", "time_end"))
                    elif(tstart == "no_tstart"):
                        config.set("mvsc", "z", "%s,%s"%(original_variables, "time_locked"))                
                    elif(tend == "no_tend"):
                        config.set("mvsc", "z", "%s,%s"%(original_variables, "time_end"))
                    else:
                        config.set("mvsc","z", original_variables)
                    logger.info('submitting %s training dag'%classifier)

                    classifiersD, mla, ovl = idq.config_to_classifiersD( config )

                    classD = classifiersD["mvsc"]

                    train_cache = dict( (classifier, idq.Cachefile(idq.cache(traindir, classifier, tag='_train%s'%usertag))) for classifier in ["mvsc"])
                    ### submit training job
                    miniconfig = classD['config']
                    newtrain_dir = train_dir+"/"+"%s-%s/"%(tstart, tend)
                    if not os.path.exists(newtrain_dir):
                        os.makedirs(newtrain_dir)
                    else:
                        os.system("rm %s/*.dag*"%newtrain_dir) 
                    (submit_dag_exit_status, dag_file) = idq.dag_train(flavor, pat,  train_cache[classifier], miniconfig, newtrain_dir, cwd)

                    if submit_dag_exit_status:
                        logger.warning("WARNING: was not able to submit %s training dag"%classifier)
                        dags[dag_file] = "submit failed"
                        if opts.force:
                            raise StandardError, "submission of %s training dag failed"%classifier
                        else:
                            dags[dag_file] = "incomplete"
                    print("Train dir: " + train_dir)
                    print("Dag File: " + dag_file)
                    print(config.get("mvsc","z"))
                    logger.info('dag file %s/%s'%(train_dir, dag_file))
            config.set("general", "usertag", usertag)
            config.set("mvsc", "z", original_variables)
        elif (N_clean >= min_num_cln):
            logger.warning("WARNING: not enough glitches in training set. skipping %s training"%"mvsc")
        elif (N_glitch >= min_num_gch):
            logger.warning("WARNING: not enough cleans in training set. skipping %s training"%"mvsc")
        else:
            logger.warning("WARNING: neither enough cleans nor enough glitches in training set. skipping %s training"%"mvsc")




"""status = (0,0)
if opts.include_time_locked:
    status[0] = 1
if opts.include_time_until_end:
    status[1] = 1
ranked_file = pat[:-3]+str(status)+"dat"
trainedforest = "/home/jessica.hyde/test/dev-branch/now/train/1117258800_1117260600/mvsc/L1_mla-1117258800-1800.spr"
gps_start_time = -np.infty
gps_end_time = np.infty
dir = "."
idq.forest_evaluate(pat, dag_file, ranked_file, miniconfig, gps_start_time, gps_end_time, dir)"""
