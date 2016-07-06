from laldetchar.idq import idq
from laldetchar.idq import time_lock as tl
import ConfigParser
from optparse import OptionParser
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import subprocess as sp
from laldetchar import git_version
import sys
import logging
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
parser.add_option('-c', '--config', default='laldetchar/idq/idq.ini', type='string', help='configuration file')

parser.add_option('-s', '--gpsstart', dest="gpsstart", default=False, type='int', help='a GPS start time for the analysis. If default, gpsstart is calculated from the current time.')
parser.add_option('-e', '--gpsstop', dest="gpsstop", default=False, type='int', help='a GPS stop time for the analysis. If default, gpsstop is calculated from the current time.')
parser.add_option('', '--min-seg-length', dest="min_seg_length", default = None, type='int', help="minimum length for included segments")
parser.add_option('-l', '--log-file', default='idq_train.log', type='string', help='log file')
parser.add_option('-t', '--threshold', default=None, dest="threshold", type='int', help="significance threshold for glitches")
parser.add_option('', '--tag', default="", dest="tag", type='string', help="a tag to label filenames with")
(opts, args) = parser.parse_args()

#===================================================================================================
### setup logger to record processes
logger = idq.setup_logger('idq_logger', opts.log_file, sys.stdout, format='%(asctime)s %(message)s')

sys.stdout = idq.LogFile(logger)
sys.stderr = idq.LogFile(logger)
#================================================================================================
### read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(opts.config)

### define the gravitational wave channel/what is a glitch
gwchannel = config.get('general', 'gwchannel')
if opts.threshold:
    gwthreshold = opts.threshold
else:
    gwthreshold = config.getfloat('general', 'gw_kwsignif_thr')

### kleineWelle config
GWkwconfigpath = config.get('data_discovery', 'GWkwconfig')
GWkwconfig = idq.loadkwconfig(GWkwconfigpath)
GWkwbasename = GWkwconfig['basename']
GWgdsdir = config.get('data_discovery', 'GWgdsdir')
GWkwstride = int(float(GWkwconfig['stride']))

col_kw = {
    'tstart': 0,
    'tstop': 1,
    'tcent': 2,
    'fcent': 3,
    'uenergy': 4,
    'nenergy': 5,
    'npix': 6,
    'signif': 7,
    }
col = col_kw.copy()
### get segments and write to xml doc
gps_start_time = opts.gpsstart
gps_end_time = opts.gpsstop
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
newsegs = []
print(opts.min_seg_length)
for i in range(len(scisegs)):
    if ((scisegs[i][1] - scisegs[i][0]) >= opts.min_seg_length):
        newsegs.append(scisegs[i])
scisegs = newsegs
print("Segs: " + str(scisegs))

trigger_dict = idq.retrieve_kwtrigs(GWgdsdir, GWkwbasename, opts.gpsstart, opts.gpsstop-opts.gpsstart, GWkwstride, sleep=0, ntrials=1, logger=logger, segments=scisegs) ### go find GW kwtrgs
if gwchannel not in trigger_dict:
    trigger_dict[gwchannel] = []
print(len(trigger_dict[gwchannel]))
# apply significance threshold to the triggers from the main channel
trigger_dict.apply_signif_threshold(channels=[gwchannel], threshold=gwthreshold)
file = open("output.txt", 'w')
tag = opts.tag
print(len(trigger_dict[gwchannel]))
total_glitches = []
for seg in scisegs:
    ratio_glitches = []
    absolute_glitches=[]
    for trig in trigger_dict[gwchannel]:
        time = trig[col['tcent']]
        if (seg[0] <= time) and (time <= seg[1]):
            tsince = time - seg[0]
            absolute_glitches.append(tsince)
            tratio = tsince/(seg[1]-seg[0])
            ratio_glitches.append(tratio)
            total_glitches.append(tratio)
    file.write("GPS Start Time: %s" % (seg[0],)+ " GPS Stop Time: %s" % (seg[1],)+ " Duration: %s" % (seg[1]-seg[0],)+ " Number of Glitches: %s" % (len(absolute_glitches),))
    if len(ratio_glitches) > 0:
        plt.hist(ratio_glitches, range = (0, 1), bins = 25)
        plt.ylabel("Glitches")
        plt.xlabel("Time Normalized Over Segment Length")
        plt.title("Normalized Glitch Distribution for Segment %s-%s" % (seg[0], seg[1]))
        plt.savefig("/home/jessica.hyde/public_html/%sNormalizedSegment%s-%s.png" % (tag, seg[0], seg[1]))
        plt.close()
    if len(absolute_glitches) > 0:
        plt.hist(absolute_glitches, range = (0, seg[1]-seg[0]), bins = 25)
        plt.ylabel("Glitches")
        plt.xlabel("Raw Time Since Start of Lock")
        plt.title("Raw Glitch Distribution for Segment %s-%s" % (seg[0], seg[1]))
        plt.savefig("/home/jessica.hyde/public_html/%sRawSegment%s-%s.png" % (tag, seg[0], seg[1]))
        plt.close()
        for time in [60, 120, 300, 600, 1800]:
            plt.hist(absolute_glitches, range = (0,time))
            plt.ylabel("Glitches")
            plt.xlabel("Raw Time Since Start of Lock")
            plt.title("Glitch Distribution for First %s Seconds of Segment %s-%s" % (time, seg[0], seg[1]))
            plt.savefig("/home/jessica.hyde/public_html/%sFirst%sSegment%s-%s.png" % (tag, time, seg[0], seg[1]))
            plt.close()
            newglitchtimes = [time - seg[1]-seg[0]-glitch for glitch in absolute_glitches]
            plt.hist(newglitchtimes, range = (0,time))
            plt.ylabel("Glitches")
            plt.xlabel("Raw Time Since %s Seconds Before Lock End" % (time,))
            plt.title("Glitch Distribution for Last %s Seconds of Segment %s-%s" % (time, seg[0], seg[1]))
            plt.savefig("/home/jessica.hyde/public_html/%sLast%sSegment%s-%s.png" % (tag, time, seg[0], seg[1]))
            plt.close()
if len(total_glitches) > 0:
    plt.hist(total_glitches, range=(0,1), bins = 25)
    plt.ylabel("Glitches")
    plt.xlabel("Time Normalized Over Segment Length")
    plt.title("Normalized Glitch Distribution for All er7 Segments")
    plt.savefig("/home/jessica.hyde/public_html/%sTotalSegments.png" % (tag))
    plt.close()
file.close()
