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
file = open("output.txt", 'w')
tag = opts.tag
print(len(trigger_dict[gwchannel]))
total_glitches = []
for seg in scisegs:
    absolute_glitches=[]
    for trig in trigger_dict[gwchannel]:
        time = trig[col['tcent']]
        signif = trig[col['signif']]
        if (seg[0] <= time) and (time <= seg[1]):
            total_glitches.append(signif)
            absolute_glitches.append(signif)
    file.write("GPS Start Time: %s" % (seg[0],)+ " GPS Stop Time: %s" % (seg[1],)+ " Duration: %s" % (seg[1]-seg[0],)+ " Number of Glitches: %s" % (len(absolute_glitches),) + "\n")
    if len(absolute_glitches) > 0:
        plt.hist(absolute_glitches, range = (50, 0), bins = 30)
        plt.ylabel("Glitches")
        plt.xlabel("Significance")
        plt.title("Glitch Counts Versus Significance for Segment %s-%s" % (seg[0], seg[1]))
        plt.savefig("/home/jessica.hyde/public_html/%sSignificanceSegment%s-%s.png" % (tag, seg[0], seg[1]))
        plt.close()
    if len(absolute_glitches) > 0:
        plt.hist(absolute_glitches, range = (50, 0), bins = 30, cumulative=True)
        plt.ylabel("Cumulative Glitches")
        plt.xlabel("Significance")
        plt.title("Cumulative Glitch Counts Versus Significance for Segment %s-%s" % (seg[0], seg[1]))
        plt.savefig("/home/jessica.hyde/public_html/%sCumulativeSignificanceSegment%s-%s.png" % (tag, seg[0], seg[1]))
        plt.close()
if len(total_glitches) > 0:
    plt.hist(total_glitches, range=(100, 0), bins = 50)
    plt.ylabel("Glitches")
    plt.xlabel("Significance")
    plt.title("Glitch Counts Versus Significance for All er7 Segments")
    plt.savefig("/home/jessica.hyde/public_html/%sSignifTotalSegments.png" % (tag))
    plt.close()
if len(total_glitches) > 0:
    plt.hist(total_glitches, range=(100, 0), bins = 50, cumulative=True)
    plt.ylabel("Cumulative Glitches")
    plt.xlabel("Significance")
    plt.title("Cumulative Glitch Counts Versus Significance for All er7 Segments")
    plt.savefig("/home/jessica.hyde/public_html/%sCumulativeSignifTotalSegments.png" % (tag))
    plt.close()
file.close()
