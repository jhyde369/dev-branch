#================================
usage="a string giving an example of how to use the script: eg: script.py [--options] args"
description="a string describing what the script does"

import os
import subprocess as sp

from optparse import OptionParser
import module
#================================

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-t", "--gps-time", default=[], type="float", action="append", help="the gps time for which we calculate the time since locked. You can specify multiple times by repeating this argument.")

opts, args = parser.parse_args()

#================================

#os.system("ligolw_print -c start_time -c end_time -t segment sample_segments.xml")
#os.system("ligolw_segment_query_dqsegdb -t https://dqsegdb5.phy.syr.edu -q -a L1:DMT-ANALYSIS_READY:1 -s 1117224000 -e 1117832400")

print opts
print len(args)


### get segments and write to xml doc
proc = sp.Popen(["ligolw_segment_query_dqsegdb", "-t", "https://dqsegdb5.phy.syr.edu", "-q", "-a", "L1:DMT-ANALYSIS_READY:1", "-s", "1117224000", "-e", "1117832400", "-o", "sample_segments.xml"], stdout = sp.PIPE)
ans = proc.communicate()[0]
### get start and end times from xml doc, write to output doc
output_file = open("output.txt", "w")
proc = sp.Popen(["ligolw_print", "-c" "start_time", "-c" "end_time", "-t" "segment", "sample_segments.xml"], stdout=output_file)
ans = proc.communicate()[0]
output_file.close()

### 
segs = module.parse_segments("output.txt")
print module.compute_time_since_locked(segs[0][0] + 100, segs)
