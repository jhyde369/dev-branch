import subprocess as sp 

def compute_time_normalized(t, parsed_segs):
""" """
    for seg in parsed_segs:
        if (t >= seg[0] and t < seg[1]):
            return((t-seg[0])/(seg[1]-seg[0]))
    else:
        #print("not in range" + str(t))
        return(-1)

def compute_time_since_locked(t, parsed_segs):
""" """
    for seg in parsed_segs:
        if (t >= seg[0] and t < seg[1]):
            return(t-seg[0])
    else:
        #print("not in range" + str(t))
        return(-1)

def compute_time_until_end(t, parsed_segs):
""" """
    for seg in parsed_segs:
        if (t >= seg[0] and t < seg[1]):
            return(seg[1]-t)
    else:
        #print("not in range" + str(t))
        return(-1)

def parse_segments(segs_doc):
""" """
    input_file = open(segs_doc, "r")
    lines = input_file.readlines()
    segs = []
    for line in lines:
        start, end = line.split(",")
        segs.append([float(start), float(end)])
    input_file.close()
    return(segs)

def get_segments(gpsstart, gpsstop):    
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

