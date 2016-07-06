def compute_time_since_locked(t, parsed_segs):
    for seg in parsed_segs:
        if (t >= seg[0] and t < seg[1]):
            return(t-seg[0])
    else:
        return(None)

def parse_segments(segs_doc):
    input_file = open(segs_doc, "r")
    lines = input_file.readlines()
    segs = []
    for line in lines:
        start, end = line.split(",")
        segs.append((float(start), float(end)))
    input_file.close()
    return(segs)
    
