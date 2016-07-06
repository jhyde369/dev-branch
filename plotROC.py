from laldetchar.idq import idq
from laldetchar.idq import idq_summary_plots as isp
import numpy as np


def generate_filenames(gpsstart, gpsstop, tag="B"):
    root = "/home/jessica.hyde/test/dev-branch/now/train/%s_%s/mvsc/" % (gpsstart, gpsstop)
    file = "L1_mla-%s-%s-%s.dat" % (gpsstart, gpsstop-gpsstart, tag)
    filenames = []
    for tstop in ["no_tend", "yes_tend"]:
        for tstart in ["no_tstart", "yes_tstart"]:
            filenames.append(root + "%s-%s/" % (tstart, tstop) + file)
    return(filenames)

def generate_lock_filenames(gpsstart, gpsstop, segstart, segstop, tag = ""):
    root = "/home/jessica.hyde/test/dev-branch/now/train/%s_%s/mvsc/" % (gpsstart, gpsstop)
    file = "L1_mla%s-%s-%s-B.dat" % (tag, int(segstart), int(segstop-segstart))
    filenames = []
    for tstop in ["no_tend", "yes_tend"]:
        for tstart in ["no_tstart", "yes_tstart"]:
            filenames.append(root + "%s-%s/" % (tstart, tstop) + file)
    return(filenames)

#filenames = ["/home/jessica.hyde/test/dev-branch/now/train/1117329856_1117331656/mvsc/no_tstart-no_tend/L1_mla-1117329856-1007080-B.dat", "/home/jessica.hyde/test/dev-branch/now/train/1117329856_1117331656/mvsc/yes_tstart-yes_tend/L1_mla-1117329856-1007080-B.dat"]
#filenames = ["/home/jessica.hyde/test/dev-branch/now/train/1117329856_1007080/mvsc/no_tstart-no_tend/L1_mla-1117329856-1007080-B.dat", "/home/jessica.hyde/test/dev-branch/now/train/1117329856_1007080/mvsc/yes_tstart-no_tend/L1_mla-1117329856-1007080-B.dat", "/home/jessica.hyde/test/dev-branch/now/train/1117329856_1007080/mvsc/no_tstart-yes_tend/L1_mla-1117329856-1007080-B.dat", "/home/jessica.hyde/test/dev-branch/now/train/1117329856_1007080/mvsc/yes_tstart-yes_tend/L1_mla-1117329856-1007080-B.dat"]
#filenames= ["/home/jessica.hyde/test/dev-branch/now/train/1117329856_1118336936/mvsc/no_tstart-no_tend/L1_mla-1117329856-1007080-B.dat", "/home/jessica.hyde/test/dev-branch/now/train/1117329856_1118336936/mvsc/yes_tstart-no_tend/L1_mla-1117329856-1007080-B.dat", "/home/jessica.hyde/test/dev-branch/now/train/1117329856_1118336936/mvsc/no_tstart-yes_tend/L1_mla-1117329856-1007080-B.dat", "/home/jessica.hyde/test/dev-branch/now/train/1117329856_1118336936/mvsc/yes_tstart-yes_tend/L1_mla-1117329856-1007080-B.dat"]

def plot_ROC_combinations(filename, gpsstart, gpsstop, segstart=None, segstop=None):
    if segstart:
        filenames = generate_lock_filenames(gpsstart, gpsstop, segstart, segstop)
    else:
        filenames = generate_filenames(gpsstart, gpsstop, tag="A")
    colors=['red', 'purple', 'green', 'blue']
    labels = ["original mvsc", "added time locked", "added time until end", "added both"]
    columns = ["GPS", "i", "rank"]
    print(filenames[0])
    output = idq.slim_load_datfiles([filenames[0]], columns=columns, skip_lines=0)
    for col in columns:
        output[col] = [float(ele) for ele in output[col]]
    r, c, g = idq.dat_to_rcg(output)
    color = colors[0]
    fig, ax = isp.rcg_to_rocFig(c, g, color=color, label = "original mvsc")
    print(filenames[0])
    i=1
    for dat in filenames[1:]:
        print(dat, colors[i], labels[i])
        output = idq.slim_load_datfiles([dat], columns=columns, skip_lines = 0)
        for col in columns:
            output[col] = [float(ele) for ele in output[col]]
        r, c, g = idq.dat_to_rcg( output ) ### gives cumulative numbers of cleans and glitches at different rank values.
        color = colors[i]
        fig, ax = isp.rcg_to_rocFig(c, g, color=color, label=labels[i], figax=(fig,ax)) ### makes ROC plot. Returns matplotlib figure and axis objects
        ax.legend(loc='best') ### add a legend
        i=i+1
    xvalues = range(1,1000)
    axis = [1.0/ele for ele in xvalues]
    print(axis[5])
    npaxis = np.array(axis)
    ax.plot(npaxis, npaxis, color="black")
    fig.savefig( "/home/jessica.hyde/public_html/%s" % (filename,) ) ### save the figure as "figure.png"
    isp.close(fig) ### close the figure
