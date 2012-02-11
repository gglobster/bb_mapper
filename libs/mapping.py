from config import segtype
from array_tetris import offset_q2r_coords
from drawing import pairwise_draw
from loaders import load_genbank
import numpy as np

def map_pw_aln(ref, genome, segs_dir, maps_dir):
    """Generate map of construct aligned to reference."""
    # set inputs and outputs
    seg_file = segs_dir+"/"+genome.name+"_"+ref.name+"_segs.txt"
    map_file = maps_dir+genome.name+"_vs_"+ref.name+".pdf"
    # start mapping
    try: open(genome.gbk)
    except IOError:
        print "WARNING: No scaffold construct to map"
    else:
        try:
            # load segments TODO: add idp-based clumping
            segdata = np.loadtxt(seg_file, skiprows=1, dtype=segtype)
        except IOError:
                msg = "\nERROR: could not load segments data"
                print msg
        except StopIteration:
                msg = "\nERROR: could not make map"
                print msg
        else:
            # offset coordinates where desired
            try:
                g_offset = genome.offset
                if g_offset[0] != 0 or g_offset[1] != 0:
                    q_len = len(load_genbank(genome.gbk).seq)
                    segdata = offset_q2r_coords(segdata, q_len, g_offset)
                # determine whether to flip the query sequence (neg. offset)
                if g_offset[1] < 0:
                    q_invert = True
                else:
                    q_invert = False
            except KeyError:
                g_offset = (0,0)
                q_invert = False
            # generate graphical map
            pairwise_draw(ref.name, genome.name, ref.gbk, genome.gbk, segdata,
                          map_file, q_invert, g_offset, 'dual', 'dual', 'm',
                          'fct', 'fct')
