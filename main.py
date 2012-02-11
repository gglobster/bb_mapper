import os, re
from sys import argv
from libs.aligning import align_mauve, iter_align
from libs.parsing import mauver_load2_k0
from libs.loaders import load_genbank, from_dir
from libs.writers import ensure_dir
from libs.array_tetris import chop_rows
from libs.mapping import map_pw_aln
from libs.classes import Noodle
from config import max_size, chop_mode


print "\n", \
      "##################################################\n", \
      "### BB-Mapper v. 0.1  (prev. PlasmiG)          ###\n", \
      "### Copyright 2012 Geraldine A. Van der Auwera ###\n", \
      "##################################################\n", \


# one-pair or multiple?
# start with one-pair, will add multi abilities later

data_dir = "data/"+argv[1]+"/"

# prep directories
seqfile_dir = data_dir+"seqfiles/"
mauve_dir = data_dir+"mauve/"
aln_segs_dir = data_dir+"seg_aln/"
maps_dir = data_dir+"maps/"

ensure_dir([seqfile_dir, mauve_dir, aln_segs_dir, maps_dir])

# load inputs
ref = Noodle(argv[2], seqfile_dir+argv[3], argv[4])
genome = Noodle(argv[5], seqfile_dir+argv[6], argv[7])

new_align = True
if len(argv) > 8:
    mode = argv[8]
    if mode == 'redraw':
        new_align = False

# set outputs
mauve_outfile = mauve_dir+genome.name+"_"+ref.name+".mauve"
segfile = aln_segs_dir+genome.name+"_"+ref.name+"_segs.txt"


if new_align:
    # prep segments file
    open(segfile, 'w').write('')

    # purge any pre-existing sslist files
    sslist_files = from_dir(seqfile_dir, re.compile(r'.*\.sslist.*'))
    for sslist in sslist_files:
        try: os.remove(seqfile_dir+sslist)
        except Exception: raise

    # do Mauve alignment
    file_list = [ref.gbk, genome.gbk]
    print "Aligning", ref.name, "and", genome.name, "..."
    align_mauve(file_list, mauve_outfile)

    try:
        # parse Mauve output (without initial clumping)
        coords = mauver_load2_k0(mauve_outfile+".backbone", 0)
        print "Segment results:", len(coords), '->',

        # chop segments that are too long
        chop_array = chop_rows(coords, max_size, chop_mode)
        print len(chop_array), 'segments <', max_size, 'bp'

        # make detailed pairwise alignments of the segments
        print "Aligning segments ..."
        ref_rec = load_genbank(ref.gbk)
        query_rec = load_genbank(genome.gbk)
        id = iter_align(chop_array, ref_rec, query_rec, aln_segs_dir, segfile)
        print "Results:", id, "% id. overall"

    except IOError:
        print "\nERROR: Mauve alignment failed"

# map of construct aligned to reference
print "Mapping....",
map_pw_aln(ref, genome, aln_segs_dir, maps_dir)
print "OK\n"