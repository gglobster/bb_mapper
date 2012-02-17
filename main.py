from sys import argv
from libs.writers import ensure_dir
from libs.scripts import align_pairwise, align_multi

print "\n", \
      "##################################################\n", \
      "### BB-Mapper v. 0.1  (prev. PlasmiG)          ###\n", \
      "### Copyright 2012 Geraldine A. Van der Auwera ###\n", \
      "##################################################\n"

#from sets.serratia import genomes
from sets.thuricins import genomes

run = argv[1]

# prep directories
data_dir = "data/"+argv[2]+"/"
dirs = {
    'seqfiles': data_dir+"seqfiles/",
    'mauve': data_dir+"mauve/",
    'aln_segs': data_dir+"seg_aln/",
    'maps': data_dir+"maps/"
}
for item in dirs.values():
    ensure_dir([item])

# determine mode (new alignment or redrawing)
new_align = True
if len(argv) > 3:
    mode = argv[3]
    if mode == 'redraw':
        new_align = False

# one-pair or multiple?
if len(genomes) < 2:
    raise Exception("ERROR: Not enough genomes to align!")
elif len(genomes) == 2:
    align_pairwise(genomes, new_align, dirs, run)
else:
    align_multi(genomes, new_align, dirs, run)



