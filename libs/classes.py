import sys
from libs.converters import fas2gbk, gbk2fas

class Noodle(object):
    """Object that holds info about a piece of DNA."""

    def __init__(self, name, file, format):
        self.name = name
        self.fas = None
        self.gbk = None
        self.offset = (0,0)

        if format == 'fas':
            self.fas = file
            self.gbk = fas2gbk(file)
        elif format == 'gbk':
            self.gbk = file
            self.fas = gbk2fas(file)
        else:
            print "ERROR in input format: FastA or Genbank required"
            sys.exit()

        