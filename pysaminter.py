import pysam

from . import seqtransform

from .seqreads import SeqReadFacade


class UnAlignedRead (pysam.AlignedRead):
    """Enables convenient initialization of an unmapped read."""
    
    def __init__ (self):
        super().__init__()
        self.tid = -1
        self.pos = -1
        self.is_unmapped = True
        self.mapq = 255
        self.cigar = None
        self.tlen = 0
        self.rnext = self.pnext = -1


def new_unaligned_pair ():
    """Creates two UnAlignedRead instances and sets mate-specific flags."""
    
    r1 = UnAlignedRead()
    r2 = UnAlignedRead()
    r1.flag = r1.flag | 73  # set is_paired, is_read1 & mate_is_unmapped bits
    r2.flag = r2.flag | 137 # set is_paired, is_read2 & mate_is_unmapped bits

    return PysamRead(r1), PysamRead(r2)
    
    
class PysamRead (SeqReadFacade):
    def __init__ (self, read_object=None, is_dna=True):
        if read_object is None:
            read_object = UnAlignedRead()
        super().__init__(read_object, is_dna)
    
    @property
    def sequence (self):
        return self.read.seq

    @property
    def quality (self):
        return self.read.qual

    @property
    def full_title (self):
        return self.read.qname

    @property
    def read_data (self):
        return self.read.qname, self.read.seq, self.read.qual

    @property
    def flag (self):
        return self.read.flag

    @property
    def tags (self):
        return self.read.tags

    def update_read(read_data):
        self.read.qname, self.read.seq, self.read.qual = read_data

    def reverse (self):
        self.read.seq, self.read.qual = (
            self.read.seq[::-1], self.read.qual[::-1]
            )

    def complement(self):
        self.read.seq = seqtransform.complement(
            self.read.seq, is_dna=self.is_dna
            )

    def reverse_complement(self):
        self.read.seq, self.read.qual = (
            seqtransform.complement(self.read.seq, is_dna=self.is_dna),
            self.read.qual[::-1]
            )
