# TO DO:
# it is inconvenient to have to call iter on ReadSplitter objects first
# (it's also kind of surprising given that they wrap iterators themselves)
#
# update documentation


"""
Provide classes and functions for processing reads of nucleic acid sequences.

A sequence read, in this context, is composed of an identifier, the actual
sequence determined in an experiment and a string of quality scores, one for
every base call.
"""


class SeqReadFacade (object):
    def __init__ (self, read_object, is_dna=True):
        self.read = read_object
        self.is_dna = is_dna
            
    def _parse_identifier (self):
        parts = self.full_title.split(None, 1)
        if len(parts) == 1:
            return parts[0], type(parts[0])()
        else:
            return parts[0], parts[1]

    @property
    def identifier (self):
        return self._parse_identifier()[0]

    @property
    def description (self):
        return self._parse_identifier()[1]
    
    @property
    def sequence (self):
        raise NotImplementedError

    @property
    def quality (self):
        raise NotImplementedError

    @property
    def full_title (self):
        raise NotImplementedError

    @property
    def read_data (self):
        return self.full_title, self.sequence, self.quality

    @property
    def flag (self):
        raise NotImplementedError

    @property
    def tags (self):
        raise NotImplementedError

    @property
    def rg_id (self):
        """Retrieve the read group id from the read object's tags.

        This default implementation assumes that tags is structured like in
        a pysam read object.
        """
        this_rg_id = None
        for tag in self.read.tags:
            if tag[0] == 'RG':
                this_rg_id = tag[1]
                break
        return this_rg_id

    def reverse (self):
        raise NotImplementedError

    def complement (self):
        raise NotImplementedError
    
    def reverse_complement (self):
        self.reverse()
        self.complement()
        
    def __len__ (self):
        return len(self.sequence)

    def __str__ (self):
        return '\n'.join(
            str(self.full_title), str(self.sequence), str(self.quality)
            )


class ReadSplitter (object):
    """Supprt grouped retrieval of reads from an iterable.

    This class wraps group_reads_by_qname and adds optional
    sanity checks for read to read group associations.
    """
    
    def __init__ (self, src, header=None, default_rg=None):
        """Initialize a ReadSplitter instance.

        src is the iterable to retrieve reads from.
        header is an optional header dictionary that, under an 'RG' key,
        lists read groups expected to be associated with the reads in iterable.
        When header is provided and lists any read groups under its 'RG' key,
        then any read found associated with a read group not declared in
        header causes an exception during iteration.
        If header is provided, but does not list any read groups under its 'RG'
        key, then all reads iterated over must be associated with the same
        read group as the first read or an exception will be raised.
        If no header is provided, then no checks on read groups are performed.
        If a default_rg is specified, it will be used instead of None as the
        read group for all reads that do not declare their read group
        explicitly. This default_rg will be substituted for missing read
        groups before any read group sanity checks.
        """
        if header and 'RG' not in header:
            raise ValueError('header needs to provide a "RG" key')
        if header is None:
            self.rg_set = None
        else:
            self.rg_set = {rg['ID'] for rg in header['RG']}
        if self.rg_set and default_rg is not None:
            raise ValueError(
                'cannot use default_rg when read groups are stated explicitly '
                'in header'
                )
        self.src = iter(src)
        self.header = header
        self.default_rg = default_rg

    def __iter__ (self):
        """
        Wrap group_reads_by_qname to provide iteration with read group checks.
        """
        it = group_reads_by_qname(self.src)
        try:
            rg_id, grouped_reads = next(it)
        except StopIteration:
            raise RuntimeError('No reads in file. Aborting.')
        if not self.rg_set:
            if self.rg_set is not None:
                rg_id = rg_id or self.default_rg
                self.rg_set = {rg_id}
        elif rg_id not in rg_set:
                # fixed RG info and not matching
                raise
        yield rg_id, grouped_reads
        for rg_id, grouped_reads in it:
            rg_id = rg_id or self.default_rg
            if self.rg_set and rg_id not in self.rg_set:
                if rg_id is None:
                    raise FormatParseError(
                        'Missing "RG" tag for read from a SAM/BAM file with multiple '
                        'read groups.',
                        help='The input file defines multiple read groups in its '
                        'header, but not all reads in the file body specify which '
                        'read group they belong to.'
                        )
                raise FormatParseError(
                    'Unknown "RG" tag "{token}" for read.',
                    token=this_rg_id,
                    help='A read in the file body claims it is belonging to a '
                    'read group that is not defined in the file header.'
                    )
            yield rg_id, grouped_reads


def group_reads_by_qname (src):
    """
    Yield primary reads from iterable src grouped by read group and identifier.

    Yields tuples of (read group, [read, ...]), in which all reads have the
    same identifier, i.e., are segments obtained from the same template
    during sequencing.
    Inspects the flag attribute of all read objects to ensure that it yields
    each unique read only once in the form of its primary alignment.
    Since this function yields a group of reads when it encounters the first
    read not belonging to the group, it requires reads sorted by read
    identifiers to be able to return all segments of all reads as single
    groups.
    """
    try:
        read1 = next(src)
    except StopIteration:
        return
    ret_type = type(read1)
    while read1:
        if read1.flag & 0x900:
            # skip reads until a primary read is found
            try:
                read1 = next(src)
                continue
            except StopIteration:
                break

        this_rg_id = read1.rg_id
        this_read_id = read1.identifier
        ret_reads = [ret_type(read1.read)]
        while True:
            try:
                readn = next(src)
            except StopIteration:
                readn = None
                break
            if this_read_id != readn.identifier:
                break
            if this_rg_id != readn.rg_id:
                break
            if readn.flag & 0x900:
                # still same read name, but not a primary alignment
                # skip and try the next read
                continue
            ret_reads.append(ret_type(readn.read))
        yield (this_rg_id, ret_reads)
        read1 = readn    
