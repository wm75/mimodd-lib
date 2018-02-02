# TO DO:
# it is inconvenient to have to call iter on FastqReader and ReadSplitter
# objects first (it's also kind of surprising given that they wrap iterators
# themselves)
#
# update all names and documentation to reflect changes


"""
Provide classes and functions for processing reads of nucleic acid sequences.

A sequence read, in this context, is composed of an identifier, the actual
sequence determined in an experiment and a string of quality scores, one for
every base call.

This information can be stored in string or object form and the module
supports conversion between both.

Conversion between serialized and OO representations of sequence reads
----------------------------------------------------------------------

FastqReader:
Use this class to parse a stream in fastq format into records or read objects.

FastqWriter:
Use this class to serialize read records or objects to fastq format.
"""

class SeqReadFacade (object):
    def __init__ (self, read_object):
        self.read = read_object

    def _parse_identifier (self):
        parts = self.full_title.split(None, 1)
        if len(parts) == 1:
            return parts[0], ''
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
        return '\n'.join(self.full_title, self.sequence, self.quality)


class SimpleSeqRead (SeqReadFacade):
    def __init__ (self, read_object=None):
        if read_object is None:
            read_object = (None, None, None)
        super().__init__(read_object)
    
    @property
    def sequence (self):
        return self.read[1]

    @property
    def quality (self):
        return self.read[2]

    @property
    def full_title (self):
        return self.read[0]

    @property
    def read_data (self):
        return self.read[0], self.read[1], self.read[2]
    
    @property
    def flag (self):
        return 4 # simple sequencing reads are assumed to be unmapped
    
    @property
    def rg_id (self):
        return None # a simple sequencing read does not belong to a read group

    def reverse (self):
        self.read = (self.read[0], self.read[1][::-1], self.read[2][::-1])

    def complement (self):
        raise NotImplementedError

    
class FastqReader (object):
    """Parse a line-based stream in fastq format into records or read objects.

    Provide iteration over sequence reads encoded in fastq format in a text or
    binary stream and yield reads as simple records of
    (identifier, sequence, quality_scores) or as read objects with the record
    contents stored in attributes.
    Yields read content of the same type as the input, i.e., bytes or text
    strings depending on the input stream.
    Fastq is a notoriously ill-specified format so this class offers very
    lenient, but fast parsing, which should be compatible with most format
    variations.
    """
    
    def __init__ (self, src, read_object=None):
        """Initialize a FastqReader instance.

        The instance will consume the iterable src, which is expected to
        yield elements of type bytes or str.
        Note that bytes are parsed faster than strings because fastq format
        stores quality scores as single bytes so these can be retrieved from
        a quality score bytes string directly without decoding.
        If a read_object is provided, it must support setting qname, seq and
        qual attributes. If read_object is an instance of a class, iteration
        over src will yield read_object repeatedly with these attributes set
        to the current record identifier, sequence and quality scores,
        respectively. If read_object is a class, it will be instantiated and
        the instance be used to communicate record contents.
        If no read_object is provided results will be generated in the form
        of tuples of (identifier, sequence, quality scores).
        The primary usecase for providing a record_object is to retrieve
        records as pysam AlignedRead objects for further processing with
        pysam.
        """
        self.src = iter(src)
        if read_object is None:
            self.robject = SimpleSeqRead()
        else:
            self.robject = read_object
        self._it = self._read_fastq_records()
        try:
            # prime the fastq record generator
            # This will yield a throw-away None, but gives
            # the generator a chance to look at the first element in src
            # to see if it is of type bytes or str and set the instance's
            # is_bytes_source attribute accordingly.
            next(self._it)
        except StopIteration:
            # src is empty.
            # We suppress the exception here because we only want to raise it
            # during iteration.
            pass
            
    def __iter__ (self):
        """Provide record or read object based iteration."""
        for record in self._it:
            self.robject.read = record
            yield self.robject

    def _read_fastq_records (self):
        """Parse a fastq format fast and robustly.

        Deal with multiline sequences and quality scores.
        Allow arbitrary numbers of empty lines anywhere in the file.
        Perform (only) the following file format checks while parsing:
        - each title line MUST start with an @ symbol
        - each record MUST have a sequence
        - sequence and quality score of each record MUST be of equal length
        No alphabet checks are done on sequences or quality scores.

        Yield tuples of (identifier, sequence, quality scores), where
        identifier includes the leading @ and none of the elements are
        processed beyond removal of trailing whitespace.
        """
        try:
            title = next(self.src)
        except StopIteration:
            # src is an empty iterator.
            # We do not know whether it was meant to yield bytes or str
            # elements.
            self.is_bytes_source = None
            return
        if isinstance(title, bytes):
            title_token = ord('@')
            sep_token = ord('+')
            glue = b''
            self.is_bytes_source = True
        elif isinstance(title, str):
            title_token = '@'
            sep_token = '+'
            glue = ''
            self.is_bytes_source = False
        else:
            raise TypeError(
                'An iterator over bytes or str elements is required. '
                'Found element of type {0}.'
                .format(type(title).__name__)
                )
        # With the source format figured out, pause here to give the
        # instance a chance to finish configuring itself.
        yield None
        
        while True:
            if not title[0] == title_token:
                # allow empty lines between records
                if not title.rstrip(): continue
                raise ValueError(
                    'Invalid format: Title line not starting with @',
                    title
                    )
            title = title[1:].rstrip()
            line_tmp = []
            try:
                # from here on any StopIteration means an incomplete record
                while True:
                    currentLine = next(self.src)
                    # we are accepting arbitrary numbers of empty lines anywhere
                    if currentLine[0] == sep_token: break
                    line_tmp.append(currentLine.rstrip())
                seq = glue.join(line_tmp)
                seqlen = len(seq)
                if seqlen == 0:
                    raise ValueError(
                        'Invalid format: Record without sequence',
                        title
                        )
                quallen = 0
                line_tmp = []
                while seqlen > quallen:
                    currentLine = next(self.src).rstrip()
                    # again we are accepting any number of empty lines
                    line_tmp.append(currentLine)
                    quallen += len(currentLine)
            except StopIteration:
                raise ValueError(
                    'Invalid format: Last record is incomplete',
                    title
                    ) from None
            if seqlen < quallen:
                raise ValueError(
                    'Invalid format: Record with inconsistent lengths of '
                    'sequence and quality score',
                    title
                    )
            qual = glue.join(line_tmp)
            yield title, seq, qual

            try:
                title = next(self.src)
            except StopIteration:
                # no more records to parse
                return


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
        self.src = src
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
