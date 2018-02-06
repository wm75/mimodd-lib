# TO DO:
#
# update documentation


"""
Provide classes for parsing fastq-formatted streams and for representing
sequenced reads encoded in them.

SimpleSeqRead:
A facade providing standardized and high-level access to a simple
(identifier, sequence, quality) tuple.

FastqReader:
Parse a stream in fastq-format into sequenced read objects.
"""


from . import seqtransform, seqreads


class SimpleSeqRead (seqreads.SeqReadFacade):
    def __init__ (self, read_object=None, is_dna=True):
        if read_object is None:
            read_object = (None, None, None)
        super().__init__(read_object, is_dna)

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
        self.read = (
            self.read[0],
            seqtransform.complement(self.read[1], is_dna=self.is_dna),
            self.read[2]
            )

    def reverse_complement (self):
        self.read = (
            self.read[0],
            seqtransform.reverse_complement(self.read[1], is_dna=self.is_dna),
            self.read[2][::-1]
            )


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
            self._seqread = SimpleSeqRead()
        else:
            self._seqread = read_object
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

    def __next__ (self):
        self._seqread.read = next(self._it)
        return self._seqread

    def __iter__ (self):
        """Provide record or read object based iteration."""
        return self

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
