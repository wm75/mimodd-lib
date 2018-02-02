import hashlib

from itertools import groupby


class FastaParseError (ValueError):
    pass


RECORD_SEP = '>'
UNSAFE_ID_CHARS = '<>[]*;=,'
_UNSAFE_CHARSET = set(UNSAFE_ID_CHARS)


def is_safe_char (c):
    return 32 < ord(c) < 127 and c not in _UNSAFE_CHARSET


def is_safe_id (identifier):
    return all(is_safe_char(c) for c in identifier)
    

def get_sanitized_id (identifier, replacement_str = None):
    amended_chars = []
    for c in identifier:
        if is_safe_char(c):
            amended_chars.append(c)
        elif replacement_str:
            # a fixed replacement string has been specified
            # => use it to replace any invalid character
            amended_chars.append(replacement_str)
        else:
            # use the capitalized percent encoding of the UTF-8
            # representation of the invalid character
            utf_bytes=c.encode('utf-8')
            for byte_val in utf_bytes:
                amended_chars.append('%')
                amended_chars.append(hex(byte_val)[2:].upper().zfill(2))
    return ''.join(amended_chars)


class FastaReader (object):
    SEP = RECORD_SEP
    
    def __init__ (self, iterable):
        self.i = self._group_on_separator(iterable, self.SEP)
        
    def __iter__ (self):
        first_header, first_seq_iter = next(self.i)
        if first_header is None:
            raise FastaParseError(
                'Input does not seem to be in fasta format '
                '(expected "{0}" as first character).'
                .format(self.SEP)
                )
        yield (
            first_header,
            self._get_parsed_seqlines(first_header, first_seq_iter)
            )
        for header, seq_iter in self.i:
            yield (
                header,
                self._get_parsed_seqlines(header, seq_iter)
                )
    
    def identifiers (self):
        for header, seq_iter in self:
            yield header
            
    def sequences (self):
        for header, seq_iter in self:
            yield header, ''.join(seq_iter)
        
    def seqlens (self):
        for header, seq_iter in self:
            seqlen = 0
            for seqline in seq_iter:
                seqlen += len(seqline)
            yield header, seqlen
            
    def md5sums (self, encoding=None):
        if encoding is None:
            encoding = getattr(self.i, 'encoding', 'ascii')
        for header, seq_iter in self:
            md5 = hashlib.md5()
            seqline = None
            for seqline in seq_iter:
                md5.update(
                    seqline.upper().encode(encoding)
                    )
            if seqline is None:
                yield header, None
            else:
                yield header, md5.hexdigest()

    def describe_records (self, encoding=None):
        if encoding is None:
            encoding = getattr(self.i, 'encoding', 'ascii')
        for header, seq_iter in self:
            seqlen = 0
            md5 = hashlib.md5()
            for seqline in seq_iter:
                seqlen += len(seqline)
                md5.update(
                    seqline.upper().encode(encoding)
                    )
            if seqlen == 0:
                yield header, seqlen, None
            else:
                yield header, seqlen, md5.hexdigest()
    
    def _get_parsed_seqlines (self, header, seq_iter):
        try:
            for n, seqline in enumerate(seq_iter, 1):
                yield self._parse_sequence_line(seqline)
        except FastaParseError as e:
            # add the record identifier and the line number within the record
            # to the error message and the offending character in the
            # exception args.
            e.args = e.args[0], e.args[1], header, n
            raise e

    def _parse_sequence_line (self, raw_seq):
        return raw_seq.strip().replace(' ', '')

    def _group_on_separator (self, iterable, separator):
        """Returns start-tagged records - tuple(header, content) - from an iterable.

        Looks for the record separator at the start of each element.
        If the speparator is found, the rest of that element is interpreted
        as the record header. All elements between the header and the next separator
        element are considered the record content which is returned as an iterator.
        Missing content is reported as [], a missing header for the first record
        as None.
        BEWARE: in the generated tuple, the content iterator is only meaningful until
        the next iteration over the generator, and MUST NOT be used afterwards.
        Intended as a building block for higher level classes."""
        
        sep_len = len(separator)
        header_tail = None
        for is_header, item in groupby(
            iterable, lambda line: line[:sep_len] == separator
            ):
            if is_header:
                header_tails = [h[sep_len:].strip() for h in item]
                for naked_header in header_tails[:-1]:
                    yield (naked_header, iter([]))
                header_tail = header_tails[-1]
            else:
                yield (header_tail, item)


class FastaWithAlphabetReader (FastaReader):
    def __init__ (self, iterable, alphabet):
        super().__init__(iterable)
        self.alphabet = alphabet

    def _parse_sequence_line (self, raw_seq):
        seq = super()._parse_sequence_line(raw_seq)
        for c in seq:
            if c not in self.alphabet:
                raise FastaParseError(
                    'Invalid letter in sequence.', c
                    )
        return seq


class FastaNucleotideReader (FastaWithAlphabetReader):
    def __init__ (self, iterable):
        valid_IUPAC_symbols = 'ACGTNKSYMWRBDHV'
        valid_IUPAC_symbols += valid_IUPAC_symbols.lower()
        super().__init__(iterable, set(valid_IUPAC_symbols))
