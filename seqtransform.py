IUPAC_DNA = 'AGRMHBSNWVDKYCT'
IUPAC_RNA = 'AGRMHBSNWVDKYCU'


# bytes translation tables can be used to translate str and bytes
DNA_TRANSLATION_TABLE = bytes.maketrans(
    bytes(ord(c) for c in IUPAC_DNA+IUPAC_DNA.lower()),
    bytes(ord(c) for c in IUPAC_DNA[::-1]+IUPAC_DNA[::-1].lower())
    )

RNA_TRANSLATION_TABLE = bytes.maketrans(
    bytes(ord(c) for c in IUPAC_RNA+IUPAC_RNA.lower()),
    bytes(ord(c) for c in IUPAC_RNA[::-1]+IUPAC_RNA[::-1].lower())
    )


def reverse (seq):
    return seq[::-1]


def complement (seq, is_dna=True):
    if is_dna:
        return seq.translate(DNA_TRANSLATION_TABLE)
    else:
        return seq.translate(RNA_TRANSLATION_TABLE)


def reverse_complement (seq, is_dna=True):
    if is_dna:
        return seq[::-1].translate(DNA_TRANSLATION_TABLE)
    else:
        return seq[::-1].translate(RNA_TRANSLATION_TABLE)
