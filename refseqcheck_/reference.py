import pysam


class Reference(object):
    """"Class for reference genome"""


    def __init__(self, filename):
        """"Constructor of Reference class"""

        self.fastafile = pysam.Fastafile(filename)


    def get_sequence(self, chrom, start, end):
        """"Retrieving the sequence of a genomic region"""

        # Checking if chromosome name exists
        goodchrom = chrom
        if not goodchrom in self.fastafile.references:
            goodchrom = 'chr' + chrom
            if not goodchrom in self.fastafile.references:
                if chrom == 'MT':
                    goodchrom = 'chrM'
                    if not goodchrom in self.fastafile.references:
                        return None
                else:
                    return None

        # Fetching data from reference genome
        if end <= start:
            return ''
        if start < 0:
            start = 0

        if pysam.__version__ in ['0.7.7', '0.7.8', '0.8.0']:
            last = self.fastafile.getReferenceLength(goodchrom)
        else:
            last = self.fastafile.get_reference_length(goodchrom)

        if end > last: end = last
        seq = self.fastafile.fetch(goodchrom, start, end)
        return seq.upper()
