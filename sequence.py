import sys
import os


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class InputFileError(Error):
    """Excpetion raise for errors in input files."""
    pass


class Sequence(object):


    def __init__(self, fasta_name):
        """Initialize a Sequence object with the fasta file.

        Attributes:
            fa_name: name of reference genome fasta file (string)
            fai_name: name of reference genome fasta index file (string)
            fai_data: fai data: chrom,length,offset,nbase/nchar (dict)

        Raises:
            InputFileError: error in file extension of file opening

        """
        if (not os.path.exists(fasta_name) or
            not os.path.isfile(fasta_name)):
            raise InputFileError('can not read input fasta "{}"\n'
                                 .format(fasta_name))
        if not fasta_name.endswith(('.fa', '.fasta')):
            raise InputFileError('input fasta without .fa/.fasta ({})\n'
                                 .format(fasta_name))
        fasta_index_name = fasta_name + '.fai'
        if (not os.path.exists(fasta_index_name) or
            not os.path.isfile(fasta_index_name)):
            raise InputFileError('can not read input fasta index "{}"\n'
                                 .format(fasta_index_name))
        self.fa_name = fasta_name
        self.fai_name = fasta_index_name
        self.fai_data, self.chrom_list = self._load_index()


    def __str__(self):
        """Print the fasta name and fai index table of a sequence object."""

        out_str = []
        out_str.append('fasta_name: {}'.format(self.fa_name))
        out_str.append('contig\tlength\toffset\tline_nbase\tline_nchar')
        for chrom in self.chrom_list:
            out_str.append('{}\t{}\t{}\t{}\t{}'
                           .format(chrom,
                                   self.fai_data[chrom]['chrom_len'],
                                   self.fai_data[chrom]['byte_offset'],
                                   self.fai_data[chrom]['line_nbase'],
                                   self.fai_data[chrom]['line_nchar']))
        return '\n'.join(out_str)


    def _load_index(self):
        """Load a fasta index file.

        Returns:
            fai_data: a dictionary with index data
                {
                    chrom (string):
                    {
                        'chrom_len': len (int),
                        'byte_offset': offset (int),
                        'line_nbase': nbase (int),
                        'line_nchar': nchar (int),
                    },
                }
            chrom_list: a sorted list of chromosomes/contigs (for ordered
                printing)

        """
        fai_data = {}
        chrom_list = []
        with open(self.fai_name) as fai_file:
            for line in fai_file:
                chrom, chrom_len, byte_offset, line_nbase, line_nchar \
                    = line.rstrip().split('\t')
                fai_data[chrom] = {
                    'chrom_len': int(chrom_len),
                    'byte_offset': int(byte_offset),
                    'line_nbase': int(line_nbase),
                    'line_nchar': int(line_nchar),
                }
                chrom_list.append(chrom)
        print(':: loaded {} contigs from {}'
              .format(len(fai_data), self.fai_name), file=sys.stderr)
        return fai_data, chrom_list


    def query_region(self, qregion):
        """Query a sequence object with a region object.

        Args:
            region: region object initialized with region.Region(...)

        Returns:
            seq_str: a string of sequence in the query region.
            return empty string if query chrom is not found in the fai index

        """
        seq_str = ''
        if qregion.chrom not in self.fai_data:
            print('*Warning* chrom "{}" is not found in ref. seq.'
                  .format(qregion.chrom), file=sys.stderr)
            return seq_str

        chrom_offset = self.fai_data[qregion.chrom]['byte_offset']
        # query_pstart query_pend are 0-based inclusive
        if (qregion.length == -1 or
            qregion.pend > self.fai_data[qregion.chrom]['chrom_len']):
            query_start = 0
            query_end = self.fai_data[qregion.chrom]['chrom_len'] - 1
        else:
            query_start = qregion.pstart - 1
            query_end = qregion.pend - 1

        line_ndiff = (self.fai_data[qregion.chrom]['line_nchar'] -
                      self.fai_data[qregion.chrom]['line_nbase'])
        line_nbase = self.fai_data[qregion.chrom]['line_nbase']
        offset_start = (chrom_offset + query_start + line_ndiff *
                        (query_start // line_nbase))
        offset_end = (chrom_offset + query_end + line_ndiff *
                      (query_end // line_nbase))
        read_nbyte = offset_end - offset_start + 1

        with open(self.fa_name) as file_fa:
            file_fa.seek(offset_start)
            seq_str = file_fa.read(read_nbyte).replace('\n', '')

        return seq_str
