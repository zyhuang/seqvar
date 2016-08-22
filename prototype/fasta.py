# -*- coding: utf-8 -*-
"""This module contains functions dealing with fasta file.

    query_region() query fasta file with a given region.

        Example:

            $ python3 fasta.py -f fasta_name -c chrom -a pstart -b pend

        This will query the fasta file "fasta_name" and return a sequence in
        region "chrom:pstart-pend", e.g.

            > 16:90000-90060
            ACCATGCCCAGCATGAGCCACTGCACCCAGATTAATTTTTGTATTTTTAGTAAAGATGAGG

"""

import sys
import os
import argparse


# =============================================================================

class Fasta(object):
    """A class for fastq file query.

    Attributes:
        fa_name (str): name of reference genome fasta file
        fai_name (str): name of reference genome fasta index file
        fai_data (dict): key = chrom, value = chrom lenth and offset in fasta
        chrom_list (list): ordered list of chromosome names

    """

    def __init__(self, fasta_name, mode='r'):
        """Initialize a Fastq object to read/write.

        Args:
            fasta_name (string): fasta file name

        """

        self.fa_name = fasta_name
        self.fai_name = fasta_name + '.fai'
        if mode == 'r':
            self.fai_data, self.chrom_list = self.read_index()
        else:
            # functions for writing a fasta file
            self.fai_data = {}
            self.chrom_list = []
            pass

# =============================================================================

    def __str__(self):
        """Print the fasta name and fai index table."""

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

# =============================================================================

    def read_index(self):
        """Read a fasta index file.

        Returns:
            fai_data: a dictionary with index data
            {
                chrom (string):
                {
                    'chrom_len': len (int), // chromosome name
                    'byte_offset': offset (int), // offset in BYTE
                    'line_nbase': nbase (int), // #base per line
                    'line_nchar': nchar (int), // #char per line (with \n)
                },
            }
            chrom_list: an ordered list of chromosomes/contigs (for printing)

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
        if not fai_data or not chrom_list:
            print('*WARNING* can not find any chromosome/contig in {}'
                  .format(self.fai_name), file=sys.stderr, flush=True)

        return fai_data, chrom_list

# =============================================================================

    def query_region(self, chrom, pstart, pend):
        """Query a sequence region in the fasta.

        Args:
            chrom (str): query chromosome name
            pstart (int): query start position (1-based inclusive, 1 if None)
            pend (int): query end position (1-based inclusive, chrom_len if
                None)

        Returns:
            region (str): query region in format chrom:pstart:pend (pstart and
                pend are 1-based inclusive)
            sequence (str): fasta sequence within the query region
            None if input error happens (unknown chrom, or illegal position
                range)
        """

        sequence = None
        if chrom not in self.fai_data:
            print('*Warning* query chrom "{}" is not found in fasta {}'
                  .format(chrom, self.fa_name), file=sys.stderr, flush=True)
            return sequence

        if not pstart:
            pstart = 1

        if not pend:
            pend = self.fai_data[chrom]['chrom_len']

        region = '{}:{}-{}'.format(chrom, pstart, pend)

        if pstart < 1 or pstart > self.fai_data[chrom]['chrom_len']:
            print('*Warning* illegal query starting position {} (must be '
                  '1-based within range [1,{}] for chrom {})'
                  .format(pstart, self.fai_data[chrom]['chrom_len'], chrom),
                  file=sys.stderr, flush=True)
            return region, sequence

        if pend < 1 or pend > self.fai_data[chrom]['chrom_len']:
            print('*Warning* illegal query ending position {} (must be '
                  '1-based within range [1,{}] for chrom {})'
                  .format(pend, self.fai_data[chrom]['chrom_len'], chrom),
                  file=sys.stderr, flush=True)
            return region, sequence

        if pend < pstart:
            print('*Warning* query region is of negative length ({}:{}-{})'
                  .format(chrom, pstart, pend),
                  file=sys.stderr, flush=True)
            return region, sequence

        chrom_offset = self.fai_data[chrom]['byte_offset']
        line_ndiff = (self.fai_data[chrom]['line_nchar'] -
                      self.fai_data[chrom]['line_nbase'])
        line_nbase = self.fai_data[chrom]['line_nbase']
        offset_start = (chrom_offset + (pstart - 1) + line_ndiff *
                        ((pstart - 1) // line_nbase))
        offset_end = (chrom_offset + (pend - 1) + line_ndiff *
                        ((pend - 1) // line_nbase))
        read_nbyte = offset_end - offset_start + 1

        with open(self.fa_name) as f:
            f.seek(offset_start)
            sequence = f.read(read_nbyte).replace('\n', '')

        return region, sequence

# =============================================================================

def main():
    """Wrapper of function fasta.query_region() with command line inputs."""

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--fasta', metavar="fasta",
                        help="input FASTA file",
                        type=argparse.FileType('r'), required=True)
    parser.add_argument('-c', '--chrom', metavar="chrom_name",
                        help="chromosome name",
                        type=str, required=True)
    parser.add_argument('-a', '--pstart', metavar="pos_start",
                        help=("starting genomic position 1-based inclusive "
                              "(default: 1)"),
                        type=int, required=False, default=1)
    parser.add_argument('-b', '--pend', metavar="pos_end",
                        help=("ending genomic position 1-based inclusive "
                              "(default: chromosome length)"),
                        type=int, required=False, default=None)

    args = parser.parse_args()

    fasta_name = args.fasta.name
    chrom = args.chrom
    pstart = args.pstart
    pend = args.pend

    fa = Fasta(fasta_name)
    reg, seq = fa.query_region(chrom, pstart, pend)
    print('> {}'.format(reg), file=sys.stdout, flush=True)
    print(seq, file=sys.stdout, flush=True)


# =============================================================================

if __name__ == '__main__':

    main()
