# -*- coding: utf-8 -*-
"""This module contains functions dealing with fastq file.

    fastq_get_rgstr() derive RG string from input fastq file.

        Example:

            $ python3 fastq_rgstr.py -f fastq_name -s sample_name -n num_line

        This will sample "num_line" lines from fastq "fastq_name" and return an
        RG string with the "sample_name" in SM field, e.g.

            @RG ID:HWI-7001457:234:1 PL:ILLUMINA PU:HWI-7001457:234:C6UK1ANXX:1
            SM:QM4

"""

import sys
from subprocess import PIPE, Popen
import json
import os
import argparse


# =============================================================================

def fastq_get_rgstr(fastq_name, sample_name=None, head_nline=1,
                    platform='ILLUMINA'):
    """Generate RG string from input fastq file.

    This function derives ID and PU information for the RG string based on the
    reads in the input fastq file. The SM field in the RG string is either
    given in input, or derived from the fastq file name (the first field
    separated by dot).

    Args:
        fastq_name (string): fastq file name
        sample_name (string): sample name (default to 1st field in fastq file
            name is not given)
        head_nline (int): first number of lines to sample in fastq file
        platform (string): platform informtion for PL field

    Returns:
        RG string in format: '@RG\tID:{}\tPL:{}\tPU:{}\tSM:{}'


    Following the convention of RG string, suppose the read ID in the input
    fastq is of the form:

        @HWI-7001446:480:C6BH4ANXX:3:1101:2304:1985 1:N:0:CCTCCT

    The first column is the sequence identifier of a typical reads from
    Illumina sequencing machine (followed by the second column is the optional
    description).  The sequence identifier typically consists of:

        instrument id: HWI-7001446
        flowcell id (run_id): 480
        flowcell barcode: C6BH4ANXX
        lane id: 3
        tile id: 1101
        cluster x coordinate: 2304
        cluster y coordinate: 1985

    Accordingly, this module generates the RG ID and PU (platform unit) string
    with

        ID = instrument id : flowcell id : ... : lane id
        PU = instrument id : flowcell id : ... : flowcell barcode : lane id

    References:
    SAM specification v1.0
        http://samtools.github.io/hts-specs/SAMv1.pdf
    GATK definition for read groups
        https://software.broadinstitute.org/gatk/guide/article?id=6472
    GATK best practice for multiplexed sequencing
        https://software.broadinstitute.org/gatk/documentation/article?id=3060

    """

    # derive sample name from fastq name
    if not sample_name:
        sample_name = os.path.basename(fastq_name).split('.')[0]

    if fastq_name.endswith('gz'):
        cmd = 'zcat {}'.format(fastq_name)
    else:
        cmd = 'cat {}'.format(fastq_name)

    proc = Popen(cmd, stdout=PIPE, shell=True, universal_newlines=True)
    line_num = 0
    pu_str_counter = {}
    for line in proc.stdout:
        line_num += 1
        if (line_num-1) % 4 != 0:
            continue
        if head_nline and line_num > head_nline:
            break

        # exclude last 3 field in read id (tile id, cluster x coordinate,
        # cluster y coordinate)
        pu_str = ':'.join(line.rstrip().split()[0].split(':')[:-3])
        if pu_str not in pu_str_counter:
            pu_str_counter[pu_str] = 0
        pu_str_counter[pu_str] += 1
    proc.kill()

    if len(pu_str_counter) > 1:
        print('*WARNING*: found multiple flowcell_id/lane in fastq {}: {}'
              .format(fastq_name, json.dumps(pu_str_counter, sort_keys=True)),
              file=sys.stderr, flush=True)

    for pu_str in sorted(pu_str_counter):
        break

    pu_str = pu_str.lstrip('@').split(':')
    id_str = pu_str[:-2] + [pu_str[-1]]
    # instrument_id : flowcell_id : ... : lane_id
    id_str = ':'.join(id_str)
    # instrument_id : flowcell_id : ... : flowcell_barcode : lane_id
    pu_str = ':'.join(pu_str)
    rg_str = ('@RG\tID:{}\tPL:{}\tPU:{}\tSM:{}'
              .format(id_str, platform, pu_str, sample_name))

    return rg_str

# =============================================================================

def main():
    """Wrapper of function fastq_get_rgstr() with command line inputs."""

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--fastq', metavar="fastq",
                        help="input FASTQ",
                        type=argparse.FileType('r'), required=True)
    parser.add_argument('-s', '--sname', metavar="sample_name",
                        help="sample name (default: first field in file name)",
                        type=str, required=False, default=None)
    parser.add_argument('-n', '--nline', metavar="num_line",
                        help="number of fastq lines to sample (default: 1000)",
                        type=int, required=False, default=1000)

    args = parser.parse_args()

    fastq_name = args.fastq.name
    sample_name = args.sname
    head_nline = args.nline

    read_group = fastq_get_rgstr(fastq_name, sample_name, head_nline)
    print(read_group, file=sys.stdout, flush=True)

# =============================================================================

if __name__ == '__main__':

    main()
