
# TODO(Zhuoyi): add docstring


class Error(Exception):
    pass


class RegionRangeError(Error):
    pass


class InternalError(Error):
    # coding level mistakes
    pass


class Region(object):
    """A class for genomic region defintion.

    Methods:
        __init__(): initialize a Region object
        __str__(): return a string version of region "chrom:pstart-pend"

    Attributes:
        chrom: chromosome name of the region
        pstart: region starting position (int) 1-based inclusive
        pend: region ending position (int) 1-based inclusive
        length: number of base-pairs in the region (int) -1 if region is empty
    """

    def __init__(self, *args):
        """Initialize a Region object.

        Initialize a Region object given optional inputs of chrom, starting and
        ending position of the region.

        Args:
            *args: number of elements is between 0 and 3
            args[0]: chromosome name (string)
            args[1]: region starting position 1-based inclusive (int)
            args[1]: region ending position 1-based inclusive (int)
            if no argument is given, init an emtpy region object with length -1
            if 1 arg is given , init a region for whole chromosome (len=-1)
            if 2 args are given , init a region of 1-base long

        Returns:

        Raises:
            InternalError: error in calling code (mostly with input args)
            RegionRangeError: invalid input region, pstart<=0 or pend<=0
                or pend < pstart
            ValueError: input argument type (string/int) is wrong
        """
        self.chrom = ''
        self.pstart = self.pend = self.length = -1
        if len(args) == 0:
            return

        if len(args) > 3:
            raise InternalError('initialize Region object with max 3 args\n')

        if ':' in args[0] or '-' in args[0]:
            raise InternalError('1st argument must be chrom/contig name\n')

        if len(args) >= 1:
            chrom = args[0]

        if len(args) == 2:
            pstart = int(args[1])
            if pstart <= 0:
                raise RegionRangeError('region starting position must be >=1'
                                       ' (found {})'.format(args[1]))
            pend = pstart
            length = 1

        elif len(args) == 3:
            pstart = int(args[1])
            pend = int(args[2])
            if pstart <= 0 or pend <= 0:
                raise RegionRangeError('region starting/ending position must '
                                       'be >=1 (found {} {})'.format(args[1],
                                                                     args[2]))
            if pend < pstart:
                raise RegionRangeError('region is empty or negative length '
                                       '(found {} {})'.format(args[1],args[2]))
            length = pend - pstart + 1

        self.chrom, self.pstart, self.pend, self.length \
            = chrom, pstart, pend, length


    def __str__(self):
        """Format a region object for print (chrom:pstart-pend).

        Returns:
            a string of format "chrom(:pstart-pend)". pstart, pend are 1-based
            inclusive (if available).
        """
        if self.length == -1:
            return self.chrom
        return '{}:{}-{}'.format(self.chrom, self.pstart, self.pend)
