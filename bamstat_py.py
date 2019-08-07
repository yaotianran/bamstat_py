#!/usr/bin/env python3
import os, sys, collections, statistics, pprint, math
import pysam
import os.path as path

SELFPATH = path.dirname(__file__)
MINIMUM_COVERAGE = 5
MINIMUM_MAP_REGION_GAP = 150

def pileup_stats(bamfile: str) -> dict:
    '''
    This function tries to guess the enriched regions depending on the coverage.
    Any loci that has more coverages than MINIMUM_COVERAGE would be considered
    as a enriched region "break point". The regions must not be closer than MINIMUM_MAP_REGION_GAP.

    Parameters:
        **bamfile**: string
            The bam file to be calculated. We assume it has been coordinate-sorted and indexed

    Returns:
        **None**: None
            None on success

    '''
    global MINIMUM_COVERAGE, MINIMUM_MAP_REGION_GAP

    loci = collections.namedtuple('loci', ['chrom', 'pos'])   # loci(str chrom, int pos)
    amplicon = collections.namedtuple('amplicon', ['start', 'end'])   # amplicon(loci start, loci end)

    coverage_dict = collections.OrderedDict() #  {loci1: converage_int1, loci2: converage_int2.......}
    map_qualtity_dict = collections.OrderedDict()  # {loci1: mean_quality1, loci2: mean_quality2.......}
    bamfile_af = pysam.AlignmentFile(bamfile, 'rb')
    s = ''
    for pileupcolumn in bamfile_af.pileup(min_base_quality= 0, ignore_overlaps= False, ignore_orphans= False):
        if pileupcolumn.reference_name != s:
            print('Searching enriched regions at contig {}...'.format(pileupcolumn.reference_name))
            s = pileupcolumn.reference_name

        current_loci = loci(chrom= pileupcolumn.reference_name, pos= pileupcolumn.reference_pos)
        try:
            current_coverage_int = pileupcolumn.nsegments
        except:
            message = 'Fail to get coverage number at "{chr}:{pos}". Skip'.format(chr= pileupcolumn.reference_name, pos= pileupcolumn.reference_pos)
            print(message)
            next

        if current_coverage_int >= MINIMUM_COVERAGE:
            coverage_dict[current_loci] = current_coverage_int
            try:
                map_qualtity_dict[current_loci] = round(statistics.mean(pileupcolumn.get_mapping_qualities()))
            except Exception as ex:
                print('An error occured when deal mapping quality at "{chr}:{pos}". Skip.'.format(chr= pileupcolumn.reference_name, pos= pileupcolumn.reference_pos))
                print(ex)
                next
    bamfile_af.close()

    if len(coverage_dict) == 0:
        raise ValueError('Fail to get any coverage info from bam file.')
    elif len(coverage_dict) == 1:
        print(coverage_dict.keys(), 'coverage: ', coverage_dict.values())
        raise ValueError('Bam file not enough coverage info.')
    else:
        loci_lst = list(coverage_dict.keys())

    # assign each ''
    start_loci_lst = [loci_lst[0]]
    end_loci_lst = []
    for i in range(1, len(loci_lst)-1):
        current_loci, next_loci = loci_lst[i], loci_lst[i+1]
        if next_loci.chrom != current_loci.chrom or next_loci.pos - current_loci.pos >= MINIMUM_MAP_REGION_GAP:   # a break point exists between current loci and next loci
            start_loci_lst.append(next_loci)
            end_loci_lst.append(current_loci)

    end_loci_lst.append(loci_lst[-1])

    # zip the two amplicons lists
    if len(start_loci_lst) == len(end_loci_lst):
        amplicons_lst = []
        for i in range(len(start_loci_lst)):
            amplicons_lst.append( amplicon(start= start_loci_lst[i], end= end_loci_lst[i]) )
    else:
        with open('error.txt', 'wt') as error_f:
            error_f.writelines('Start locis:\n')
            for loc in start_loci_lst:
                message = '\t{chr}:{pos}\n'.format(chr=loc.chrom, pos=loc.pos)
                error_f.writelines(message)
            error_f.writelines('End locis:\n')
            for loc in end_loci_lst:
                message = '\t{chr}:{pos}\n'.format(chr=loc.chrom, pos=loc.pos)
                error_f.writelines(message)
        message = 'The start locis (N={}) do not match end loci (N={}), check error.txt for detail.'.format(len(start_loci_lst), len(end_loci_lst))
        raise ValueError(message)

    # calculate the map quality and coverage quantile, and print out the results
    print('\nChr\tStart\tEnd\tLength\tMap_Qual_Avg.\tCov_Avg.\tCov_Min\tCov_25\tCov_Med.\tCov_75\tCov_max')
    for amp in amplicons_lst:
        map_quality_lst = []
        coverage_lst = []
        for pos_int in range(amp.start.pos, amp.end.pos+1):
            try:
                map_quality_lst.append(map_qualtity_dict[loci(chrom=amp.start.chrom, pos= pos_int)])
            except KeyError:
                message = 'Fail to get map quality at {chr}:{pos}. Skip.'.format(chr= amp.start.chrom, pos= pos_int)
                #print(message)
                next

            try:
                coverage_lst.append( coverage_dict[loci(chrom=amp.start.chrom, pos= pos_int)] )
            except KeyError:
                message = 'Fail to get coverage info at {chr}:{pos}. Skip.'.format(chr= amp.start.chrom, pos= pos_int)
                #print(message)
                next

        if len(map_quality_lst) == 0 or len(coverage_lst) == 0:
            message = 'Amplicon ({chr}:{start}-{end}) contains neither map quality nor coverage info. Skip.'.format(chr= amp.start.chrom,
                                                                                                                    start= amp.start.pos,
                                                                                                                    end= amp.end.pos)
            print(message)
            next
        else:
            map_quality_lst.sort()
            coverage_lst.sort()

        try:
            message = '{chr}\t{start:>9}\t{end:>9}\t{length}\t{map_avg}\t{cov_avg}\t{cov_min}\t{cov_25}\t{cov_50}\t{cov_75}\t{cov_max}'.format(chr= amp.start.chrom,
                                                                                                                                     start= amp.start.pos,
                                                                                                                                     end= amp.end.pos,
                                                                                                                                     length= amp.end.pos-amp.start.pos+1,
                                                                                                                                     map_avg= round(statistics.mean(map_quality_lst)),
                                                                                                                                     cov_avg= round(statistics.mean(coverage_lst)),
                                                                                                                                     cov_min= min(coverage_lst),
                                                                                                                                     cov_25= coverage_lst[math.floor(len(coverage_lst)*0.25)],
                                                                                                                                     cov_50= statistics.median(coverage_lst),
                                                                                                                                     cov_75= coverage_lst[math.floor(len(coverage_lst)*0.75)],
                                                                                                                                     cov_max= max(coverage_lst) )
            print(message)
        except statistics.StatisticsError:
            message = 'An error occured when statisticing amplicon ({chr}:{start}-{end}). Skip.'.format(chr= amp.start.chrom,
                                                                                                        start= amp.start.pos,
                                                                                                        end= amp.end.pos)
            print(message)
            next

    return None


def bed_stats(bamfile: str, bedfile: str) -> dict:
    '''
    This function simply print out some coverage statistics for the regions in a bed file.

    Parameters:
        **bamfile**: string
            The bam file to be calculated. We assume it has been coordinate-sorted and indexed

        **bedfile**: string
            The regions in this bed file will be checked. Only first three columns will be used.

    Returns:
        **None**:
            Return None on success
    '''

    loci = collections.namedtuple('loci', ['chrom', 'pos'])  # loci(str chrom, int pos)
    amplicon = collections.namedtuple('amplicon', ['start', 'end'])   # amplicon(loci start, loci end)
    amplicons_lst = []
    with open(bedfile, 'rt') as bedfile_f:
        for line_str in bedfile_f.readlines():
            try:
                if line_str.strip() != '':
                    line_lst = line_str.strip().split('\t')
            except:
                next

            try:
                amplicons_lst.append( amplicon(start=loci(chrom=line_lst[0], pos=int(line_lst[1])), end=loci(chrom=line_lst[0], pos=int(line_lst[2]))) )
            except IndexError:
                message = 'BED file format error. BED file must have at least three columns. Skip\n{}'.format(line_str)
                print(message)
                next
            except ValueError:
                message = 'BED file format error. The second and third column must be integers. Skip\n{}'.format(line_str)
                print(message)
                next
            except:
                message = 'An error occur when parse BED file. Skip\n{}'.format(line_str)
                print(message)
                next

    if len(amplicons_lst) == 0:
        message = 'BED file contains no amplicon region. Check the BED file.'
        raise ValueError(message)

    # calculate map quality and coverage
    print('Chr\tStart\tEnd\tLength\tMap_Qual_Avg.\tCov_Avg.\tCov_Min\tCov_25\tCov_Med.\tCov_75\tCov_max')
    bamfile_af = pysam.AlignmentFile(bamfile, 'rb')
    for amp in amplicons_lst:
        #print('Checking amplicon ({chr}:{start}-{end})...'.format(chr= amp.start.chrom, start= amp.start.pos, end= amp.end.pos), end='\r', flush=True)
        coverage_lst = []
        map_quality_lst = []
        for pileupcolumn in bamfile_af.pileup(contig= amp.start.chrom, start= amp.start.pos, stop= amp.end.pos, min_base_quality= 0, ignore_overlaps= False, ignore_orphans= False):
            try:
                coverage_lst.append(pileupcolumn.nsegments)
            except:
                pass

            try:
                map_quality_lst.append(round(statistics.mean(pileupcolumn.get_mapping_qualities())))
            except:
                pass

        if len(map_quality_lst) == 0 or len(coverage_lst) == 0:
            message = 'Amplicon ({chr}:{start}-{end}) contains neither map quality nor coverage info. Skip.'.format(chr= amp.start.chrom,
                                                                                                                    start= amp.start.pos,
                                                                                                                    end= amp.end.pos)
            print(message)
            next
        else:
            map_quality_lst.sort()
            coverage_lst.sort()

        try:
            message = '{chr}\t{start:>9}\t{end:>9}\t{length}\t{map_avg}\t{cov_avg}\t{cov_min}\t{cov_25}\t{cov_50}\t{cov_75}\t{cov_max}'.format(chr= amp.start.chrom,
                                                                                                                                         start= amp.start.pos,
                                                                                                                                         end= amp.end.pos,
                                                                                                                                         length= amp.end.pos-amp.start.pos+1,
                                                                                                                                         map_avg= round(statistics.mean(map_quality_lst)),
                                                                                                                                         cov_avg= round(statistics.mean(coverage_lst)),
                                                                                                                                         cov_min= min(coverage_lst),
                                                                                                                                         cov_25= coverage_lst[round(len(coverage_lst)*0.25)],
                                                                                                                                         cov_50= statistics.median(coverage_lst),
                                                                                                                                         cov_75= coverage_lst[round(len(coverage_lst)*0.75)],
                                                                                                                                         cov_max= max(coverage_lst) )
            print(message)
        except statistics.StatisticsError:
            message = 'An error occured when statisticing amplicon ({chr}:{start}-{end}). Skip.'.format(chr= amp.start.chrom,
                                                                                                        start= amp.start.pos,
                                                                                                        end= amp.end.pos)
            print(message)
            next

    bamfile_af.close()
    return None

def segment_stats(bamfile: str) -> dict:






    return




def main(argvList = sys.argv, argv_int = len(sys.argv)):
    '''
    This is the function explanation.

    Parameters:
        **parameter 1**: type
            its function

        **parameter 2**: type
            its function

    Returns:
        **value**: type
            its meaning
    '''

    if argv_int == 2:
        pileup_stats(bamfile= path.realpath(path.expanduser(argvList[1])) )
    elif argv_int == 3:
        bed_stats(bamfile= path.realpath(path.expanduser(argvList[1])), bedfile= path.realpath(path.expanduser(argvList[2])))
    else:
        print('Wrong parameters:\n{}'.format(argvList))
        sys.exit(0)

    return

main()
'''
bamfile = '/home/user/H04135D-T1_sort.bam'
#bamfile = '/home/user/S1_region.bam'
bedfile = '/home/user/test.bed'

#pileup_stats(bamfile)
bed_stats(bamfile, bedfile)

print()


bamfile_af = pysam.AlignmentFile(bamfile, 'rb')
#for segment in bamfile_af.fetch():
#    print(segment.query_name)

for pileupcolumn in bamfile_af.pileup(min_base_quality=0 , ignore_overlaps= False, ignore_orphans= False):
    print()
    print(pileupcolumn.reference_name, pileupcolumn.reference_pos, pileupcolumn.get_query_positions(), pileupcolumn.get_mapping_qualities())
    print(pileupcolumn.get_query_sequences(add_indels = True))
    i = 0
    seq_lst = []
    map_lst = []
    for pileupread in pileupcolumn.pileups:
        if pileupread.is_del:
            print('is_del')
        elif pileupread.is_refskip:
            print(pileupread.aligment.query_name, 'is_refskip')
        else:
            i += 1
            seq_lst.append(pileupread.alignment.query_sequence[pileupread.query_position])
            #map_lst.append(pileupread.alignment.query_sequence[pileupread.query_position])
    print(i, seq_lst)

bamfile_af.close()
'''