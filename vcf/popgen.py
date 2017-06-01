'''
popgen v0.1

A suite of functions that calculate useful population genetics statistics
between a pair of VCF records - ideally SNPs with single ALTs. Deviations from this
won't break the script, but will throw up (overrideable) warning messages and 
should probably be done with caution.

If doing some exploratory work, it may be worthwhile to load in subsets of records
using reclist() and running functions on records out of those lists. 

If working with known sites, reclook() allows for easy access of records based on 
their genomic positions and can be entered as input to the pop gen functions. 

This is an extremely early draft of what might be put together into a full
Python library down the line.

AH - 05/2017
'''

import vcf

def snpchecker(record1, record2):
    '''Given two VCF records, checks whether they are SNPs with single ALTs each.
    Helper function that just throws warning messages if necessary.'''
    if record1.is_snp == False:
        print('caution: first record is not a SNP')
    elif record2.is_snp == False:
        print('caution: second record is not a SNP')
    if len(record1.ALT) > 1:
        print('caution: multiple alts in record 1 - ', record1.ALT)
    if len(record2.ALT) > 1:
        print('caution: multiple alts in record 2 - ', record2.ALT)

def freqscalc(record1, record2, snpcheck = True):
    '''Given two VCF records, returns observed haplotype frequencies.
    Will check that records are single-ALT SNPs unless snpcheck = False.'''
    if snpcheck == True:
        snpchecker(record1, record2)
    elif snpcheck == False:
        pass
    p = 1 - record1.aaf[0]
    q = 1 - record2.aaf[0]
    p2 = record1.aaf[0]
    q2 = record2.aaf[0]
    print(record1.CHROM, record1.POS, '- ref', record1.REF, 'alt', record1.ALT[0])
    print(record2.CHROM, record2.POS, '- ref', record2.REF, 'alt', record2.ALT[0])
    print('p1 ', p, 'p2 ', p2)
    print('q1 ', q, 'q2 ', q2)
    # get samples + check for same samples b/w both records
    strainlist = [record1.samples[i].sample for i in range(len(record1.samples))] 
    assert strainlist == [record2.samples[i].sample for i in range(len(record2.samples))]    
    # score haplotypes
    haplist = []
    for strain in strainlist:
        gt1 = record1.genotype(strain)['GT']
        gt2 = record2.genotype(strain)['GT']
        if gt1 == '.' or gt2 == '.':
            continue
        if gt1 == '1' and gt2 == '1':
            outgt = str(record1.ALT[0]) + str(record2.ALT[0]) 
        elif gt1 == '1' and gt2 == '0':
            outgt = str(record1.ALT[0]) + record2.REF
        elif gt1 == '0' and gt2 == '1':
            outgt = record1.REF + str(record2.ALT[0])
        elif gt1 == '0' and gt2 == '0':
            outgt = record1.REF + record2.REF
        haplist.append(outgt) # create list of observed genotypes
    uniques = set(haplist)
    for hap in uniques:
        print(hap, round(haplist.count(hap)/len(haplist), 5))   
        
        
def dcalc(record1, record2, snpcheck = True):
    '''Calculates D statistic between two VCF records.
    Will check that records are single-ALT SNPs unless snpcheck = False.'''
    if snpcheck == True:
        snpchecker(record1, record2)
    elif snpcheck == False:
        pass
    p = 1 - record1.aaf[0]
    q = 1 - record2.aaf[0]
    p2 = record1.aaf[0]
    q2 = record2.aaf[0]
    refcount = 0
    altcount = 0
    strainlist = [record1.samples[i].sample for i in range(len(record1.samples))] 
    assert strainlist == [record2.samples[i].sample for i in range(len(record2.samples))]    
    for strain in strainlist:
        gt1 = record1.genotype(strain)['GT']
        gt2 = record2.genotype(strain)['GT']
        if gt1 == '.' or gt2 == '.':
            continue
        elif gt1 == '0' and gt2 == '0':
            refcount = refcount + 1
        else:
            altcount = altcount + 1
    ABfreq = refcount/(refcount + altcount)
    d = ABfreq - (p * q)
    d = round(d, 5)
    return d      

def dprimecalc(record1, record2, snpcheck = True):
    """Calculates Lewontin's D' statistic (D/Dmax) between two VCF records.
    Will check that records are single-ALT SNPs unless snpcheck = False."""
    if snpcheck == True:
        snpchecker(record1, record2)
    elif snpcheck == False:
        pass
    p = 1 - record1.aaf[0]
    q = 1 - record2.aaf[0]
    p2 = record1.aaf[0]
    q2 = record2.aaf[0]
    d = dcalc(record1, record2)
    if d >= 0:
        dmax = min(p * q2, p2 * q)
        out = d/dmax
    elif d < 0:
        dmin = max(-1 * p * q, -1 * p2 * q2)
        out = d/dmin
    elif d == 0:
        print("Error - SNPs in linkage equilibrium")
    out = round(out, 4)
    return out

def r2calc(record1, record2, snpcheck = True):
    '''Calculates r^2 (correlation) between two VCF records.
    Will check that records are single-ALT SNPs unless snpcheck = False.'''
    if snpcheck == True:
        snpchecker(record1, record2)
    elif snpcheck == False:
        pass
    p = 1 - record1.aaf[0]
    q = 1 - record2.aaf[0]
    p2 = record1.aaf[0]
    q2 = record2.aaf[0]
    dsquared = dcalc(record1, record2)**2
    out = dsquared/(p * q * p2 * q2)
    out = round(out, 4)
    return out

def ldstats(record1, record2, snpcheck = True, freqs = False):
    '''Convenience function - returns D, D prime, and r2.
    If freqs = True, will also return an allele frequency report.
    Will check that records are single-ALT SNPs unless snpcheck = False.'''
    if snpcheck == True:
        snpchecker(record1, record2)
    elif snpcheck == False:
        pass
    print('D =', dcalc(record1, record2, snpcheck = False))
    print("D'=", dprimecalc(record1, record2, snpcheck = False)) 
    print("r2 =", r2calc(record1, record2, snpcheck = False))
    if freqs == True:
        print('\nFrequencies report:')
        freqscalc(record1, record2, snpcheck = False) 
    
def reclist(vcf_file, chrom = None, pos = None, snpsonly = False):
    '''Returns records in given gzipped VCF file as a list.
    If given chrom, will fetch just chrom; if given both chrom 
    and pos (in the format 'start-end') will fetch just that 
    subset of the VCF.'''
    vcfin = vcf.Reader(filename = vcf_file, compressed = True)
    if chrom is not None and pos is not None:
        try:
            assert isinstance(pos, str) 
            pos = pos.split('-')
            assert len(pos) == 2
            start = int(pos[0]) - 1
            end = int(pos[1])
            snippet = vcfin.fetch(chrom = chrom, start = start, end = end) 
            if snpsonly == True:
                reclist = [record for record in snippet if record.is_snp == True]
            elif snpsonly == False:
                reclist = [record for record in snippet]
        except:
            print('Error in pos parameter.')
            print('Please enter positions in a start-end format with no spaces.')
    elif chrom is not None and pos is None:
        snippet = vcfin.fetch(chrom = chrom)
        if snpsonly == True:
            reclist = [record for record in snippet if record.is_snp == True]
        elif snpsonly == False:
            reclist = [record for record in snippet]
    elif chrom is None and pos is not None:
        print('Error - pos supplied without chrom specification.')
    else:
        if snpsonly == True:
            reclist = [record for record in vcfin if record.is_snp == True]
        elif snpsonly == False:
            reclist = [record for record in vcfin]
    return reclist

def reclook(pos, reclist):
    '''Convenience function. If VCF records are saved to a list 
    (ie via reclist()), will return the record at an input 
    position value. Can be used as input to pop gen functions.'''
    for record in reclist:
        if record.POS == pos:
            return record
        else:
            pass


