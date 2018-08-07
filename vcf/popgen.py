'''
popgen v0.6
(requires pyVCF)

Functions that calculate useful population genetics statistics
between a pair of VCF records - ideally SNPs with single ALTs. Deviations from this
won't break anything, but will throw up (overrideable) warning messages and 
should probably be done with caution.

If doing some exploratory work, it may be worthwhile to load in subsets of records
using reclist() and run functions on records out of the resultant list objects. 

If working with known sites, reclook() allows for easy access of records based on 
their genomic positions and can be entered as input to the pop gen functions. 

These functions are only compatible with haploid genomes (for now). 
Use with higher ploidies at your own risk.

This is an extremely early draft of what might be put together into a full
Python library down the line. 

AH - 06/2017
'''

import vcf

def snpchecker(record1, record2):
    '''Given two VCF records, checks whether they are SNPs with single ALTs each.
    Helper function that just throws warning messages if necessary.
    '''
    if len(record1.REF) == 1 and len(record1.ALT) == 1 and len(record1.ALT[0]) == 1 and len(record1.alleles) == 2:
        pass
    elif len(record1.REF) != 1 or len(record1.ALT) != 1 or len(record1.ALT[0]) != 1 or len(record1.alleles) != 2:
        print('caution: first record is not a biallelic SNP')
        print(record1.REF, record1.ALT, record1.alleles)
    if len(record2.REF) == 1 and len(record2.ALT) == 1 and len(record2.ALT[0]) == 1 and len(record2.alleles) == 2:
        pass
    elif len(record2.REF) != 1 or len(record2.ALT) != 1 or len(record2.ALT[0]) != 1 or len(record2.alleles) == 2:
        print('caution: second record is not a biallelic SNP')
        print(record2.REF, record2.ALT, record2.alleles)

def straingetter(record1, record2, GQ_threshold = 30):
    '''Given two VCF records, returns a list of individuals in the population that
    contain calls at both sites. Helper function for LD calculations.
    Will filter strains for GQ (default 30).
    '''
    rec1set = set([record1.samples[i].sample for i in range(len(record1.samples)) if record1.samples[i]['GT'] != '.'])
    rec2set = set([record2.samples[i].sample for i in range(len(record2.samples)) if record2.samples[i]['GT'] != '.'])
    strainlist = list(rec1set.intersection(rec2set))
    # filter for GQ
    garbage = []
    for strain in strainlist:
        if record1.genotype(strain)['GQ'] < GQ_threshold or record2.genotype(strain)['GQ'] < GQ_threshold:
            garbage.append(strain)
    if len(garbage) > 0:
        for trash_strain in garbage:
            strainlist.remove(trash_strain)
    return strainlist

def freqsgetter(record1, record2, snpcheck = True):
    '''Helper function for LD statistic calculations. Returns, in order: 1. a dict containing
    haplotype frequencies; 2. a dict containing allele frequency values (p1, p2, q1, q2), and
    3. a dict containing haplotype frequencies with A/B notation. Has more accurate allele
    freq calculations than freqscalc() - handles missing data better (i.e. when biallelic genotype
    not present).
    '''
    if snpcheck:
        snpchecker(record1, record2)
    elif not snpcheck:
        pass
    # check strains b/w compared records are identical
    strainlist = straingetter(record1, record2)
    # parse through VCF calls
    haplist = []
    pcount = 0
    qcount = 0
    totcalls = 0
    for strain in strainlist:
        outgt = ''
        gt1 = record1.genotype(strain)['GT']
        gt2 = record2.genotype(strain)['GT']
        if gt1 == '.' or gt2 == '.':
            return
        elif gt1 == '1' and gt2 == '1': # ab
            outgt = str(record1.ALT[0]) + str(record2.ALT[0])
            totcalls = totcalls + 1
        elif gt1 == '1' and gt2 == '0': # aB
            qcount = qcount + 1
            outgt = str(record1.ALT[0]) + record2.REF
            totcalls = totcalls + 1
        elif gt1 == '0' and gt2 == '1': # Ab
            pcount = pcount + 1
            outgt = record1.REF + str(record2.ALT[0])
            totcalls = totcalls + 1
        elif gt1 == '0' and gt2 == '0': # AB
            pcount = pcount + 1
            qcount = qcount + 1
            outgt = record1.REF + record2.REF
            totcalls = totcalls + 1
        haplist.append(outgt)
    # assign allele freq values
    values = {}
    if totcalls == 0:
        values['p1'] = 0
        values['q1'] = 0
    else:
        values['p1'] = pcount/totcalls
        values['q1'] = qcount/totcalls
    values['p2'] = 1 - values['p1']
    values['q2'] = 1 - values['q1']
    # haplotype frequencies
    uniques = dict.fromkeys(set(haplist))
    for hap in uniques:
        uniques[hap] = haplist.count(hap)/len(haplist)
    # create tuples w/ actual haps and corresponding AB notation
    homref = record1.REF + record2.REF, 'AB'
    homalt = str(record1.ALT[0]) + str(record2.ALT[0]), 'ab'
    het1 = record1.REF + str(record2.ALT[0]), 'Ab'
    het2 = str(record1.ALT[0]) + record2.REF, 'aB'
    haps = {}
    # use tuples to assign hap freqs to AB notation for downstream use
    genotypes = [homref, homalt, het1, het2]
    for genotype in genotypes:
        if genotype[0] in uniques.keys():
            haps[genotype[1]] = uniques[genotype[0]]
        else:
            haps[genotype[1]] = 0
    return uniques, values, haps

def dcalc(record1, record2, snpcheck = True):
    '''Calculates D statistic between two VCF records.
    Will check that records are biallelic (single-ALT) SNPs unless snpcheck = False.
    '''
    if snpcheck:
        snpchecker(record1, record2)
    elif not snpcheck:
        pass    
    haps = freqsgetter(record1, record2, snpcheck = False)[2]
    try:
        LHS = haps['AB'] * haps['ab']
    except KeyError: # either hap missing
        LHS = 0
    try:
        RHS = haps['Ab'] * haps['aB']
    except KeyError:
        RHS = 0
    d = LHS - RHS
    # d = round(d, 5)
    return d 

def dprimecalc(record1, record2, snpcheck = True):
    '''Calculates Lewontin's D' statistic (D/Dmax) between two VCF records.
    Will check that records are biallelic (single-ALT) SNPs unless snpcheck = False.
    '''
    if snpcheck:
        snpchecker(record1, record2)
    elif not snpcheck:
        pass
    values = freqsgetter(record1, record2, snpcheck = False)[1] # get allele frequencies
    d = dcalc(record1, record2, snpcheck = False)
    if d >= 0:
        dmax = min(values['p1'] * values['q2'], values['p2'] * values['q1'])
        if dmax == 0:
            out = 0
        else:
            out = d/dmax
    elif d < 0:
        dmin = max(-1 * values['p1'] * values['q1'], -1 * values['p2'] * values['q2'])
        if dmin == 0:
            out = 0
        else:
            out = d/dmin
    elif round(d, 6) == 0:
        out = 0
    return out

def r2calc(record1, record2, snpcheck = True):
    '''Calculates r^2 (correlation) between two VCF records.
    Will check that records are biallelic (single-ALT) SNPs unless snpcheck = False.
    '''
    if snpcheck:
        snpchecker(record1, record2)
    elif not snpcheck:
        pass
    values = freqsgetter(record1, record2, snpcheck = False)[1]
    if values['p1'] == 0 or values['q1'] == 0 or values['p2'] == 0 or values['q2'] == 0:
        out = 0
    else:
        dsquared = dcalc(record1, record2, snpcheck = False)**2
        out = dsquared/(values['p1'] * values['q1'] * values['p2'] * values['q2'])
        # out = round(out, 4)
    return out

def ldstats(record1, record2, snpcheck = True, freqs = False):
    '''Convenience function - returns D, D prime, and r2.
    If freqs = True, will also return an allele frequency report.
    Will check that records are biallelic (single-ALT) SNPs unless snpcheck = False.
    '''
    if snpcheck:
        snpchecker(record1, record2)
    elif not snpcheck:
        pass
    print('D =', dcalc(record1, record2, snpcheck = False))
    print("D'=", dprimecalc(record1, record2, snpcheck = False)) 
    print("r2 =", r2calc(record1, record2, snpcheck = False))
    if freqs == True:
        print('\nFrequencies report:')
        freqscalc(record1, record2, snpcheck = False) 

### exploratory functions        
        
def reclist(vcf_file, chrom = None, start = None, end = None, snpsonly = False):
    '''Returns records in given bgzipped VCF file as a list.
    If given chrom, will fetch just chrom; if given both chrom 
    and start-end coordinates, will fetch just that 
    subset of the VCF.
    '''
    
    vcfin = vcf.Reader(filename = vcf_file, compressed = True)
    outlist = []
    
    # filters
    def hardsnpcheck(record):
        if len(record.REF) == 1 and len(record.ALT) == 1 and len(record.ALT[0]) == 1 and len(record.alleles) == 2:
            return True
        elif len(record.REF) != 1 or len(record.ALT) != 1 or len(record.ALT[0]) != 1 or len(record.alleles) != 2:
            return False
    
    def issingleton(record): # ensure not singleton
        if type(record.INFO['AN']) == list:
            count = record.INFO['AN'][0] - record.INFO['AC'][0]
        else:
            count = record.INFO['AN'] - record.INFO['AC'][0]
        if count == 1:
            return True
        if record.INFO['AC'][0] == 1:
            return True
        else:
            return False
        
    def isinvariant(record):
        if record.INFO['AC'][0] == 0:
            return True
        if record.INFO['AF'][0] == 1.0:
            return True
        else:
            return False
    
    # make list
    if chrom and start and end:
        snippet = vcfin.fetch(chrom = chrom, start = start - 1, end = end) 
        if snpsonly:
            outlist = [record for record in snippet if hardsnpcheck(record) and not issingleton(record) and not isinvariant(record)]
        elif not snpsonly:
            outlist = [record for record in snippet]
    elif chrom and not start and not end:
        snippet = vcfin.fetch(chrom = chrom)
        if snpsonly:
            outlist = [record for record in snippet if hardsnpcheck(record) and not issingleton(record) and not isinvariant(record)]
        elif not snpsonly:
            outlist = [record for record in snippet]
    elif not chrom and start and end:
        print('Error - positions supplied without chromosome specification.')
    else:
        if snpsonly:
            outlist = [record for record in vcfin if hardsnpcheck(record) and not issingleton(record) and not isinvariant(record)]
        elif not snpsonly:
            outlist = [record for record in vcfin]
    return outlist

def reclook(reclist, pos):
    '''Convenience function. If VCF records are saved to a list 
    via reclist(), will return the record at an input 
    position value. Can be used as input to pop gen functions or ldstats().
    '''
    for record in reclist:
        if record.POS == pos:
            return record
        else:
            pass

def freqscalc(record1, record2, snpcheck = True, aaf = False):
    '''Exploratory convenience function. Given two VCF records, returns 
    observed haplotype frequencies. Will check that records are biallelic (single-ALT) SNPs 
    unless snpcheck = False. aaf will return AF values hardcoded in the VCF
    itself, while aaf = False (default) will make freqscalc calculate them instead
    (more accurate option).
    '''
    if snpcheck:
        snpchecker(record1, record2)
    elif not snpcheck:
        pass
    # get allele frequencies
    if aaf == True:
        p = 1 - record1.aaf[0]
        q = 1 - record2.aaf[0]
        p2 = record1.aaf[0]
        q2 = record2.aaf[0]
    elif aaf == False:
        values = freqsgetter(record1, record2)[1]
        p = values['p1']
        q = values['q1']
        p2 = values['p2']
        q2 = values['q2']
    print(record1.CHROM, record1.POS, '- ref', record1.REF, 'alt', record1.ALT[0])
    print(record2.CHROM, record2.POS, '- ref', record2.REF, 'alt', record2.ALT[0])
    print('p1 ', p, 'p2 ', p2)
    print('q1 ', q, 'q2 ', q2)
    # get samples + check for same samples b/w both records
    strainlist = straingetter(record1, record2)   
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
        
def singlegtcounts(record, showlist = False):
    '''Exploratory convenience function. For a given record, returns
    counts of refs, alts, and missing calls in the population. Can take
    a reclook function as input.
    '''
    samplelen = len(record.samples)
    gtlist = [record.samples[i]['GT'] for i in range(samplelen)]
    print('ref:', gtlist.count('0'))
    print('alt:', gtlist.count('1'))
    print('missing:', gtlist.count('.'))
    print('total:', samplelen)
    if showlist == True:
        print(gtlist)
    
def doublegtcounts(record1, record2, freqs = True, missing = True):
    '''Exploratory convenience function. For a given pair of records,
    prints all observed haplotypes in the population.
    '''
    # check strains b/w compared records are identical
    strainlist = [record1.samples[i].sample for i in range(len(record1.samples))] 
    assert strainlist == [record2.samples[i].sample for i in range(len(record2.samples))]
    if freqs:
        freqscalc(record1, record2)
    elif not freqs:
        pass
    # parse through VCF calls
    for strain in strainlist:
        gt1 = record1.genotype(strain)['GT']
        gt2 = record2.genotype(strain)['GT']
        if gt1 == '.' and gt2 == '.' and missing:
            print('--')
            continue
        elif gt1 == '.' and missing:
            if gt2 == '0':
                print('-B', '-' + record2.REF)
            elif gt2 == '1':
                print('-b', '-' + str(record2.ALT[0]))
            else:
                continue
        elif gt1 == '0':
            if gt2 == '.' and missing:
                print('A-', record1.REF + '-')
            elif gt2 == '0':
                print('AB', record1.REF + str(record2.REF))
            elif gt2 == '1':
                print('Ab', record1.REF + str(record2.ALT[0]))
        elif gt1 == '1':
            if gt2 == '.' and missing:
                print('a-', str(record1.ALT[0]) + '-')
            elif gt2 == '0':
                print('aB', str(record1.ALT[0]) + record2.REF)
            elif gt2 == '1':
                print('ab', str(record1.ALT[0]) + str(record2.ALT[0]))
                          
        
### experimental functions - LD at three loci (based on Genetic Data Analysis by Weir)

def triplefreqsgetter(record1, record2, record3):
    
    def triplestraingetter(record1, record2, record3):
        rec1set = set([record1.samples[i].sample for i in range(len(record1.samples)) if record1.samples[i]['GT'] != '.'])
        rec2set = set([record2.samples[i].sample for i in range(len(record2.samples)) if record2.samples[i]['GT'] != '.'])
        rec3set = set([record3.samples[i].sample for i in range(len(record3.samples)) if record3.samples[i]['GT'] != '.'])
        strainlist = list(rec1set.intersection(rec2set).intersection(rec3set))
        return strainlist
    
    # check strains b/w compared records are identical
    strainlist = triplestraingetter(record1, record2, record3)
    # parse through VCF calls
    haplist = []
    pcount = 0 # rec 1
    qcount = 0 # rec 2
    rcount = 0 # rec 3
    totcalls = 0
    for strain in strainlist:
        outgt = ''
        gt1 = record1.genotype(strain)['GT']
        gt2 = record2.genotype(strain)['GT']
        gt3 = record3.genotype(strain)['GT']
        if gt1 == '.' or gt2 == '.' or gt3 == '.':
            return
        elif gt1 == '1' and gt2 == '1' and gt3 == '1': # abc
            outgt = str(record1.ALT[0]) + str(record2.ALT[0]) + str(record3.ALT[0])
            totcalls = totcalls + 1
        elif gt1 == '1' and gt2 == '0' and gt3 == '1': # aBc
            qcount = qcount + 1
            outgt = str(record1.ALT[0]) + record2.REF + str(record3.ALT[0])
            totcalls = totcalls + 1
        elif gt1 == '0' and gt2 == '1' and gt3 == '1': # Abc
            pcount = pcount + 1
            outgt = record1.REF + str(record2.ALT[0]) + str(record3.ALT[0])
            totcalls = totcalls + 1
        elif gt1 == '0' and gt2 == '0' and gt3 == '1': # ABc
            pcount = pcount + 1
            qcount = qcount + 1
            outgt = record1.REF + record2.REF + str(record3.ALT[0])
            totcalls = totcalls + 1
        elif gt1 == '1' and gt2 == '1' and gt3 == '0': # abC
            rcount = rcount + 1
            outgt = str(record1.ALT[0]) + str(record2.ALT[0]) + record3.REF
            totcalls = totcalls + 1
        elif gt1 == '1' and gt2 == '0' and gt3 == '0': # aBC
            qcount = qcount + 1
            rcount = rcount + 1
            outgt = str(record1.ALT[0]) + record2.REF + record3.REF
            totcalls = totcalls + 1
        elif gt1 == '0' and gt2 == '1' and gt3 == '0': # AbC
            pcount = pcount + 1
            rcount = rcount + 1
            outgt = record1.REF + str(record2.ALT[0]) + record3.REF
            totcalls = totcalls + 1
        elif gt1 == '0' and gt2 == '0' and gt3 == '0': # ABC
            pcount = pcount + 1
            qcount = qcount + 1
            rcount = rcount + 1
            outgt = record1.REF + record2.REF + record3.REF
            totcalls = totcalls + 1        
        haplist.append(outgt)
    # assign allele freq values
    values = {}
    if totcalls == 0:
        values['p1'] = 0
        values['q1'] = 0
        values['r1'] = 0
    else:
        values['p1'] = pcount/totcalls
        values['q1'] = qcount/totcalls
        values['r1'] = rcount/totcalls
    values['p2'] = 1 - values['p1']
    values['q2'] = 1 - values['q1']
    values['r2'] = 1 - values['r1']
    # haplotype frequencies
    uniques = dict.fromkeys(set(haplist))
    for hap in uniques:
        uniques[hap] = haplist.count(hap)/len(haplist)
    # create tuples w/ actual haps and corresponding AB notation
    homrefC = record1.REF + record2.REF + record3.REF, 'ABC'
    homrefc = record1.REF + record2.REF + record3.REF, 'ABc'
    homaltC = str(record1.ALT[0]) + str(record2.ALT[0]) + str(record3.ALT[0]), 'abC'
    homaltc = str(record1.ALT[0]) + str(record2.ALT[0]) + str(record3.ALT[0]), 'abc'
    het1c = record1.REF + str(record2.ALT[0]) + str(record3.ALT[0]), 'Abc'
    het2c = str(record1.ALT[0]) + record2.REF + str(record3.ALT[0]), 'aBc'
    het1C = record1.REF + str(record2.ALT[0]) + record3.REF, 'AbC'
    het2C = str(record1.ALT[0]) + record2.REF + record3.REF, 'aBC'
    haps = {}
    # use tuples to assign hap freqs to AB notation for downstream use
    genotypes = [homrefC, homrefc, homaltC, homaltc, het1c, het2c, het1C, het2C]
    for genotype in genotypes:
        if genotype[0] in uniques.keys():
            haps[genotype[1]] = uniques[genotype[0]]
        else:
            haps[genotype[1]] = 0

    return uniques, values, haps


def triple_dcalc(record1, record2, record3):
    # Dabc = Fabc - p1*Dbc - q1*Dac - r1*Dab - p1*q1*r1
    
    from popgen import freqsgetter
    
    def quick_dcalc(record1, record2):
        haps = freqsgetter(record1, record2, snpcheck = False)[2]
        try:
            LHS = haps['AB'] * haps['ab']
        except KeyError: # either hap missing
            LHS = 0
        try:
            RHS = haps['Ab'] * haps['aB']
        except KeyError:
            RHS = 0
        d = LHS - RHS
        # d = round(d, 5)
        return d

    freqs = triplefreqsgetter(record1, record2, record3)[1]
    haps = triplefreqsgetter(record1, record2, record3)[2]

    try:
        Fabc = haps['ABC']
    except KeyError:
        Fabc = 0

    second_term = freqs['p1'] * quick_dcalc(record2, record3)
    third_term = freqs['q1'] * quick_dcalc(record1, record3)
    fourth_term = freqs['r1'] * quick_dcalc(record1, record2)
    fifth_term = freqs['p1'] * freqs['q1'] * freqs['r1']

    dabc = Fabc - second_term - third_term - fourth_term - fifth_term

    # dabc = round(dabc, 5)
    return dabc
