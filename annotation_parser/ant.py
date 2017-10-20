'''
Early draft of a parser (modeled heavily after pyVCF) to read in a custom annotation table.

Uses generator objects + lazy loading of records.

AH - 10/2017
'''

import gzip

class _Record(object):
    '''A record object - stores all information at a single row in the annotation table.
    '''
    def __init__(self, chromosome, position, reference_base, genic, exonic, intronic, intergenic, utr5,
        utr3, fold0, fold4, fold2, fold3, CDS, mRNA, rRNA, tRNA, feature_names, feature_types,
        feature_ID, cds_position, strand, frame, codon, aa, degen, FPKM, rho, FAIRE, recombination):

        def type_make(ob, type):
            if type == 'bool':
                if ob == 0:
                    return False
                elif ob == 1:
                    return True
                elif ob == '.':
                    return 'NA'
            elif type == 'int':
                try:
                    ob = int(ob)
                except ValueError:
                    ob = '.'
                return ob
            elif type == 'float':
                try:
                    ob = float(ob)
                except ValueError:
                    ob = '.'
                return ob
            elif type == 'str':
                try:
                    ob = str(ob)
                except:
                    ob = '.'
                return ob

        self.chrom = type_make(chromosome, 'str')
        self.pos = type_make(position, 'int')
        self.ref = type_make(reference_base, 'str')
        self.is_genic = type_make(genic, 'bool')
        self.is_exonic = type_make(exonic, 'bool')
        self.is_intronic = type_make(intronic, 'bool')
        self.is_intergenic = type_make(intergenic, 'bool')
        self.is_utr5 = type_make(utr5, 'bool')
        self.is_utr3 = type_make(utr3, 'bool')
        self.is_fold0 = type_make(fold0, 'bool')
        self.is_fold4 = type_make(fold0, 'bool')
        self.is_fold2 = type_make(fold0, 'bool')
        self.is_fold3 = type_make(fold0, 'bool')
        self.is_in_CDS = type_make(CDS, 'bool')
        self.is_in_mRNA = type_make(mRNA, 'bool')
        self.is_rRNA = type_make(rRNA, 'bool')
        self.is_tRNA = type_make(tRNA, 'bool')
        self.feature_names = list(feature_names)
        self.feature_types = list(feature_types)
        self.feature_ID = type_make(feature_ID, 'str')
        self.cds_position = type_make(cds_position, 'int')
        self.strand = type_make(strand, 'str')
        self.frame = type_make(frame, 'int')
        self.codon = type_make(codon, 'str')
        self.aa = type_make(aa, 'bool')
        self.degeneracy = type_make(degen, 'int')
        self.FPKM = type_make(FPKM, 'float')
        self.rho = type_make(rho, 'float')
        self.FAIRE = type_make(FAIRE, 'float')
        self.map_rho = type_make(recombination, 'float')


class Reader(object):
    '''The actual parser.
    
    Usage: 
    import ant
    parser = ant.Reader([annotation table filename])
    records = [r for r in parser.reader] # make sure to have this second .reader attribute to access generator
    
    Tabix compatibility unavailable at present, but actively being worked on.
    '''
    def __init__(self, filename = None, compressed = None):
        
        super(Reader, self).__init__

        if not filename:
            raise Exception('Error: filename not provided.')
        elif filename:
            if compressed is None:
                compressed = filename.endswith('.gz')
            self._reader = open(filename, 'rb' if compressed else 'rt')
        self.filename = filename
        if compressed:
            self._reader = gzip.GzipFile(fileobj = self._reader)

        self.reader = (line for line in self._reader) # init generator

        # 'burn' header from generator + set aside if user needs
        line = next(self.reader)
        header = []
        header.append(line)
        while line.startswith('##'):
            line = next(self.reader)
            header.append(line)

        assert line.startswith('#chromosome') # make sure header has been read in

        def _line_to_rec(line):
            '''Converts lines in annotation table to Record objects.
            Helper function to ensure easy fetching of actual Records and not just
            tab-split lines.
            '''
            row = line.rstrip().split('\t')
            assert len(row) == 30

            chromosome = row[0]
            position = row[1]
            reference_base = row[2]
            genic = row[3]
            exonic = row[4]
            intronic = row[5]
            intergenic = row[6]
            utr5 = row[7]
            utr3 = row[8]
            fold0 = row[9]
            fold4 = row[10]
            fold2 = row[11]
            fold3 = row[12]
            CDS = row[13]
            mRNA = row[14]
            rRNA = row[15]
            tRNA = row[16]
            feature_names = row[17]
            feature_types = row[18]
            feature_ID = row[19]
            cds_position = row[20]
            strand = row[21]
            frame = row[22]
            codon = row[23]
            aa = row[24]
            degen = row[25]
            FPKM = row[26]
            rho = row[27]
            FAIRE = row[28]
            recombination = row[29]

            record = _Record(chromosome, position, reference_base, genic, exonic, intronic, intergenic, utr5,
            utr3, fold0, fold4, fold2, fold3, CDS, mRNA, rRNA, tRNA, feature_names, feature_types,
            feature_ID, cds_position, strand, frame, codon, aa, degen, FPKM, rho, FAIRE, recombination)

            return record

        # generator without header
        self.reader = (_line_to_rec(line) for line in self.reader)

        self.cols = line.split('#')[1].split('\t') # get column names
        self.header = list(header)

    def __iter__(self):
        return self

    def metadata(self): # parser.metadata returns column names
        return self.cols

    def head(self): # parser.head returns entire header
        return self.header

    def next(self): # very similar to _line_to_rec above
        '''Return next record in file.'''
        line = next(self.reader)
        row = line.rstrip().split('\t')
        assert len(row) == 30

        chromosome = row[0]
        position = row[1]
        reference_base = row[2]
        genic = row[3]
        exonic = row[4]
        intronic = row[5]
        intergenic = row[6]
        utr5 = row[7]
        utr3 = row[8]
        fold0 = row[9]
        fold4 = row[10]
        fold2 = row[11]
        fold3 = row[12]
        CDS = row[13]
        mRNA = row[14]
        rRNA = row[15]
        tRNA = row[16]
        feature_names = row[17]
        feature_types = row[18]
        feature_ID = row[19]
        cds_position = row[20]
        strand = row[21]
        frame = row[22]
        codon = row[23]
        aa = row[24]
        degen = row[25]
        FPKM = row[26]
        rho = row[27]
        FAIRE = row[28]
        recombination = row[29]

        record = _Record(chromosome, position, reference_base, genic, exonic, intronic, intergenic, utr5,
        utr3, fold0, fold4, fold2, fold3, CDS, mRNA, rRNA, tRNA, feature_names, feature_types,
        feature_ID, cds_position, strand, frame, codon, aa, degen, FPKM, rho, FAIRE, recombination) # most args I've ever written...

        return record
