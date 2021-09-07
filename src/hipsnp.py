import io
import os
import glob
import shutil
import subprocess
import requests
import pandas as pd
from datalad import api as datalad

def ensembl_human_rsid(rsid):
    """
    make a REST call to ensemble and return json info of a variant given a rsid
    """
    url = 'http://rest.ensembl.org/variation/human/' + rsid + '?content-type=application/json'
    response = requests.get(url)
    return response


def datalad_get_chromosome(c,
        source=None,
        imputationdir='imputation',
        path=None):
    """
    get a particular chromosome's (imputed) data
    """
    if source is None or source == '':
        source="ria+http://ukb.ds.inm7.de#~genetic"

    if path is None or path == '':
        path = os.path.join('/tmp', 'genetic')

    ds = datalad.clone(source=source, path=path)
    files = glob.glob(os.path.join(ds.path, imputationdir, '*_c' + str(c) + '_*'))
    ds.get(files)
    return files, ds


def rsid2chromosome(rsids):
    if isinstance(rsids, str) and os.path.isfile(rsids):
        rsids = pd.read_csv(rsids, header=None)
        rsids = list(rsids.iloc[:,0])
    elif isinstance(rsids, str):
        rsids = [rsids]

    chromosomes = [None] * len(rsids)
    for rs in range(len(rsids)):
        ens = ensembl_human_rsid(rsids[rs])
        ens = ens.json()
        ens = ens['mappings']
        for m in range(len(ens)):
            if ens[m]['ancestral_allele'] is not None:
                chromosomes[rs] = ens[m]['seq_region_name']

    df = pd.DataFrame()
    df['chromosomes'] = chromosomes
    df['rsids'] = rsids
    return df


def rsid2vcf(rsids, outdir,
        datalad_source="ria+http://ukb.ds.inm7.de#~genetic",
        qctool=None,
        datalad_drop=True,
        tmpdir='/tmp'):
    # check if qctool is available
    if qctool is None:
        qctool = shutil.which('qctool')

    if qctool is None:
        print('qctool is not available')
        raise

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if os.listdir(outdir):
        print('the output directory must be empty')
        raise

    # get chromosome of each rsid
    ch_rs = rsid2chromosome(rsids)
    uchromosomes = pd.unique(ch_rs['chromosomes'])
    print('chromosomes needed: ' + str(uchromosomes) + '\n')
    for c in range(len(uchromosomes)):
        ch = uchromosomes[c]
        ind = [i for i, x in enumerate(ch_rs['chromosomes']) if x == uchromosomes[c]]
        rs_ch = [rsids[i] for i in ind]
        print('chromosome ' + str(ch) + ' with ' + str(len(rs_ch)) + ' rsids\n')
        if len(rs_ch) < 11:
            print('rsids: ' + str(rs_ch) + '\n')
        # get the data
        print('datalad: getting files')
        files, ds = datalad_get_chromosome(ch, source=datalad_source)
        # find the bgen and sample files
        file_bgen = None
        file_sample = None
        for fl in files:
            name, ext = os.path.splitext(fl)
            if ext == '.bgen':
                assert file_bgen is None
                file_bgen = fl
            elif ext == '.sample':
                assert file_sample is None
                file_sample = fl

        assert file_bgen is not None and file_sample is not None
        file_rsids = os.path.join(outdir, 'rsids_chromosome' + str(ch) + '.txt')
        df = pd.DataFrame(rs_ch)
        df.to_csv(file_rsids, index=False, header=False)

        file_vcf = os.path.join(outdir, 'chromosome' + str(ch) + '.vcf')
        cmd = qctool + ' -g ' + file_bgen + ' -s ' + file_sample \
              + ' -incl-rsids ' + file_rsids  + ' -og ' + file_vcf
        print('running qctool: ' + cmd  + '\n')
        os.system(cmd)

        if datalad_drop:
            print('datalad: dropping files')
            ds.drop(files)

        print('done with chromosome ' + str(ch) + '\n')

    return ch_rs

def read_vcf(path):
    """
    taken shameless from: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
    Thanks Daichi Narushima
    """
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def vcf2genotype(vcf, th=0.9, snps=None, samples=None):
    """
    given a vcf file path or a pandas df from read_vcf returns genotypes
    """
    if isinstance(vcf ,str):
        vcf = read_vcf(vcf)
    elif isinstance(vcf, pd.DataFrame):
        pass
    else:
        print("don't know how to handle the input")
        raise

    format = pd.unique(vcf['FORMAT'])
    if len(format) != 1 or format[0] != 'GP':
        print('I can only deal with the GP format')
        raise

    nsnp = vcf.shape[0]
    ncol = vcf.shape[1]
    if samples is None:
        samples = [vcf.columns[i] for i in range(9, ncol)]
    else:
        assert all(sam in list(vcf.columns) for sam in samples)


    if snps is None:
        snps = list(vcf['ID'])
    else:
        assert all(snp in list(vcf['ID']) for snp in snps)

    labels = pd.DataFrame(index=range(len(snps)), columns=range(len(samples)))
    labels.index = snps
    labels.columns = samples
    snps_index = [snps.index(snp) for snp in snps]
    for snp in snps_index:
        REF = vcf['REF'][snp]
        ALT = vcf['ALT'][snp]
        for sam in samples:
            GP = vcf[sam][snp]
            GP = [float(x) for x in GP.split(',')]
            f = lambda i: GP[i]
            GT = max(range(len(GP)), key=f)
            if GP[GT] >= th:
                if GT == 0:
                    labels[sam][snp] = REF + REF
                elif GT == 1:
                    labels[sam][snp] = REF + ALT
                else:
                    labels[sam][snp] = ALT + ALT

    return labels

