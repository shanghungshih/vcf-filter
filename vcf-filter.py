from multiprocessing import Pool
import argparse
from textwrap import dedent
import logging
import gzip
import json
import re
import time

logger_stdout = logging.getLogger(__name__+'stdout')
logger_stdout.setLevel(logging.INFO)
stdout_handler = logging.StreamHandler()
stdout_handler.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
logger_stdout.addHandler(stdout_handler)

VERSION = (0, 1, 0)
__version__ = '.'.join(map(str, VERSION[0:3])) + ''.join(VERSION[3:])

class vcf_filter:
    def __init__(self, vcf, hookers):
        self.vcf = vcf
        self.hookers = hookers['hookers']
        self.results = {'total': 0}
        for hooker in self.hookers:
            self.results[hooker] = 0
        self.counter()
    
    def variant_generator(self):   # chr, start, end, ref, alt
        file = gzip.open(self.vcf) if self.vcf.endswith('.gz') else open(self.vcf, 'rb')
        for line in file:
            line = line.decode('utf-8')
            if line.startswith('#') is False:
                variant = {}
                sep = line.strip('\n').split('\t')
                variant['chr'] = sep[0][3:] if sep[0].startswith('chr') else sep[0]
                variant['start'] = int(sep[1])
                variant['end'] = int(sep[1]) + int(len(sep[3])) - 1
                variant['ref'] = sep[3]
                variant['alt'] = sep[4]
                variant['id'] = sep[2]
                variant['PASS'] = sep[6]
                annotations = sep[7].split(';')
                for annotation in annotations:
                    if '=' in annotation:
                        anno = annotation.split('=')
                        try:
                            variant[anno[0]] = float(anno[1])
                        except ValueError:
                            variant[anno[0]] = anno[1]
                variant['others'] = sep[5:]
                yield variant
    
    def counter(self):
        generator = self.variant_generator()
        for variant in generator:
            if self.results['total'] == 10000:
                break
            self.results['total'] += 1
            for hooker in self.hookers:
                try:
                    if self.hookers[hooker]['type'] == '==':
                        if variant[self.hookers[hooker]['key']] == self.hookers[hooker]['value']:
                            self.results[hooker] += 1
                    elif self.hookers[hooker]['type'] == '>=':
                        if variant[self.hookers[hooker]['key']] >= self.hookers[hooker]['value']:
                            self.results[hooker] += 1
                    elif self.hookers[hooker]['type'] == '<=':
                        if variant[self.hookers[hooker]['key']] <= self.hookers[hooker]['value']:
                            self.results[hooker] += 1
                    elif self.hookers[hooker]['type'] == '>':
                        if variant[self.hookers[hooker]['key']] > self.hookers[hooker]['value']:
                            self.results[hooker] += 1
                    elif self.hookers[hooker]['type'] == '<':
                        if variant[self.hookers[hooker]['key']] < self.hookers[hooker]['value']:
                            self.results[hooker] += 1
                    elif self.hookers[hooker]['type'] == 'in':
                        if variant[self.hookers[hooker]['key']] in self.hookers[hooker]['value']:
                            self.results[hooker] += 1
                    elif self.hookers[hooker]['type'] == 'not in':
                        if variant[self.hookers[hooker]['key']] not in self.hookers[hooker]['value']:
                            self.results[hooker] += 1
                except TypeError:
                    logger_stdout.info('ERROR [%s] contain >= two alt' %(variant))
                    exit()
                except KeyError:
                    logger_stdout.info('ERROR key [%s] for [%s] not in vcf' %(self.hookers[hooker]['key'], hooker))
                    exit()
            
def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Unsupported value encountered.')
        
def vcf_handler(data):
    start = time.time()
    v = vcf_filter(vcf=data[0], hookers=data[1])
    end = time.time()
    logger_stdout.info(v.results)
    logger_stdout.info('done for [%s] (%.2f min).' %(data[0], (end-start)/60))

def pool_handler(pool_size, vcfs, hookers, write2file):
    pool = Pool(pool_size)
    queue = [[vcf, hookers] for vcf in vcfs.strip().split(',')]
    pool.map(vcf_handler, queue)
    pool.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
Testing environment: Python 3
Inputs:
1. vcf: Specify the vcf file with the -v or --vcf argument.
2. hookers: Specify the json file which contains the info of counters to count with the -h or --hookers argument.
3. write2file: Specify the -w or --write2file argument if you want to write to a new vcf.
4. thread: Specify the number of threads with the -t or --thread argument.

Usage:
1. Import files with start and end into database
    python3 vcf_filter.py -v /home/nas210/TaiwanBiobank_data/joint-calling-vcf/taiwan_biobank_joint_calling_993_SNP_INDEL.recaled.decompose.normalize.vcf.gz -H hookers.json
"""))
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument('-v', '--vcfs', type=str, help='the vcfs file which seperate by \',\'')
    required.add_argument('-H', '--hookers', type=str, help='the info of counters')
    optional.add_argument('-w', '--write2file', type=str2bool, default=False, help='write to new vcf')
    optional.add_argument('-t', '--thread', type=int, default=1, help='pool size for multi-thread importing (default: 1)')
    optional.add_argument('-V', '--version', action='version', version='%(prog)s ' + __version__)
    parser._action_groups.append(optional)
    args = parser.parse_args()

    with open(args.hookers) as f:
        hookers = json.loads(f.read())
    
    logger_stdout.info('[vcf] : %s' %(args.vcfs))
    logger_stdout.info('[hookers file] : [%s]' %(args.hookers))
    logger_stdout.info('[hookers] : %s' %(hookers.keys()))
    logger_stdout.info('[write2file] : [%s]' %(args.write2file))
    logger_stdout.info('[threads] : [%s]' %(args.thread))
    
    if args.thread > 20:
        logger_stdout.info('[Terminated] Please specify a smaller number of [threads] <= 20.')
        exit()

    pool_handler(pool_size=args.thread, vcfs=args.vcfs, hookers=hookers, write2file=args.write2file)    