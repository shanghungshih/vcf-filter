from multiprocessing import Pool
import argparse
from textwrap import dedent
import logging
import gzip
import json
import re
import os
import time

VERSION = (0, 1, 0)
__version__ = '.'.join(map(str, VERSION[0:3])) + ''.join(VERSION[3:])

class VCF_filter:
    def __init__(self, vcf, anchors, write2file, anchors_name):
        self.vcf = vcf
        self.anchors = anchors
        self.write2file = write2file
        self.anchors_name = anchors_name
        if self.write2file is True:
            if os.path.exists('passed-vcfs') is False:
                os.makedirs('passed-vcfs')
            vcf_name = self.vcf.split('/')[-1] if self.vcf[-2:] != 'gz' else self.vcf.split('/')[-1][:-2]
            self.output = open('passed-vcfs/%s-%s'%(self.anchors_name, vcf_name), 'w')
        self.results = {'total': 0, 'pass_anchors': 0}
        for anchor in self.anchors:
            if self.anchors[anchor]['type'] == 'in' or self.anchors[anchor]['type'] == 'not in':
                self.results['details_%s'%anchor] = {}
            self.results[anchor] = 0
            
        self.counter()
    
    def variant_generator(self):   # chr, start, end, ref, alt
        file = gzip.open(self.vcf) if self.vcf.endswith('.gz') else open(self.vcf, 'rb')
        for line in file:
            line = line.decode('utf-8')
            if line.startswith('#') is True:
                if self.write2file is True:
                    self.output.write(line)
            elif line.startswith('#') is False:
                variant = {}
                sep = line.strip('\n').split('\t')
                variant['chr'] = sep[0][3:] if sep[0].startswith('chr') else sep[0]
                variant['start'] = int(sep[1])
                variant['end'] = int(sep[1]) + int(len(sep[3])) - 1
                variant['ref'] = sep[3]
                variant['alt'] = sep[4]
                variant['variant'] = {'chr': variant['chr'], 'start': variant['start'], 'end': variant['end'], 'ref': variant['ref'], 'alt': variant['alt']}
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
                variant['line'] = line
                yield variant
    
    def counter(self):
        generator = self.variant_generator()
        for variant in generator:
            flag = True
            self.results['total'] += 1
            for anchor in self.anchors:
                text_details = variant[self.anchors[anchor]['key']] if type(variant[self.anchors[anchor]['key']]) is not dict else '%s %s %s %s %s'%(variant[self.anchors[anchor]['key']]['chr'], variant[self.anchors[anchor]['key']]['start'], variant[self.anchors[anchor]['key']]['end'], variant[self.anchors[anchor]['key']]['ref'], variant[self.anchors[anchor]['key']]['alt'])
                if self.anchors[anchor]['type'] == '==':
                    if variant[self.anchors[anchor]['key']] == self.anchors[anchor]['value']:
                        self.results[anchor] += 1
                    else:
                        flag = False
                elif self.anchors[anchor]['type'] == '>=':
                    if variant[self.anchors[anchor]['key']] >= self.anchors[anchor]['value']:
                        self.results[anchor] += 1
                    else:
                        flag = False
                elif self.anchors[anchor]['type'] == '<=':
                    if variant[self.anchors[anchor]['key']] <= self.anchors[anchor]['value']:
                        self.results[anchor] += 1
                    else:
                        flag = False
                elif self.anchors[anchor]['type'] == '>':
                    if variant[self.anchors[anchor]['key']] > self.anchors[anchor]['value']:
                        self.results[anchor] += 1
                    else:
                        flag = False
                elif self.anchors[anchor]['type'] == '<':
                    if variant[self.anchors[anchor]['key']] < self.anchors[anchor]['value']:
                        self.results[anchor] += 1
                    else:
                        flag = False
                elif self.anchors[anchor]['type'] == 'in':
                    if variant[self.anchors[anchor]['key']] in self.anchors[anchor]['value']:
                        self.results[anchor] += 1
                        if text_details not in self.results['details_%s'%anchor]:
                            self.results['details_%s'%anchor][text_details] = 1
                        else:
                            self.results['details_%s'%anchor][text_details] += 1
                    else:
                        flag = False
                elif self.anchors[anchor]['type'] == 'not in':
                    if variant[self.anchors[anchor]['key']] not in self.anchors[anchor]['value']:
                        self.results[anchor] += 1
                        if text_details not in self.results['details_%s'%anchor]:
                            self.results['details_%s'%anchor][text_details] = 1
                        else:
                            self.results['details_%s'%anchor][text_details] += 1
                    else:
                        flag = False
            if flag is True:
                self.results['pass_anchors'] += 1
                if self.write2file is True:
                    self.output.write(variant['line'])
            
class Queue_processor:
    def __init__(self, pool_size, vcfs, anchors, write2file):
        self.pool_size = pool_size
        self.vcfs = vcfs.strip().split(',')
        self.anchors = anchors
        self.write2file = write2file
        self.queue = []

        self.gen_queue()
        self.pool_handler()

    def vcf_handler(self, data):
        start = time.time()
        v = VCF_filter(vcf=data['vcf'], anchors=data['anchors'], write2file=data['write2file'], anchors_name=data['anchors_name'])
        self.queue[data['index']]['results'] = v.results
        end = time.time()
        logger_stdout.info('[%s]-[%s] (time used: %.2f min)' %(data['anchors_name'], data['vcf'], (end-start)/60))
        sorted_details = []
        
        for i in v.results:
            if i.startswith('details') is False:
                logger_stdout.info('[%s] : %s' %(i, v.results[i]))
            else:
                sorted_detail = {}
                sorted_detail['name'] = i
                sorted_detail['details'] = sorted(v.results[i], key=lambda x: x[1], reverse=False)
                sorted_details.append(sorted_detail)
        for i in sorted_details:
            logger_stdout.info('\t- %s:' %(i['name']))
            for j in i['details']:
                logger_stdout.info('\t[%s]: %s' %(j, v.results[i['name']][j]))
        logger_stdout.info('')

    def gen_queue(self):
        ct = 0
        for anchor in self.anchors:
            for vcf in self.vcfs:
                self.queue.append({'vcf': vcf, 'anchors_name': anchor, 'anchors': self.anchors[anchor], 'index': ct, 'results': None, 'write2file': self.write2file})
                ct += 1

    def pool_handler(self):
        pool = Pool(self.pool_size)
        pool.map(self.vcf_handler, self.queue)
        pool.close()
    
def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Unsupported value encountered.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
Testing environment: Python 3
Inputs:
1. vcf: Specify the vcf file with the -v or --vcf argument.
2. anchors: Specify the json file which contains the info of counters to count with the -h or --anchors argument.
3. write2file: Specify the -w or --write2file argument if you want to write to a new vcf.
4. thread: Specify the number of threads with the -t or --thread argument.

Usage:
1. without write2file
    python3 vcf-filter.py -v sample1.hg19_multianno.vcf,sample2.hg19_multianno.vcf -a multi-anchors/anchors-basic.json -t 2
2. with write2file
    python3 vcf-filter.py -w true -v sample1.hg19_multianno.vcf,sample2.hg19_multianno.vcf -a multi-anchors/anchors-PG-853variant.json
"""))
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument('-v', '--vcfs', type=str, help='the vcfs file which seperate by \',\'')
    required.add_argument('-a', '--anchors', type=str, help='the info of counters')
    optional.add_argument('-w', '--write2file', type=str2bool, default=False, help='write to new vcf')
    optional.add_argument('-t', '--thread', type=int, default=1, help='pool size for multi-thread importing (default: 1)')
    optional.add_argument('-V', '--version', action='version', version='%(prog)s ' + __version__)
    parser._action_groups.append(optional)
    args = parser.parse_args()

    with open(args.anchors) as f:
        anchors = json.loads(f.read())
    if os.path.exists('logs') is False:
        os.makedirs('logs')
    logging.basicConfig(level = logging.INFO, filename = "logs/%s.log" %(args.anchors.split('/')[-1].split('.')[0]))
    logger_stdout = logging.getLogger('stdout')
    logger_stdout.setLevel(logging.INFO)
    stdout_handler = logging.StreamHandler()
    stdout_handler.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger_stdout.addHandler(stdout_handler)
    
    logger_stdout.info('[vcf] : %s' %(args.vcfs.split(',')))
    logger_stdout.info('[anchors file] : [%s]' %(args.anchors))
    logger_stdout.info('[anchors] : %s' %([i for i in anchors.keys()]))
    logger_stdout.info('[write2file] : [%s]' %(args.write2file))
    logger_stdout.info('[threads] : [%s]' %(args.thread))
    
    if args.thread > 20:
        logger_stdout.info('[Terminated] Please specify a smaller number of [threads] <= 20.')
        exit()

    p = Queue_processor(pool_size=args.thread, vcfs=args.vcfs, anchors=anchors, write2file=args.write2file)
