# vcf-filter
multi-hookers(key for filtering) &amp; multi-threads implementation for vcf filtering

## Prerequisite
* Python 3
* vcf
    * decomposed & normalized by [vt](https://genome.sph.umich.edu/wiki/Vt) or bcftools
    * annotated by annovar

## Installation
``` shell
git clone https://github.com/shanghungshih/vcf-filter.git
```

## Parameters
- `-v` `--vcfs` : the vcfs file which seperate by ','
- `-H` `--hookers` :  the information of the counters
- `-t` `--thread` :  pool size for multi-thread importing (default: 1)
- `--write2file` : `to be update`

## Notes
- configure the hookers.json before you run the program, and make sure the key of each hooker appear in your vcf annotation.
- multi-threads is for multiple vcfs


## Quick start
``` python
python3 vcf-filter.py -v sample1.hg19_multianno.vcf,sample2.hg19_multianno.vcf -H hookers.json
```
