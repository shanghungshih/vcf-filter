# vcf-filter
multi-hookers (key for filtering, count if variant pass all hookers) &amp; multi-threads implementation for vcf filtering

## Prerequisite
* Python 3
* vcf
    * decomposed & normalized by [vt](https://genome.sph.umich.edu/wiki/Vt) or bcftools
    * annotated by annovar (.vcf, not .txt)

## Installation
``` shell
git clone https://github.com/shanghungshih/vcf-filter.git
```

## Parameters
- `-v` `--vcfs` : the vcfs file which seperate by ','
- `-H` `--hookers` :  the information of the counters
- `-t` `--thread` :  pool size for multi-thread importing (default: 1)
- `--write2file` : `to be update`

## Quick start
``` python
python3 vcf-filter.py -v sample1.hg19_multianno.vcf,sample2.hg19_multianno.vcf -H hookers.json
```

## hookers configuration
1. in hookers, define every `hooker name` which will show in results.
2. for each hooker, please define:
    - `key`: keys presents in info column of annovar-annotated vcf (ex. Func.refGene=TP53;AF=0.001;)
    - `type`: operator to perform comparison (valid types: `==`, `>=`, `<=`, `>`, `<`, `in`, `not in`)
    - `value`: operand to compare with vcf

```js
{
    "hookers": {
        "PASS": {
            "key": "PASS",
            "type": "==",
            "value": "PASS"
        },
        "AF<0.1": {
            "key": "AF",
            "type": "<",
            "value": 0.1
        },
        "Func.refGene": {
            "key": "Func.refGene",
            "type": "==",
            "value": "exonic"
        },
        "ExonicFunc.refGene": {
            "key": "ExonicFunc.refGene",
            "type": "==",
            "value": "nonsynonymous"
        },
        "ACMG_list": {
            "key": "Gene.refGene",
            "type": "in",
            "value": ["ACTA2", "ACTC1", "APC", "APOB", "ATP7B", "BMPR1A", "BRCA1", "BRCA2", "CACNA1S", "COL3A1", "DSC2", "DSG2", "DSP", "FBN1", "GLA", "KCNH2", "KCNQ1", "LDLR", "LMNA", "MEN1", "MLH1", "MSH2", "MSH6", "MUTYH", "MYBPC3", "MYH11", "MYH7", "MYL2", "MYL3", "NF2", "OTC", "PCSK9", "PKP2", "PMS2", "PRKAG2", "PTEN", "RB1", "RET", "RYR1", "RYR2", "SCN5A", "SDHAF2", "SDHB", "SDHC", "SDHD", "SMAD3", "SMAD4", "STK11", "TGFBR1", "TGFBR2", "TMEM43", "TNNI3", "TNNT2", "TP53", "TPM1", "TSC1", "TSC2", "VHL", "WT1"]
        }
    }
}
```

- output
    - `total`: # of total variants
    - `pass_hookers`: # of variants that pass all hookers
```
INFO     [vcf] : ['sample1.hg19_multianno.vcf', 'sample2.hg19_multianno.vcf']
INFO     [hookers file] : [hookers.json]
INFO     [hookers] : ['PASS', 'AF<0.1', 'Func.refGene', 'ExonicFunc.refGene', 'ACMG_list']
INFO     [write2file] : [False]
INFO     [threads] : [1]
INFO     {'total': 20, 'pass_hookers': 0, 'PASS': 0, 'AF<0.1': 0, 'Func.refGene': 1, 'ExonicFunc.refGene': 1, 'ACMG_list': 1}
INFO     done for [sample1.hg19_multianno.vcf] (0.00 min).
INFO     {'total': 20, 'pass_hookers': 0, 'PASS': 0, 'AF<0.1': 0, 'Func.refGene': 1, 'ExonicFunc.refGene': 1, 'ACMG_list': 1}
INFO     done for [sample2.hg19_multianno.vcf] (0.00 min).
```

## Notes
- count if variant pass all hookers
- for input file, only .vcf will be accepted.
- configure the hookers.json before you run the program, and make sure the key of each hooker appear in your vcf annotation.
- multi-threads is for multiple vcfs.
