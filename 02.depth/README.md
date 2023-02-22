# Depth
I extracted the number of raw reads per individual from the summary files using `grep -r 'seq_retained_reads' ./*.summary >> ../../capture/paleo_summary.txt`. 

#### Using ANGSD
I used the ANGSD option ``-doDepth`` to calculate coverage of individual samples to assess whether any individual had extremely low coverage. All the quality filters were kept the same as when making the ``.beagle`` and ``.geno`` files. I used the following script:

```bash
#!/bin/bash
#SBATCH -J depth_
#SBATCH --output=depthind.out
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=40
#SBATCH --mem=600gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH -t 4-00:00:00

STARTTIME=$(date +"%s")

module load angsd/0.933-gcc8

angsd -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist98.txt -ref /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -doDepth 1 -out 01.output/depth_ind -doCounts 1 -r chr1: -minMapQ 15 -minQ 20 -nThreads 4 -sites /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bait.sites -baq 1 -remove_bads 1 -uniqueOnly 1 -C 50 -only_proper_pairs 0  -minInd 79

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60, TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))

timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)

echo "Job ended at $(date). Took $timestamp hours:minutes:seconds to complete."
```

The output is given as individual and whole population depth. 
Coverage high in some genes, representing TEs?

>Output with mandatory maximum depth means no accurate representation of highest depth, or distribution of regions with high depth -> problematic because no full picture of highest depth.

- tried different max depths to make sure the last column is always 0 (higher depth gets counted as maxdepth gene)
- tried: 100 (fail, most regions have higher depth)
- tried: 1000000000 (fail, i.e. files unopenable with vim)
- tried: 1000000 (works, files openable with vim, last column always 0).

-> use maxvim = 1000000, but not elegant.

I adapted Vitalis scripts for quick plotting.

>-dodepth calculates perbase and you cant sum over regions. 
> continue with bedtools.