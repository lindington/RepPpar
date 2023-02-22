### EMU-PCA
Because of the large amount of missing data in some individuals, I tried to additionally make an EMU-PCA, which imputes missing data better for extremely large amounts, designed especially for target/capture data.
However, EMU-PCA uses plink style ``.bed`` files and cannot be performed on genotype likelihood files (e.g. ``.beagle`` extention gl files). Since there is no pipeline for using ANGSD files for emu-pca input, I tried different workarounds.

1. I used the beta option ``-doBcf 1`` in the same script that i produced the ``.beagle`` file with, to see if i can produce an adequate VCF file for downstream analyses. This then has to be converted to plink input format.

2. Beagle has some utilities i can use in combination to produce a vcf file with this pipeline: 
   - ``gprobs2beagle.jar`` which calls genotypes into a ``.beagle`` file
   - ``beagle2vcf.jar`` which converts the genotype file into a ``.vcf`` file
   - plink has an option to convert ``.vcf`` files into plink style ``.bed`` files.

3. I can check the options for using ANGSD directly to call genotypes using ``-doGeno 2`` (from ANGSD: ``"2: write the called genotype encoded as -1,0,1,2, -1=not called"``)
 
> skipped, too much effort 
> and overkill