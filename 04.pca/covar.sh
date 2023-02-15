#gzip -d .geno file

#count the number of variable sites
#zcat ../03.beagle/01.output/219_out_maxdepth.mafs.gz | tail -n+2 | wc -l
#371665 variable sites

/dss/dsshome1/lxc0A/ru67vil/programs/ngsTools/ngsPopGen/ngsCovar -probfile ../03.beagle/01.output/219_out_maxdepth.geno -outfile 219_out_maxdepth.covar -nind 219 -nsites 371665 -call 0 -norm 0

