# Call Nucleosomes

[![DOI](https://zenodo.org/badge/77785061.svg)](https://zenodo.org/badge/latestdoi/77785061)

Starting point - bed file
1. Create a wig file from bed file, extending only ±30 bp from fragment center. I use 120-160 bp fragments to make this wig file. Usage:
  
  perl bed2nucwig.pl \<text file with list of bed files> \<out file name> \<minimum fragment length> \<maximum fragment length>

2. Create a TSS file. I have included the mouse one I used (ori.tss_d.qlist). Here is the format:
  \<spacer text> \<chromosome> \<Gene Name> \<TSS> \<Strand> \<TES>

3. Call the peaks using gene_peak_call.pl :
  
  perl gene_peak_call.pl \<wig file generated in step 1> \<TSS file from step 2> \<JNK - spaceholder> \<bp window in which you want to call peaks 4500 for example will call peaks in TSS ± 4500 bp> > \<peak list file> 
