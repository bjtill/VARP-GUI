# VARP-GUI
VCF Allele Ratio Plotter, Graphical Interface Version. A GUI tool to plot coverage from single or multi-sample variant call format (VCF) files.
_______________________________________________________________________________________________________________________________________________

Use at your own risk. I cannot provide support. All information obtained/inferred with this script is without any implied warranty of fitness for any purpose or use whatsoever.

ABOUT:  

Allele ratios can be calculated from the alleleic depth (AD) values reported in VCFs.  Deviations from expected allele ratios can be indicative of biologically significant phenomenon (e.g. loss of heterozygosity, aneuploidy) and also technical issues that may affect variant call accuracy (e.g. incorrect genome annotations, mapping errors).  VARP calculates allele ratios and plots ratios and coverage in several different ways to enable the user to visually identify potentially interesting deviations.  

INPUTS: 

1. A single or multi-sample VCF.  CAUTION: This file should be filtered for biallelic SNPs as the presence of more alleles may affect ratio calculations. 
2. A user-defined minimum number of total reads (coverage) as reported in the AD section of the VCF.
3. A user-defined choice to analyze all chromosomes or a subset of chromosomes. Using a subset is advised for genome annotations that contain many contigs.   
4. A user-defined choice to analyze all or a subset of samples in a multi-sample VCF.   


OUTPUTS:

1. Three data tables.  Ratio_Data.csv, contains information on the sample, chromosome, variant position, allele ratio calculated from AD, the coverage as the sum of AD, and a ratio bin letter code.  The following letters are assigned for different ranges of allele ratios:  

A: less than or equal to 20 percent reads supporting the reference allele (considered homozygous alternative as 80 percent or more reads support this call)

B: between 20 and  40 percent 

C: greater than or equal to 40 and less than or equal to 60 percent (considered an unambiguous heterozygous call)

D: greater than 60 and less than 80 percent 

E: greater than or euqal to 80 percent (considered homozygous reference)



This binning is applied to better visualize trends.  For example, excessive ratios in the B and D bins versus the C bin may indicate a change in ploidy, The file TotalCountRatioBin.csv contains the total count of each ratio bin. The file HetCountRatioBin.csv contains the count and percentage of allele ratio bins B,C,D.  This is to better distinguish variations in putative heterozygous calls.

2. Data plots of allele ratio bins on a per chromosome basis.  These are found in the directory PerChrom_AlleleRatioBin.  A separate plot is produced for each selected chromosome that contains data for all samples.

3. Data plots of allele ratio bins on a per sample basis.  These are found in the directory PerSample_AlleleRatioBin. A plot of all bins and a plot of only potentially heterozygous calls (bins B,C,D) are provided.  Data is plotted as a percentage of the total.

4. Data plots of allele ratios and coverage for each sample and each chromosome individually.  These are found in the directory PerSamp_PerChrom_Coverage_Ratio.  Depending on the size of the analysis, there can be many plots.  It may be best to look at the other plots first and then investigate any interesting samples/chromosomes at the individual level.  Coverage plots may aid in differentiating between, for example, changes in allele ratios due to copy number variation versus selective pressure.

5. VCFs created by the program are retained in the directory titled InputFiles.

6. A log file.  
 

REQUIREMENTS:  

Bash, YAD, Zenity, datamash, bcftools, awk, bgzip, tabix, R, ggplot2(R), cowplot(R)

TO RUN:

This program was built to run on Linux and tested on Ubuntu 20.04 and 22.04.  In theory it can be on macOS by installing the various dependencies (e.g. using Homebrew). However, I experienced issues with installing YAD.  Zenity installed okay, and so one could convert the YAD inputs to Zenity, or create a command line version.  No testing has been done with Windows and a Bash emulator.  

Download the .sh file and give it permission to run on your computer.  Open a Linux terminal and type chmod +x  VARP_V1_1.sh (or whatever the file is named).  Launch by typing ./VARP_V1_1.sh .  A window should appear where you can select input files and set various parameters. 

