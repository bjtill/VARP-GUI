#!/bin/bash 
#BT, July 8, 2024
#Version 1.1 adds organization of plots into directories. 
wget https://ucdavis.box.com/shared/static/4862en3ms2h2yj05a9tb99f2newuwwgb.jpeg
mv 4862en3ms2h2yj05a9tb99f2newuwwgb.jpeg VARPlogo.jpeg

YADINPUT=$(yad --width=1200 --title="VCF Allele Ratio Plotter (VARP)" --image=VARPlogo.jpeg --text="Version 1.1

ABOUT: A GUI tool to plot coverage and allele ratios from single or multi-sample VCFs.

CAUTION: Input VCF (single or multi-sample) must be filtered for bi-allelic SNPs. Presence of more alleles may affect ratio calculations. 
 
OUTPUTS: 
1) A table of all positions meeting the user-supplied allelic depth threshold.  This file, containing the name Ratio_Data.csv, contains information on the sample, chromosome, variant position, allele ratio calculated from AD, the coverage as the sum of AD, and a ratio bin letter code.  The following letters are assigned for different ranges of allele ratios:  

A: less than or equal to 20 percent reads supporting the reference allele (considered homozygous alternative as 80 percent or more reads support this call)
B: between 20 and  40 percent 
C: greater than or equal to 40 and less than or equal to 60 percent (considered an unambiguous heterozygous call)
D: greater than 60 and less than 80 percent 
E: greater than or euqal to 80 percent (considered homozygous reference)

2) Two additional data tables. The file TotalCountRatioBin.csv contains the total count of each ratio bin. The file HetCountRatioBin.csv contains the count and percentage of allele ratio bins B,C,D.  This is to better distinguish variations in putative heterozygous calls.  

3) Data plots of allele ratios by chromsome, sample, and both chromosome and sample.

DEPENDENCIES:  Bash, yad, Zenity, datamash, bcftools, awk, bgzip, tabix, R, ggplot2(R), cowplot(R)

VERSION INFORMATION: July 9, 2024 BT

LICENSE: MIT License, Copyright (c) 2024 Bradley John Till, see www.github.com/bjtill for full license information." --form --field="Your Initials for the log file" "Enter" --field="Optional notes for log file" "Enter" --field="Minimum allelic depth to retain for plotting (Click to edit):CBE" '20!5!10!50!100!150' --field="Percent of computer CPUs to use (Click to edit):CBE" '20!5!10!50!100' --field="Evaluate all chromosomes and contigs or a subset:CB" 'All!Subset' --field="Select the chromosome subset file (leave blank if evaluating all):FL" --field="Select the VCF file:FL" --field="Select a list of samples to subset the VCF. Leave blank if no sample subsetting required:FL" --field="Name for new directory. Your data will be in here. CAUTION-No spaces or symbols" )
echo $YADINPUT |  tr '|' '\t' | datamash transpose | head -n -1  > VARPparm1


#######################################################################################################################
#Check that the user provided the required information
a=$(awk 'NR==7 {if ($1=="") print "Missing"; else print "OK"}' VARPparm1)
cp VARPparm1 ${a}.vcfanswer
b=$(awk 'NR==9 {if ($1=="") print "Missing"; else print "OK"}' VARPparm1) 
cp VARPparm1 ${b}.diranswer
c=$(awk 'NR==5 {if ($1=="All") print "All"; else print "Subset"}' VARPparm1)
cp VARPparm1 ${c}.subanswer
d=$(awk 'NR==6 {if ($1=="") print "Missing"; else print "OK"}' VARPparm1)
cp VARPparm1 ${d}.sublist
if [ -f "Subset.subanswer" ] && [ -f "Missing.sublist" ] ;
then 
cp VARPparm1 NOSUBLIST
fi
if [ -f "Missing.vcfanswer" ] || [ -f "Missing.diranswer" ] || [ -f NOSUBLIST ] ; 
then

zenity --width 1200 --warning --text='<span font="32" foreground="red">VCF file, directory name and/or chromosome subset list not supplied. </span> \n Please close and start over.' --title="INFORMATION ENTRY FAILURE" 
rm *.vcfanswer VARPlogo.jpeg *.diranswer *.subanswer *.sublist NOSUBLIST
exit
fi 
rm *.vcfanswer VARPlogo.jpeg *.diranswer *.subanswer *.sublist NOSUBLIST

#######################################################################################################################
#Enter the directory
e=$(awk 'NR==9 {print $1}' VARPparm1)
mkdir ${e}
mv VARPparm1 ./${e}/
cd ${e}
#######################################################################################################################
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>VARPt.log 2>&1
now=$(date)  
echo "VCF Allele Ratio Plotter (VARP) Version 1.1
Script Started $now" 
(#Start
echo "# Checking if VCF file is compressed and indexed and taking these actions if it is not."; sleep 2

#Check if the vcf is compressed and if there is a bcftools index associated.  Create these if they don't exist. 
a=$(awk -F'.' 'NR==7 {if ($NF=="gz") print "YESGZ"; else print "NOGZ"}' VARPparm1)
awk 'NR==1 {print "answer"}' VARPparm1 > ${a}.answer1
if [ -f "NOGZ.answer1" ]; 
then
b=$(awk 'NR==7 {print $1}' VARPparm1)
bgzip $b 
#Assume no csi 
#percentage of CPUs
c=$(awk 'NR==7 {print $1".gz"}' VARPparm1)
e=$(awk 'NR==4 {print $1}' VARPparm1)
d=$(lscpu | grep "CPU(s):" | head -1 | awk '{print $2}' | awk -v var=$e '{print ($1*var)/100}' | awk '{printf "%.f\n", int($1+0.5)}')
bcftools index --threads $d $c
fi
a=$(awk 'NR==7 {print $1 }' VARPparm1 | sed 's/.gz//g' | awk '{print $1".gz.csi"}' )
 if [ ! -f $a ]; 
then
c=$(awk 'NR==7 {print $1 }' VARPparm1 | sed 's/.gz//g' | awk '{print $1".gz"}')
e=$(awk 'NR==4 {print $1}' VARPparm1)
d=$(lscpu | grep "CPU(s):" | head -1 | awk '{print $2}' | awk -v var=$e '{print ($1*var)/100}' | awk '{printf "%.f\n", int($1+0.5)}')
bcftools index --threads $d $c
fi 
echo "10"
echo "# Processing input VCF. This may take some time."; sleep 2 
#######################################################################################################################
#This section takes user inputs for chromosome and sample subsetting and takes appropriate action to create single-sample VCFs for later processing. 
A=$(awk 'NR==5 {if ($1=="All") print "All"; else print "Subset"}' VARPparm1)
cp VARPparm1 ${A}.subanswer

B=$(awk 'NR==8 {if ($1=="") print "No"; else print "Subset"}' VARPparm1)
cp VARPparm1 ${B}.sampanswer

if [ -f "Subset.subanswer"  ] && [ -f "Subset.sampanswer"  ] ;   
then 
C=$(awk 'NR==6 {print $1}' VARPparm1)
E=$(datamash transpose < $C | tr ' ' ',' | tr '\t' ',' )
D=$(awk 'NR==7 {print $1 }' VARPparm1 | sed 's/.gz//g' | awk '{print $1".gz"}')
g=$(awk 'NR==4 {print $1}' VARPparm1) #percent cpu
h=$(lscpu | grep "CPU(s):" | head -1 | awk 'NR==1{print $2}' | awk -v var=$g '{print ($1*var)/100}' | awk '{printf "%.f\n", int($1+0.5)}')
bcftools view --threads $h $D --regions ${E} -o ChromSub.vcf.gz
bcftools index --threads $h ChromSub.vcf.gz
e=$(awk 'NR==8 {print $1}' VARPparm1)
cp ${e} samples.txt
while IFS= read -r line; do
f=$(echo "$line" | awk '{print $1}') 
g=$(awk 'NR==4 {print $1}' VARPparm1) #percent cpu
h=$(lscpu | grep "CPU(s):" | head -1 | awk 'NR==1{print $2}' | awk -v var=$g '{print ($1*var)/100}' | awk '{printf "%.f\n", int($1+0.5)}')
bcftools view --threads $h -s $f ChromSub.vcf.gz -o ${f}.vcf
done < samples.txt
mkdir InputFiles
mv ChromSub* ./InputFiles/
rm samples.txt
fi 


if [ -f "Subset.subanswer"  ] && [ ! -f "Subset.sampanswer"  ] ;   
then 
C=$(awk 'NR==6 {print $1}' VARPparm1)
E=$(datamash transpose < $C | tr ' ' ',' | tr '\t' ',' )
D=$(awk 'NR==7 {print $1 }' VARPparm1 | sed 's/.gz//g' | awk '{print $1".gz"}')
g=$(awk 'NR==4 {print $1}' VARPparm1) #percent cpu
h=$(lscpu | grep "CPU(s):" | head -1 | awk 'NR==1{print $2}' | awk -v var=$g '{print ($1*var)/100}' | awk '{printf "%.f\n", int($1+0.5)}')
bcftools view --threads $h $D --regions ${E} -o ChromSub.vcf.gz
bcftools index --threads $h ChromSub.vcf.gz
bcftools query -l ChromSub.vcf.gz > samples.txt
while IFS= read -r line; do
f=$(echo "$line" | awk '{print $1}') 
g=$(awk 'NR==4 {print $1}' VARPparm1) #percent cpu
h=$(lscpu | grep "CPU(s):" | head -1 | awk 'NR==1{print $2}' | awk -v var=$g '{print ($1*var)/100}' | awk '{printf "%.f\n", int($1+0.5)}')
bcftools view --threads $h -s $f ChromSub.vcf.gz -o ${f}.vcf
done < samples.txt
mkdir InputFiles
mv ChromSub* ./InputFiles/
rm samples.txt
fi 

if [ ! -f "Subset.subanswer"  ] && [ -f "Subset.sampanswer"  ] ;   
then 
e=$(awk 'NR==8 {print $1}' VARPparm1)
cp ${e} samples.txt

while IFS= read -r line; do
f=$(echo "$line" | awk '{print $1}') 
D=$(awk 'NR==7 {print $1 }' VARPparm1 | sed 's/.gz//g' | awk '{print $1".gz"}')
g=$(awk 'NR==4 {print $1}' VARPparm1) #percent cpu
h=$(lscpu | grep "CPU(s):" | head -1 | awk 'NR==1{print $2}' | awk -v var=$g '{print ($1*var)/100}' | awk '{printf "%.f\n", int($1+0.5)}')
bcftools view --threads $h -s $f $D -o ${f}.vcf
done < samples.txt
rm samples.txt
mkdir InputFiles
fi 

if [ ! -f "Subset.subanswer"  ] && [ ! -f "Subset.sampanswer"  ] ;   
then 
D=$(awk 'NR==7 {print $1 }' VARPparm1 | sed 's/.gz//g' | awk '{print $1".gz"}')

bcftools query -l $D > samples.txt
while IFS= read -r line; do
f=$(echo "$line" | awk '{print $1}') 
D=$(awk 'NR==7 {print $1 }' VARPparm1 | sed 's/.gz//g' | awk '{print $1".gz"}')
g=$(awk 'NR==4 {print $1}' VARPparm1) #percent cpu
h=$(lscpu | grep "CPU(s):" | head -1 | awk 'NR==1{print $2}' | awk -v var=$g '{print ($1*var)/100}' | awk '{printf "%.f\n", int($1+0.5)}')
bcftools view --threads $h -s $f $D -o ${f}.vcf
done < samples.txt
rm samples.txt
mkdir InputFiles
fi 
echo "60"
echo "# Retrieving coverage and calculating allele ratios."; sleep 2 

#######################################################################################################################

ls *.vcf > vcflist
while IFS= read -r line; do
i=$(awk 'NR==3 {print $1}' VARPparm1) #depth
j=$(echo "$line" | awk '{print $1}')
bcftools query -f '%CHROM %POS [%AD\n]' $j | tr ',' ' '  | awk '{print $0, $3+$4}' | awk -v var=$i '{if ($5 >=var) print $1, $2, $3/$5, $5}' | awk -v var=$j '{if ($4 > 200) print var, $1, $2, $3, "200"; else print var, $0}' | sed 's/.vcf//g' | awk '{if ($4<=0.2) print $0, "A"; else if ($4>0.2 && $4<0.4) print $0, "B"; else if ($4>=0.4 && $4<=0.6) print $0, "C"; else if ($4>0.6 && $4<0.8) print $0, "D"; else if ($4>=0.8 && $4<=1) print $0, "E"}' >> tmpalldata  #downsamples to 200x, this could be adjusted further
done < vcflist 
#add line numbers- which I think was to control for R and then try to plot things out.  
awk '{print $1"_"$2, $0}' tmpalldata | awk '{print $4, $0}'  | awk 'BEGIN{print "P2", "Sample_Chrom", "Sample", "Chromosome", "Position", "Ratio", "Coverage", "Ratio_Bin"}1' | tr ' ' ',' > Ratio_Data.csv
#Position is repeated above because it is used in a match/replace 
#Plot by chromosome: re-generate a chrom list in case there was something bogus on the user list that could cause problems.  
awk '{print $2}' tmpalldata | awk '!visited[$1]++' > chrplotlist
awk '{print $1}' tmpalldata | awk '!visited[$1]++' > samplotlist

echo "70"
echo "# Creating per chromosome allele ratio bin plots."; sleep 2 
#create per chromosome allele ratio bin plots
while IFS= read -r line; do
k=$(echo "$line" | awk '{print $1}') #samplemarkerlist, note that grepping does not work to distinguish between F2BT and core, and so need to double nest or something
#issue: with few variants (like radseq) you get few positions and an unreadable figure. Therefore you want to add line numbers again with a list

head -1 Ratio_Data.csv > head

tail -n +2 Ratio_Data.csv | grep $k | awk -F,  'FNR==0{print $5,"linenumber";next} {print $5,FNR-1}' | tr ' ' ',' > tmpkey
tail -n +2 Ratio_Data.csv | grep $k | awk -F, 'NR==FNR{a[$1]=$2;next}{if (a[$1]) print a[$1], $2, $3, $4, $5, $6, $7, $8}' tmpkey - | tr ' ' ','  | cat head -  > tmpfp1.csv

printf 'library(ggplot2) \ng<-read.csv("tmpfp1.csv") \np <- ggplot(g, aes(x=P2, y=Sample, fill=factor(Ratio_Bin))) + geom_tile() + scale_fill_manual(values=c("A"= "blue", "B"="orange", "C"="green4", "D"="darkorange", "E"="darkblue", "F"="black")) + theme(axis.text.x = element_text(angle = 90, size =6, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 6)) + xlab("Position")\np2 <- p + labs(title= "Allele Ratio Bins For Chromosome: %s", fill="Allele Ratio Bin") + xlab ("Variant Number") \nggsave(plot = p2, filename="%s_AFbins.jpeg")' $k $k > AFbin.R
Rscript AFbin.R
done < chrplotlist
#move everything into separate dir to keep things clean
mkdir PerChrom_AlleleRatioBin
mv *.jpeg ./PerChrom_AlleleRatioBin/
rm tmpkey tmpfp1.csv AFbin.R

########################### 
#Per sample stats:  Try for a bar plot of the total counts of each ratio bin type.  The idea is that you care about the ratio of "true hets" to some mapping issue.  Could do this by chromosome as well
echo "80"
echo "# Creating per sample allele ratio bin plots."; sleep 2 
while IFS= read -r line; do
k=$(echo "$line" | awk '{print $1}') #samplemarkerlist, note that grepping does not work to distinguish between F2BT and core, and so need to double nest or something
#issue: with few variants (like radseq) you get few positions and an unreadable figure. Therefore you want to add line numbers again with a list
tail -n +2 Ratio_Data.csv | grep $k | awk -F, -v var=$k '{count[$8]++} END {for (word in count) print var, word, count[word]}' >> totalcouont
#and then go for just percentage of hets
tail -n +2 Ratio_Data.csv | grep $k | awk -F, '{if ($8=="B" || $8=="C" || $8=="D") print $0}' | awk -F, -v var=$k '{count[$8]++} END {for (word in count) print var, word, count[word]}' > tmp 
awk 'NR==FNR{a = a + $3;next} {c = ($3/a)*100;print $1,$2,$3,c }' tmp tmp >> hetcount
done < samplotlist

#all samples plot done outside of loop

awk 'BEGIN{print "Sample", "Bin", "VariantCount", "Percent"}1' hetcount | tr ' ' ',' > HetCountRatioBin.csv
awk 'BEGIN{print "Sample", "Bin", "Count"}1' totalcouont | tr ' ' ',' > TotalCountRatioBin.csv
printf 'library(ggplot2) \ng<-read.csv("TotalCountRatioBin.csv") \np <- ggplot(g, aes(x=Sample, y=Count, fill=Bin)) + geom_bar(position="stack", stat="identity") + scale_fill_manual(values=c("A"= "blue", "B"="orange", "C"="green4", "D"="darkorange", "E"="darkblue", "F"="black")) + theme(axis.text.x = element_text(angle = 90, size =6, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 6)) + xlab("Sample")\np2 <- p + labs(title= "Allele Ratio Bins by Sample and Count", fill="Allele Ratio Bin") \nggsave(plot = p2, filename="TotalCount_AlleleRatioBin.jpeg")' > t1.R
Rscript t1.R

printf 'library(ggplot2) \ng<-read.csv("HetCountRatioBin.csv") \np <- ggplot(g, aes(x=Sample, y=Percent, fill=Bin)) + geom_bar(position="stack", stat="identity") + scale_fill_manual(values=c("A"= "blue", "B"="orange", "C"="green4", "D"="darkorange", "E"="darkblue", "F"="black")) + theme(axis.text.x = element_text(angle = 90, size =6, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 6)) + xlab("Sample")\np2 <- p + labs(title= "Allele Ratio Bins by Sample and Count: Heterozygous", fill="Allele Ratio Bin") \nggsave(plot = p2, filename="HetCount_AlleleRatioBin.jpeg")' > t2.R
Rscript t2.R

## Rm *.R totalcouont hetcount tmp 
mkdir PerSample_AlleleRatioBin
mv *.jpeg ./PerSample_AlleleRatioBin/
echo "90"
echo "# Creating per sample and per chromosome coverage and allele ratio plots."; sleep 2 

###Now section where for each chromosome you have a coverage stacked with ratio plot: This probably ugly as well with limited data points.  
#a double nest as each sample and each chromosome get the treatment here
#start with sample then loop to chrom 
while IFS= read -r line; do
k=$(echo "$line" | awk '{print $1}') 
tail -n +2 Ratio_Data.csv | grep $k > tmp1
while IFS= read -r line; do

l=$(echo "$line" | awk '{print $1}') 
grep $l tmp1 | awk 'BEGIN{print "P2,Sample_Chrom,Sample,Chromosome,Position,Ratio,Coverage,Ratio_Bin"}1' > fplot.csv
printf 'library(ggplot2) \nlibrary(cowplot) \ng<-read.csv("fplot.csv") \np1 <-ggplot(g, aes(x = Position, y = Coverage)) + geom_line(linewidth=0.1, alpha=0.9, color= "darkorange") + theme(axis.text.x = element_text(angle = 90, size =3), axis.text.y = element_text(size = 3)) + labs(title= "Coverage Downsampled to 200, Sample %s, Chromosome %s") \np <- ggplot(g, aes(x = Position, y = Ratio)) +  geom_point(size=0.8, alpha=0.1, color= "deepskyblue3") + theme(axis.text.x = element_text(angle = 90, size =3), axis.text.y = element_text(size = 3)) + labs(title="Allele Ratios for Sample %s and Chromosome %s") \npa <- p + geom_hline(yintercept=0.20, linetype="dashed", alpha=0.5, color="firebrick") \npb <- pa + geom_hline(yintercept=0.40, linetype="dashed", alpha=0.5, color="firebrick") \npc <- pb + geom_hline(yintercept=0.60, linetype="dashed",alpha=0.5, color="firebrick") \npd <- pc + geom_hline(yintercept=0.80, linetype="dashed", alpha=0.5, color="firebrick") \np2 <- plot_grid(p1, pd, labels = "AUTO", ncol = 1)\nggsave(plot = p2, filename="CovRatio_%s_%s.jpeg", width=10, height=5, units=c("in"))' $k $l $k $l $k $l > bsplot.r
Rscript bsplot.r
done < chrplotlist
done < samplotlist
rm fplot.csv
mkdir PerSamp_PerChrom_Coverage_Ratio 
mv *.jpeg ./PerSamp_PerChrom_Coverage_Ratio/
echo "95"
echo "# Tidying."; sleep 2 
find . -maxdepth 1 -name "*.jpeg" -print0 | xargs -0 mv -t ./InputFiles/ 
rm *.r chrplotlist head hetcount samplotlist *.sampanswer *.subanswer tmp tmp1 tmpalldata totalcouont vcflist *.answer1 *.R
mkdir DataTables
c=$(date +"%m_%d_%y")
for i in *.csv; do 
mv $i ./DataTables/${i%.*}_${c}.csv; done 
cd PerChrom_AlleleRatioBin
c=$(date +"%m_%d_%y")
for i in *.jpeg; do 
mv $i ${i%.*}_${c}.jpeg; done
cd ..
cd PerSample_AlleleRatioBin
c=$(date +"%m_%d_%y")
for i in *.jpeg; do 
mv $i ${i%.*}_${c}.jpeg; done
cd ..
cd PerSamp_PerChrom_Coverage_Ratio
c=$(date +"%m_%d_%y")
for i in *.jpeg; do 
mv $i ${i%.*}_${c}.jpeg; done
cd ..
) | zenity --width 800 --title "VARP PROGRESS" --progress --auto-close
now=$(date)
echo "Script Finished" $now

printf 'Initials of person who ran the program: \nUser notes: \nUser-defined minimum allelic coverage: \nUser-defined percentage of CPU allocation: \nUser choice for subsetting by chromosomes: \nPath to chromosome list if subset chosen: \nPath to VCF file: \nPath to sample list if sample subsetting chosen: \nName of directory created for this analysis:' > plog
paste plog VARPparm1 | tr '\t' ' ' > plog2
c=$(date +"%m_%d_%y")
grep -v "Saving" VARPt.log | grep -v "geom_line" | grep -v "aesthetic" | grep -v "Warning" | grep -v "No shared levels" | grep -v "fill values." | cat - plog2 > VARP_${c}.log
rm plog plog2 VARPt.log rm VARPparm1
########END OF PROGRAM#################################################################################################
