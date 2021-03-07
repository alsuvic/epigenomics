# Epigenomics course

## 3.1 EN‐TEx ChIP‐seq data: how to navigate the portal and run the chipnf pipeline

We run the container
```
sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
```

We clone the repository
```
git clone https://github.com/bborsari/epigenomics_uvic
cd epigenomics_uvic
cd ChIP-seq; ls
```

We download the metadata file.
```
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&assembly=GRCh38&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assay_slims=DNA+binding" 
```

We explore the structure of the metadata file
```
head -1 metadata.tsv
```

We retrieve the FASTQ IDs
```
head -1 metadata.tsv | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++){print $i, i}}'
```

We check how many FASTQ files are available for this experiment.
```
grep -F H3K4me3 metadata.tsv | grep -F sigmoid_colon | awk 'BEGIN{FS="\t"}$2=="fastq"{n++}END{print n}'
```

As we can see, there are 4 files, so we download them.
```
grep -F H3K4me3 metadata.tsv |\
grep -F sigmoid_colon |\
awk 'BEGIN{FS=OFS="\t"}$2=="fastq"{print $1}' |\
while read filename; do 
 wget -P data/fastq.files "https://www.encodeproject.org/files/$filename/@@download/$filename.fastq.gz";
done
```

We download control samples' FASTQ files. 
```
echo -e "ENCFF102SFU\nENCFF599MFK\nENCFF187FWF\nENCFF549PVM\nENCFF876EUR\nENCFF950GED" |\
while read filename; do 
 wget -P data/fastq.files "https://www.encodeproject.org/files/$filename/@@download/$filename.fastq.gz";
done
```

```
mypath=$(pwd)
```

We prepare two FASTQ files for the sample (H3K4me3 ChIP):
```
grep -F H3K4me3 metadata.tsv |\
grep -F sigmoid_colon |\
awk -v mypath="$mypath" 'BEGIN{FS=OFS="\t";n=0}$2=="fastq"{n++;split($22, a, "-"); print $10"_"a[1]"_"$34, $10"_"a[1]"_"$34"_run"n, mypath"/data/fastq.files/"$1".sub.fastq", "sigmoid_colon_control_1", a[1]}' |\
head -2 > chip-nf/pipeline.index.tsv
```

and the two FASTQ files for the control:
```
echo -e "ENCFF102SFU\nENCFF599MFK\nENCFF187FWF\nENCFF549PVM\nENCFF876EUR\nENCFF950GED" |\
head -2 |\
awk -v mypath="$mypath" 'BEGIN{FS=OFS="\t";n=0}{n++; print "sigmoid_colon_control_1", "sigmoid_colon_control_1_run"n, mypath"/data/fastq.files/"$1".sub.fastq", "sigmoid_colon_control_1", "input"}' >> chip-nf/pipeline.index.tsv
```

Now, we leave the container, install Nextflow and run the pipeline.
```
exit
cd
curl -fsSL get.nextflow.io | bash
sudo ./nextflow run guigolab/chip-nf -bg -r v0.2.3 --index epigenomics/epigenomics_uvic/ChIP-seq/chip-nf/pipeline.index.tsv --genome epigenomics/epigenomics_uvic/ChIP-seq/chip-nf/refs/GRCh38.primary_assembly.genome.chr19.fa --genomeSize hs -with-docker > pipeline.log.txt
```

When it finishes, we move the files.
```
sudo cp pipeline.log.txt epigenomics/epigenomics_uvic/ChIP-seq/chip-nf/pipeline.log.txt
sudo cp chipseq-pipeline.db epigenomics/epigenomics_uvic/ChIP-seq/chip-nf/chipseq-pipeline.db
```

We reactivate the docker container
```
cd epigenomics
sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
cd epigenomics_uvic/ChIP-seq
```

We can see the log file
```
cat chip-nf/pipeline.log.txt
```

and the output file.
```
column -t chip-nf/chipseq-pipeline.db 
```

## 3.2. EN‐TEx ChIP‐seq data: downstream analyses

We create a folder analyses.
```
mkdir analyses
```

Then, we prepare two folders to store bigBed peak calling and bigWig FC signal files.
```
grep -F H3K4me3 metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $10, $22}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.txt
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done

grep -F H3K4me3 metadata.tsv |\
grep -F "bigWig" |\
grep -F "fold_change_over_control" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $10, $22}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigWig.FC.ids.txt
cut -f1 analyses/bigWig.FC.ids.txt |\
while read filename; do
  wget -P data/bigWig.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigWig"
done
```

We check their integrity.
```
for file_type in bigBed bigWig; do

  # retrieve original MD5 hash from the metadata
  ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,45 > data/"$file_type".files/md5sum.txt

  # compute MD5 hash on the downloaded files 
  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  # make sure there are no files for which original and computed MD5 hashes differ
  awk '$2!=$3' data/"$file_type".files/md5sum.txt

done
```

We create a folder, download the  Gencode reference file and decompress it.
```
mkdir annotation
wget https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz
gunzip annotation/gencode.v24.primary_assembly.annotation.gtf.gz
```

We convert the gtf annotation file to a BED format.
```
awk '$3=="gene"' annotation/gencode.v24.primary_assembly.annotation.gtf |\
grep -F "protein_coding" |\
cut -d ";" -f1 |\
awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $10, 0, $7, $10}' |\
sed 's/\"//g' |\
awk 'BEGIN{FS=OFS="\t"}$1!="chrM"{$2=($2-1); print $0}' > annotation/gencode.v24.protein.coding.gene.body.bed
```

We prepare a text file with the IDs for the gene expression matrices for total RNA-seq.
```
echo -e "ENCFF268RWA\tsigmoid_colon\nENCFF918KPC\tstomach" > analyses/tsv.totalRNASeq.ids.txt
```

We create a folder tsv.files inside data, and download there the expression matrices
```
mkdir data/tsv.files
cut -f1 analyses/tsv.totalRNASeq.ids.txt |\
while read filename; do
  wget -P data/tsv.files "https://www.encodeproject.org/files/$filename/@@download/$filename.tsv"
done
```

We subset protein-coding genes.
```
cut -f1 analyses/tsv.totalRNASeq.ids.txt |\
while read filename; do 
  ../bin/selectRows.sh <(cut -f4 annotation/gencode.v24.protein.coding.gene.body.bed) <(cut -f1,6 data/tsv.files/"$filename".tsv) > tmp
  mv tmp data/tsv.files/"$filename".tsv
done
```

We select the 1000 most expressed genes in each of the two tissues.
```
cat analyses/tsv.totalRNASeq.ids.txt |\
while read filename tissue; do
  sort -k2,2gr data/tsv.files/"$filename".tsv |\
  head -1000 |\
  cut -f1 > analyses/"$tissue".1000.most.expressed.genes.txt
done
```

We select the 1000 least expressed genes in each of the two tissues.
```
cat analyses/tsv.totalRNASeq.ids.txt |\
while read filename tissue; do 
  sort -k2,2gr data/tsv.files/"$filename".tsv |\
  tail -1000 |\
  cut -f1 > analyses/"$tissue".1000.least.expressed.genes.txt
done
```

We prepare BED files for the 1000 least expressed genes in the two tissues:
```
for tissue in stomach sigmoid_colon; do
  ../bin/selectRows.sh analyses/"$tissue".1000.least.expressed.genes.txt <(awk 'BEGIN{FS=OFS="\t"}{print $4, $0}' annotation/gencode.v24.protein.coding.gene.body.bed) |\
  cut -f2- > annotation/"$tissue".1000.least.expressed.genes.bed
done
```

We prepare BED files for the 1000 most expressed genes in the two tissues:
```
for tissue in stomach sigmoid_colon; do
  ../bin/selectRows.sh analyses/"$tissue".1000.most.expressed.genes.txt <(awk 'BEGIN{FS=OFS="\t"}{print $4, $0}' annotation/gencode.v24.protein.coding.gene.body.bed) |\
  cut -f2- > annotation/"$tissue".1000.most.expressed.genes.bed
done
```

We create a folder to store the aggregated signal over the TSS of the selected genes.
```
cd analyses/
mkdir aggregation.plot
```

We run bwtool aggregate for both samples.
```
bwtool aggregate 2000:2000 -starts -keep-bed annotation/stomach.1000.most.expressed.genes.bed data/bigWig.files/ENCFF391KDD.bigWig analyses/aggregation.plot/stomach.1000.most.expressed.genes.aggregate.tsv
bwtool aggregate 2000:2000 -starts -keep-bed annotation/stomach.1000.least.expressed.genes.bed data/bigWig.files/ENCFF391KDD.bigWig analyses/aggregation.plot/stomach.1000.least.expressed.genes.aggregate.tsv
bwtool aggregate 2000:2000 -starts -keep-bed annotation/sigmoid_colon.1000.most.expressed.genes.bed data/bigWig.files/ENCFF391KDD.bigWig analyses/aggregation.plot/sigmoid_colon.1000.most.expressed.genes.aggregate.tsv
bwtool aggregate 2000:2000 -starts -keep-bed annotation/sigmoid_colon.1000.least.expressed.genes.bed data/bigWig.files/ENCFF391KDD.bigWig analyses/aggregation.plot/sigmoid_colon.1000.least.expressed.genes.aggregate.tsv
```

We create an aggregation plot for each tissue.
```
for tissue in stomach sigmoid_colon; do 
  Rscript ../bin/aggregation.plot.R --most analyses/aggregation.plot/"$tissue".1000.most.expressed.genes.aggregate.tsv --least analyses/aggregation.plot/"$tissue".1000.least.expressed.genes.aggregate.tsv --tissue "$tissue" --output analyses/aggregation.plot/aggregation.plot."$tissue".pdf
done
```

We download the matrix of H3K4me3 FC signals over promoter regions.
```
wget https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/H3K4me3.matrix.tsv
```

We create a folder to store the study of the correlation coefficient between expression and H3K4me3 levels-
```
mkdir scatterplot.correlation
cd ..
```

We check that the row order (list of genes) is the same in the expression matrices of the two tissue.
```
diff <(cut -f1 data/tsv.files/ENCFF268RWA.tsv) <(cut -f1 data/tsv.files/ENCFF918KPC.tsv)
```

We get the expression matrix.
```
paste data/tsv.files/ENCFF268RWA.tsv <(cut -f2 data/tsv.files/ENCFF918KPC.tsv) |\
awk 'BEGIN{FS=OFS="\t"; print "sigmoid_colon", "stomach"}{print}' > analyses/expression.matrix.tsv
```

We check that the row order is the same for expression and H3K4me3 matrices.
```
diff <(cut -f1 analyses/H3K4me3.matrix.tsv) <(cut -f1 analyses/expression.matrix.tsv)
```

We produce a scatterplot of expression (x axis) vs. H3K4me3 (y axis).
```
for tissue in sigmoid_colon stomach; do
  Rscript ../bin/scatterplot.correlation.R --expression analyses/expression.matrix.tsv --mark analyses/H3K4me3.matrix.tsv --tissue "$tissue" --output analyses/scatterplot.correlation/scatterplot.correlation."$tissue".pdf
done
```

We create some folders.
```
cd analyses/
mkdir peaks.analysis
cd ..
mkdir data/bed.files
```

We convert bigBed files of H3K4me3 peaks to BED files.
```
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

We download the list of promoters ([-2 kb, +2 Kb] from TSS) of protein-coding genes.
```
cd annotation/
wget https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/gencode.v24.protein.coding.non.redundant.TSS.bed
cd ..
```

We retrieve genes with peaks of H3K4me3 at the promoter region in each tissue.
```
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -b data/bed.files/"$filename".bed -u |\
  cut -f7 |\
  sort -u > analyses/peaks.analysis/genes.with.peaks."$tissue".H3K4me3.txt
done
```

Genes marked by H3K4me3 in both tissues.
```
../bin/selectRows.sh analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt |\
cut -d "." -f1 > analyses/peaks.analysis/genes.marked.both.tissues.H3K4me3.txt
```

Genes with sigmoid colon-specific marking.
```
../bin/discardRows.sh analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt |\
cut -d "." -f1 > analyses/peaks.analysis/genes.with.sigmoid_colon.specific.peaks.H3K4me3.txt
```

Genes with stomach-specific marking.
```
../bin/discardRows.sh analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt |\
cut -d "." -f1 > analyses/peaks.analysis/genes.with.stomach.specific.peaks.H3K4me3.txt
```

Genes not marked in any of the two tissues.
```
../bin/discardRows.sh <(cat analyses/peaks.analysis/genes.marked.both.tissues.H3K4me3.txt analyses/peaks.analysis/genes.with.stomach.specific.peaks.H3K4me3.txt analyses/peaks.analysis/genes.with.sigmoid_colon.specific.peaks.H3K4me3.txt) <(cut -f7 annotation/gencode.v24.protein.coding.gene.body.bed |\
cut -d "." -f1) > analyses/peaks.analysis/genes.not.marked.H3K4me3.txt
```

GO enrichment analysis.
```
cut -f7 annotation/gencode.v24.protein.coding.gene.body.bed |\
cut -d "." -f1 > analyses/peaks.analysis/universe.genes.txt
```

We compare the distribution of expression values between the four sets of genes.
```
Rscript ../bin/boxplot.expression.R --expression analyses/expression.matrix.tsv --marked_both_tissues analyses/peaks.analysis/genes.marked.both.tissues.H3K4me3.txt --stomach_specific analyses/peaks.analysis/genes.with.stomach.specific.peaks.H3K4me3.txt --sigmoid_colon_specific analyses/peaks.analysis/genes.with.sigmoid_colon.specific.peaks.H3K4me3.txt --not_marked analyses/peaks.analysis/genes.not.marked.H3K4me3.txt --output analyses/peaks.analysis/boxplot.expression.pdf
```

We retrieve bigBed peak calling IDs for POLR2A from the metadata.
```
grep -F POLR2A-human metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_IDR_thresholded_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $10, $22}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u >> analyses/bigBed.peaks.ids.txt
```

We download the bigBed files.
```
awk '$3=="POLR2A-human"{print $1}' analyses/bigBed.peaks.ids.txt |\
while read filename; do 
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

We check the integrity of the downloaded bigBed files.
```
# Retrieve MD5 hashes of the files from the metadata
../bin/selectRows.sh <(awk '$3=="POLR2A-human"{print $1}' analyses/bigBed.peaks.ids.txt) metadata.tsv | cut -f1,45 > data/bigBed.files/tmp

# Compute MD5 hashes on the downloaded files
cat data/bigBed.files/tmp |\
while read filename original_md5sum; do 
  md5sum data/bigBed.files/"$filename".bigBed |\
  awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}'
done >> data/bigBed.files/md5sum.txt
rm data/bigBed.files/tmp

# Make sure there are no files for which original and computed MD5 hashes differ
awk '$2!=$3' data/bigBed.files/md5sum.txt
```

We convert the bigBed files to BED files.
```
awk '$3=="POLR2A-human"{print $1}' analyses/bigBed.peaks.ids.txt |\
while read filename; do 
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

Genes with peaks of POLR2A in each tissue.
```
grep -F POLR2A analyses/bigBed.peaks.ids.txt |\
cut -f-2 |\
while read filename tissue; do 
  bedtools intersect -a annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -b data/bed.files/"$filename".bed -u |\
  cut -f7 |\
  sort -u > analyses/peaks.analysis/genes.with.peaks."$tissue".POLR2A.txt
done
```

We make the Venn Diagram between the sets of genes with peaks of H3K4me3 and POLR2A in the two tissues.
```
Rscript ../bin/VennDiagram.4groups.R --setA analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt --setB analyses/peaks.analysis/genes.with.peaks.stomach.POLR2A.txt --setC analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt --setD analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.POLR2A.txt --output analyses/peaks.analysis/Venn.Diagram.H3K4me3.POLR2A.png
```

## 4. EN‐TEx ATAC‐seq data: downstream analyses

**Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files. Make sure the files are organized in a consistent way as done for ChIP-seq.**
```
cd ..
cd ATAC-seq
mkdir analyses
mkdir data
mkdir data/bigBed.files
```

**Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections. Make sure your md5sum values coincide with the ones provided by ENCODE.**

We get the link to download the metadate file from the experiments in ENCODE associated the following characteristics:

Assay type: DNA accessibility

Assay title: ATAC-seq

Status: released

Genome assembly: GRCh38

Biosample term name: stomach AND sigmoid colon

Then, we download the file.
```
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_slims=DNA+accessibility&assay_title=ATAC-seq&biosample_ontology.term_name=stomach&biosample_ontology.term_name=sigmoid+colon" 
```

Now, we get the bigBed narrow, pseudoreplicated peaks and assembly GRCh38 for stomach and sigmoid_colon.
```
cat metadata.tsv | grep -F "bigBed_narrowPeak" |grep -F "pseudoreplicated_peaks" |grep -F "GRCh38" |awk 'BEGIN{FS=OFS="\t"}{print $1, $10, $22}' |sort -k2,2 -k1,1r |sort -k2,2 -u > analyses/bigBed.peaks.ids.txt
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

We can check their integrity by verifying their MD5 hash. As we can see in the image, they are correct.
```
for file_type in bigBed; do

  # retrieve original MD5 hash from the metadata
  ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,45 > data/"$file_type".files/md5sum.txt

  # compute MD5 hash on the downloaded files 
  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  # make sure there are no files for which original and computed MD5 hashes differ
  awk '$2!=$3' data/"$file_type".files/md5sum.txt

done
cat data/bigBed.files/md5sum.txt
```

![image](https://user-images.githubusercontent.com/80123456/110252760-34ea1280-7f87-11eb-99eb-06a5eb921d75.png)

**For each tissue, run an intersection analysis using BEDTools: report 1) the number of peaks that intersect promoter regions, 2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions).**

BEDTools need BED files, so we have to convert bigBed files of peaks to BED files with the bigBedToBed command.
```
mkdir data/bed.files
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

First, we download the list of promoters ([-2 kb, +2 Kb] from TSS) of protein-coding genes.
```
mkdir annotation
cd annotation/
wget https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/gencode.v24.protein.coding.non.redundant.TSS.bed
cd ..
```

Then, we retrieve genes with peaks at the promoter region in each tissue.
```
mkdir analyses/peaks.analysis
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -b data/bed.files/"$filename".bed -u |\
  cut -f7 |\
  sort -u > analyses/peaks.analysis/genes.with.peaks."$tissue".txt
done
```

Now, we can see the new files.
```
head analyses/peaks.analysis/genes.with.peaks.stomach.txt
```

So we count the number of peaks that intersect promoter regions.
```
cat analyses/peaks.analysis/genes.with.peaks.stomach.txt | wc -l
cat analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.txt | wc -l
```

![image](https://user-images.githubusercontent.com/80123456/110247048-77e9bd00-7f6a-11eb-85dc-eb90bfc38c86.png)

As we can see, 15029 peaks intersect promoter regions in stomach and 14273 peaks intersect promoter regions in sigmoid colon.

For the second part, we need the genomic coordinates of the genes. Therefore, we download the Gencode annotation version 24.
```
cd annotation
wget https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz
cd ..
gunzip annotation/gencode.v24.primary_assembly.annotation.gtf.gz
```

We convert the gtf annotation file to a BED format in three steps:

retrieve gene body coordinates of protein-coding genes (chr, start, end, strand).

remove mitochondrial genes (i.e. those located on chrM).

move from a 1-based to a 0-based coordinate system.

```
awk '$3=="gene"' annotation/gencode.v24.primary_assembly.annotation.gtf |\
grep -F "protein_coding" |\
cut -d ";" -f1 |\
awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $10, 0, $7, $10}' |\
sed 's/\"//g' |\
awk 'BEGIN{FS=OFS="\t"}$1!="chrM"{$2=($2-1); print $0}' > annotation/gencode.v24.protein.coding.gene.body.bed
```

We retrieve the number of peaks that fall outside gene coordinates in each tissue.
```
cut -f-2 analyses/bigBed.peaks.ids.txt |while read filename tissue; do    bedtools intersect -a annotation/gencode.v24.protein.coding.gene.body.bed -b data/bed.files/"$filename".bed -v | sort -u > analyses/peaks.analysis/genes.with.peaks.outside."$tissue".txt; done

cat analyses/peaks.analysis/genes.with.peaks.outside.stomach.txt | wc -l
cat analyses/peaks.analysis/genes.with.peaks.outside.sigmoid_colon.txt | wc -l
```

![image](https://user-images.githubusercontent.com/80123456/110249995-5a701f80-7f79-11eb-8a15-214255e3fba3.png)

Therefore, 4836 peaks fall outside gene coordinates in sigmoid colon and 4036 peaks fall outside gene coordinates in stomach.

## 5. Distal regulatory activity

**Task 1: Create a folder regulatory_elements inside epigenomics_uvic. This will be the folder where you store all your subsequent results.**

```
cd ..
mkdir regulatory_elements
cd regulatory_elements/
```

**Task 2: Distal regulatory regions are usually found to be flanked by both H3K27ac and H3K4me1. From your starting catalogue of open regions in each tissue, select those that overlap peaks of H3K27ac AND H3K4me1 in the corresponding tissue. You will get a list of candidate distal regulatory elements for each tissue. How many are they?**

First, we create some folders.
```
mkdir analyses
mkdir data
mkdir data/bigBed.files
```

We download the metadata file using the script download.metadata.sh.
```
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&assembly=GRCh38&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assay_slims=DNA+binding" 
```

We will start with H3K27ac, so we first parse the metadata file to retrieve the corresponding IDs for H3K27ac, and then download the files in the forlder bigBed.files.
```
grep -F H3K27ac metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $10, $22}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.txt

cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

We can check their integrity by verifying their MD5 hash. As we can see in the image, they are correct.
```
for file_type in bigBed; do

  # retrieve original MD5 hash from the metadata
  ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,45 > data/"$file_type".files/md5sum.txt

  # compute MD5 hash on the downloaded files 
  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  # make sure there are no files for which original and computed MD5 hashes differ
  awk '$2!=$3' data/"$file_type".files/md5sum.txt

done
cat data/bigBed.files/md5sum.txt
```

![image](https://user-images.githubusercontent.com/80123456/110252867-c8bbde80-7f87-11eb-99d8-c990339f5d01.png)

We create some folders.
```
mkdir analyses/peaks.analysis
mkdir data/bed.files
mkdir annotation
```

We convert bigBed files of H3K27ac peaks to BED files with the bigBedToBed command
```
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

We convert the previous catalogue of open regions in each tissue to bed files in the current directory.
```
for tissue in stomach sigmoid_colon; do   cat ../ATAC-seq/analyses/peaks.analysis/genes.with.peaks.outside."$tissue".txt > annotation/open.outside.genes."$tissue".bed; done
```

We retrieve regions with peaks in these regions in each tissue.
```
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a annotation/open.outside.genes."$tissue".bed -b data/bed.files/"$filename".bed -u |\
  sort -u > analyses/peaks.analysis/open.regions.H3K27ac."$tissue".bed
done
```

We follow the same steps for the H3K4me1.
```
grep -F H3K4me1 metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $10, $22}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.txt

cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

```
for file_type in bigBed; do

  # retrieve original MD5 hash from the metadata
  ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,45 > data/"$file_type".files/md5sum.txt

  # compute MD5 hash on the downloaded files 
  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  # make sure there are no files for which original and computed MD5 hashes differ
  awk '$2!=$3' data/"$file_type".files/md5sum.txt

done
cat data/bigBed.files/md5sum.txt
```

![image](https://user-images.githubusercontent.com/80123456/110252933-0f113d80-7f88-11eb-914d-b750531e4a56.png)

Their MD5 hashes are also correct.

```
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

The difference is now we retrieve regions that overlaping peaks of H3K27ac in the corresponding tissue.
```
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a analyses/peaks.analysis/open.regions.H3K27ac."$tissue".bed -b data/bed.files/"$filename".bed -u |\
  sort -u > analyses/peaks.analysis/candidate.enhancers."$tissue".txt
done
```

![image](https://user-images.githubusercontent.com/80123456/110253032-762ef200-7f88-11eb-886d-7c31157cf875.png)

As we can see in the image, we got 209 candidate enhancers in stomach and 321 in sigmoid colon.

**Task 3: Focus on regulatory elements that are located on chromosome 1 (hint: to parse a file based on the value of a specific column, have a look at what we did here), and generate a file regulatory.elements.starts.tsv that contains the name of the regulatory region (i.e. the name of the original ATAC-seq peak) and the start (5') coordinate of the region.**

We filter by the chromosome 1 and keep only the name of the regulatory region and the start (5') coordinate of the region. So we have 106 candidates.
```
for tissue in stomach sigmoid_colon; do
  awk '$1=="chr1"'  analyses/peaks.analysis/candidate.enhancers."$tissue".txt | awk 'BEGIN{FS=OFS="\t"}{if ($6=="+"){start=$2} else {start=$3}; print $4, start}' >> regulatory.elements.starts.tsv
done
head regulatory.elements.starts.tsv
cat regulatory.elements.starts.tsv | wc -l
```

![image](https://user-images.githubusercontent.com/80123456/110258205-b058bd80-7fa1-11eb-9970-20d8abb00525.png)

**Task 4: Focus on protein-coding genes located on chromosome 1. From the BED file of gene body coordinates that you generated here, prepare a tab-separated file called gene.starts.tsv which will store the name of the gene in the first column, and the start coordinate of the gene on the second column (REMEMBER: for genes located on the minus strand, the start coordinate will be at the 3').**

We copy the file to the current directory.
```
cp ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed annotation/
```

We get the file.
```
cat annotation/gencode.v24.protein.coding.gene.body.bed | awk 'BEGIN{FS=OFS="\t"}{if ($6=="+"){start=$2} else {start=$3}; print $4, start}' > gene.starts.tsv
```

**Task 5: Download or copy this python script inside the epigenomics_uvic/bin folder. This script takes as input two distinct arguments: 1) --input corresponds to the file gene.starts.tsv (i.e. the file you generated in Task #4); 2) --start corresponds to the 5' coordinate of a regulatory element. Complete the python script so that for a given coordinate --start the script returns the closest gene, the start of the gene and the distance of the regulatory element.**

We download the script.
```
wget https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/get.distance.py
```

We modify the script.
```
nano get.distance.py
cat get.distance.py
```

![image](https://user-images.githubusercontent.com/80123456/110257074-733dfc80-7f9c-11eb-99be-2a25a29ab78c.png)

**Task 6. For each regulatory element contained in the file regulatory.elements.starts.tsv, retrieve the closest gene and the distance to the closest gene using the python script you created above.**

We get the closest gene and the distance to the closest gene.
```
cat regulatory.elements.starts.tsv | while read element start; do 
   python get.distance.py --input gene.starts.tsv --start $start; 
done > regulatoryElements.genes.distances.tsv
```

**Task 7: Use R to compute the mean and the median of the distances stored in regulatoryElements.genes.distances.tsv.**












