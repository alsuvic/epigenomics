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








