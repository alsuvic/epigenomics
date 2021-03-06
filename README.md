# Epigenomics course

## 3.1 

curl -fsSL get.nextflow.io | bash


/nextflow pull guigolab/chip-nf

nextflow run chip-nf --help
git clone https://github.com/guigolab/chip-nf.git


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

We can see the log file.
```
cat chip-nf/pipeline.log.txt
```

and the output file.
```
column -t chip-nf/chipseq-pipeline.db 
```














