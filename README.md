# Acrocentric short arm mutations and recombination

Code and associated files used in analyze de novo mutations and recombination on acrocentric short arms. 
We recreated the assembly with verkko (v2.2.1) with extra HiFi and UL-ONT data for 23 samples in G2 (NA12877 and NA12878), G3 (NA12879, NA12881, NA12882, NA12883, NA12884, NA12885, NA12886, NA12887, 200080, 200100) and G4 (200081, 200082, 200084, 200085, 200086, 200087, 200101, 200102, 200103, 200104 and 200106).


## Used tools

moddotplot (v0.9.8), minimap2 (v2.28), rustybam (v0.1.33), seqtk (v1.4), wfmash (v0.13.0), deeptools(v3.5.5), mashmap (v3.1.3)

## Analysis

### Identification of pq-scatigs

The pipeline ```snakemakes/acro_parm.smk``` is used to identify pq-scatigs. It uses the configuration file ```snakemakes/liftover.yaml```. 

The configuration contains the following inputs:
1. Targeted regions are in ```snakemakes/roi.bed```. 
2. The manifest file ```master_asm.tab``` contains the assembly path.


This pipeline will create the following files:
1. ```results/{region}/fasta/{sample}_pq_contig.fa```: pq-scatig sequence in fasta format.
2. ```results/{region}/moddoplot/{sample}```: moddotplot for each pq-scatig

### Diversity

The all-vs-all alignments are created for the following comparisons with ```minimap2 -x asm20 -t 6 -c --eqx -D -P --dual=no all_seq.fa all_seq.fa > ava_mm2.paf```
1. The diversity of same chromosome among different samples in this family (e.g. Figure S4a).
2. Allelic and non-allelic comparison of chromosomes within individual genome (e.g. Figure S4b). This is only done for NA12877 and NA12878.

#### Allelic and non-allelic diversity

The alignments are divided into 10 kbp bins to calculate the sequence identity between allelic and non-allelic sequences. The results are used to create Figure 2c. 

```rb liftover --bed 10kb_bins.bed aln.paf | rb stats --paf > aln.stats```

### Transmissions and recombination

Align offspring to parents with wfmash (v0.13.0)
``
```wfmash -s 50k -l 150k -p 90 -n 1 -H 0.001 -t 10 offspring_qry.fa parent_ref.fa```

Parse alignment to find transmitted bases with ```transmission/AcroMain.py```. This script dose two things:
1. Find transmitted bases by identifying 1-to-1 homologous alignment (>=1Mbp and >=99% sequence identity) between child and parent. The results are used to create Figure 3a with function ```figure3a()``` in ```transmission/CreateFig.py```.
2. Identify potential recombination on child haplotypes. 


#### Recombination visualization

The recombination summary figure4a is created by function ```figure4a()``` in ```transmission/CreateFig.py```

The all-vs-all alignments are created for each transmission and recombination with ```minimap2 -x asm20 -t 6 -c --eqx -D -P --dual=no all_seq.fa all_seq.fa > ava_mm2.paf```.

#### Recombination validation

Align offspring's HiFi reads to parental haplotypes.
```minimap2 -y -a --eqx --cs -x map-hifi -I8g -t```

Create the coverage track on parental haplotypes with deeptools. The output is used to create coverage profile in Figure 5b as well as for other 18 q-arm recombination 

```bamCoverage -b hifi.bam -of bg -bs 10000 --minMappingQuality 30```

#### Approximate recombination breakpoint

Created a 10 kbp sliding aligner based on the all-vs-all alignments of each recombination event. The sliding aligner used the parental haplotype as reference and the 10kbp bins are created for each reference haplotype.

```rb liftover --bed 10kb_bins.bed ava_mm2.paf | rb -q stats --paf > aln.stats```

The exact 10 kbp matches are identified to define the approximate breakpoint location. The exact matches require 10 kbp segment aligned with 100% sequence identity. 
These exact matches are used to create the tracks in Figure 5c as well as for other recombination events.

The slider aligner and read depth track used in Figure5 are created by function ```figure5b()``` and ```figure5c()``` in ```transmissions/CreateFig.py```. 


#### Refine breakpoint

Find overlapping alignment around the approximate breakpoint location on child haplotype

```mashmap -s 10000 -r parent_ref.fa -q offspring_qry.fa```

Extend 30kbp on each side of the approximate breakpoint location to extrac the sequence from child and parental haplotypes. Coordinates used to extract sequences are listed in Table S2.

Create the multiple sequence alignment (MSA) for these sequences and find parental paralog sequence variants in the MSA.

#### Recombination significance

The test is done by function ```test_rc_bp()``` in ```transmission/TestRC.py```.

##### q-arm vs. euchromatin regions

Randomly shuffle 5Mbp regions on euchromatin regions and count the number of recombination events. 
The results are used to create the distribution in Figure S6b. We excluded acrocentric short arms and 5Mbp on q-arm.
```
for i in {1..5000}; do
bedtools shuffle -i qarm_5Mb.bed -g chm13.chrlen.txt -excl excl_regions.bed > tmp/shuffle_"$i".bed

bedtools intersect -c -a tmp/shuffle_"$i".bed -b ../all_autosome_bp.bed > tmp/shuffle_"$i"_bp_cnt.txt
done
```

##### acrocentric p-arm vs. metacentric p-arm

Randomly shuffle 38Mbp regions on p-arm regions and count the number of recombination events. 
The results are used to create the distribution in Figure S6c. We excluded pericentromeric regions (5Mbp within centromere) on metacentric p-arms

```
for i in {1..5000}; do
bedtools shuffle -i proximal_regions.bed -g chm13.chrlen.txt -excl excl_regions.bed > tmp/shuffle_"$i".bed
bedtools intersect -c -a tmp/shuffle_"$i".bed -b ../all_autosome_bp.bed > tmp/shuffle_"$i"_bp_cnt.txt
done
```

### De novo mutations

Variants are detected with ```snakemakes/var_detect.smk``` based on the transmitted haplotype alignments created by wfmash.

#### Mutation validation

HiFi read and ONT reads from offspring and parents are aligned to the parental reference haplotype used in variant calling.

```
## HiFi alignment
minimap2 -y -a --eqx --cs -x map-hifi -I8g -t
## ONT alignment
minimap2 -a -x map-ont --eqx -t
```

```mutation/Snakefile``` is used to validate detected variants. This pipeline will find de novo SNVs supported by HiFi and ONT reads of offspring.
The configuration file ```mutation/config.json``` requires the manifest for aligned BAM files in ```mutation/validate_var.tab```.

The IGV screenshots are created for variants validated by this pipeline (Data S3).

#### Mutation spectrum

We created the mutation spectrum for validated de novo SNVs with the pipeline ```mutation/spectrum/Snakefile```. 
The configuration file ```mutation/spectrum/config.json``` requires the manifest for aligned BAM files in ```mutation/spectrum/spectrum_var.tab```. 

The significance test is done with function ```figure6d()``` in ```mutations/MutFigs.py```

