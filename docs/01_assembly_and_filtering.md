# 01. Assembly and filtering

Processing of 3RADseq data for the bottlenose wedgefish (*Rhynchobatus australiae*, RA), sequenced on one NovaSeq S4 run (paired-end, 2 x 150 bp) together with a second wedgefish species.

Steps marked **[HPC]** were run on a cluster and **[R]** in R. Software: FastQC 0.11.5, Cutadapt 1.9.1, Stacks 2.59, VCFtools 0.1.15.

## Final dataset

| SNPs | n | Missingness | Filters |
|---|---|---|---|
| 27 562 | 71 | 6.32% | R80Het50MAC3, minGQ20minDP5miss0.8, singleSNP |

## Settings

```bash
SPECIES=RA
WORKDIR=$HOME/radseq                         # edit
READS=$WORKDIR/${SPECIES}_demultiplex/reads
CORE=$WORKDIR/${SPECIES}_assembly/core
ADAPTER_1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
ADAPTER_2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```

## Step 1. Raw data and quality control [HPC]

File integrity was checked with `md5sum -c md5sum.txt`. Each plate was sequenced with two dual-indexing combinations, which were merged before QC.

```bash
cat ${PLATE}-1_*_L001_R1_001.fastq.gz ${PLATE}-2_*_L001_R1_001.fastq.gz > ${PLATE}_L001_R1_001.fastq.gz
cat ${PLATE}-1_*_L001_R2_001.fastq.gz ${PLATE}-2_*_L001_R2_001.fastq.gz > ${PLATE}_L001_R2_001.fastq.gz

fastqc ${PLATE}_L001_R1_001.fastq.gz ${PLATE}_L001_R2_001.fastq.gz
```

## Step 2. Demultiplexing [HPC]

Reads were demultiplexed with `process_radtags` (double-inline barcodes; MspI and BamHI). Untrimmed reads, Cutadapt-trimmed reads and untrimmed reads with the adapters passed to `process_radtags` were compared in R by reads retained per sample, and the last was selected.

```bash
process_radtags \
  -1 ${PLATE}_L001_R1_001.fastq.gz -2 ${PLATE}_L001_R2_001.fastq.gz \
  -i gzfastq -b barcodes.txt --inline_inline \
  -o $READS \
  --adapter_1 $ADAPTER_1 --adapter_2 $ADAPTER_2 --adapter_mm 2 \
  -c -q -r -t 140 -w 0.15 -s 10 \
  --renz_1 mspI --renz_2 bamHI
```

## Step 3. Parameter optimisation (RADstackshelpR) [HPC, then R]

Twenty samples were used (about equal numbers per locality, close to average read depth), with remainder reads (`*.rem.*.fq.gz`) moved out of the directory. Each parameter was varied in turn (m = 3 to 7, M = 1 to 8, n = M-1, M and M+1) and the results visualised in R. See [`scripts/01_mMn_optimisation.R`](../scripts/01_mMn_optimisation.R) (RADstackshelpR); the optimal values were **m = 3, M = 3, n = 4**.

```bash
OPTIM=$WORKDIR/${SPECIES}_assembly/optim

for m in {3..7}; do                                  # m
  mkdir -p $OPTIM/denovo.m${m}
  denovo_map.pl --samples $READS --popmap $OPTIM/popmap_optim.txt \
    --out-path $OPTIM/denovo.m${m} --paired -m $m -T 12 \
    -X "populations: --vcf" > $OPTIM/denovo.m${m}/denovo_map.log 2>&1
done

for M in {1..8}; do                                  # M, with m = 3
  mkdir -p $OPTIM/denovo.bigM${M}
  denovo_map.pl --samples $READS --popmap $OPTIM/popmap_optim.txt \
    --out-path $OPTIM/denovo.bigM${M} --paired -m 3 -M $M -T 12 \
    -X "populations: --vcf" > $OPTIM/denovo.bigM${M}/denovo_map.log 2>&1
done

for n in {2..4}; do                                  # n = M-1, M, M+1, with m = 3, M = 3
  mkdir -p $OPTIM/denovo.n${n}
  denovo_map.pl --samples $READS --popmap $OPTIM/popmap_optim.txt \
    --out-path $OPTIM/denovo.n${n} --paired -m 3 -M 3 -n $n -T 12 \
    -X "populations: --vcf" > $OPTIM/denovo.n${n}/denovo_map.log 2>&1
done
```

## Step 4. De novo assembly [HPC, then R]

```bash
OUT=$CORE/denovo_m3_M3_n4; mkdir -p $OUT

denovo_map.pl --samples $READS --popmap $CORE/popmap_full.txt \
  --out-path $OUT --paired -m 3 -M 3 -n 4 -T 12 \
  -X "populations: --vcf" > $OUT/denovo_map.log 2>&1
```

Assembly checks (problematic samples were flagged in R):

```bash
stacks-dist-extract denovo_map.log cov_per_sample | grep -v '^#'          # coverage after ustacks
grep -B 2 -A 3 '^Genotyped' gstacks.log                                   # loci and coverage after gstacks
stacks-dist-extract gstacks.log.distribs effective_coverages_per_sample | grep -v '^#' | cut -f 1-2,4-5,8
```

## Step 5. Filtering

### 5.1 Base filter with `populations` [HPC]

All individuals were treated as one population. Loci had to be present in at least 80% of samples, with maximum observed heterozygosity 0.5 and a minor allele count of at least 3 (R80/Het50/MAC3). Additional base combinations tested (R60_Het50_MAC3/ R80_Het50_MAC2).

```bash
populations --in-path $OUT --out-path $CORE/filtering/base_R80_Het50_MAC3/ \
  --popmap popmap_onepop.txt \
  --min-samples-overall 0.80 --max-obs-het 0.5 --min-mac 3 --vcf
```

### 5.2 Genotype filters with VCFtools [HPC, then R]

Four filter sets were compared in R, and the one with the best balance of missingness and SNP retention (FS1) was used.

| Set | min GQ | min DP | max missing |
|---|---|---|---|
| FS1 | 20 | 5 | 0.8 |
| FS2 | 20 | 8 | 0.8 |
| FS3 | 20 | 8 | 0.9 |
| FS4 | 25 | 10 | 0.9 |

```bash
# FS1 shown; the other sets change only these three values
vcftools --vcf populations.snps.vcf --minGQ 20 --minDP 5 --max-missing 0.8 \
  --recode --recode-INFO-all --out FS1_GQ20_DP5_MISS0.8
```

### 5.3 Remove replicates, high-missingness individuals and outliers [HPC, then R]

Replicates (keeping the replicate with the lowest missingness), individuals with more than 40% missing data and outliers identified in a preliminary PCA in R were removed.

```bash
vcftools --vcf FS1_GQ20_DP5_MISS0.8.recode.vcf --remove remove_list.txt \
  --recode --recode-INFO-all --out FS1_filtered
```

### 5.4 Final dataset [HPC]

The Stacks locus IDs of the filtered SNPs were extracted from the CHROM column of the filtered VCF as a whitelist, and `populations` was rerun on the original assembly to keep one SNP per locus and write the VCF and other formats. The whitelist must contain only the locus IDs, one per line.

```bash
grep -v '^#' FS1_filtered.recode.vcf | cut -f1 | sort -n | uniq > whitelist.txt

populations --in-path $OUT --out-path $CORE/final/popgen1 \
  --popmap popmap_final.txt --whitelist whitelist.txt \
  --write-single-snp --vcf --structure --radpainter
```

Final QC:

```bash
vcftools --vcf populations.snps.vcf --missing-indv
awk 'NR>1 {sum+=$5} END {print sum/(NR-1)}' out.imiss     # mean individual missingness
grep -vc '^#' populations.snps.vcf                          # number of SNPs
```
