<!--
INSTRUCTOR NOTES (not rendered on GitHub)
This exercise expects the following prepared files under /data/courses/courseB/EFNB3/:
  ref/UU_Cfam_GSD_1.0.chr5.fa (+ .fai, .dict, bwa index)       NCBI chr5 = NC_049226.1
  ref/EFNB3_region.interval_list                               e.g. NC_049226.1:27800000-37800000
  bams/WE0025.dedup.bam, WE0028.dedup.bam, WE0029.dedup.bam (+ .bai), region subsets
  bams/WE0028.sorted.bam (+ .bai), un-deduplicated, for the MarkDuplicates step
  gvcf/*.g.vcf.gz (+ .tbi)                                     fallback copies
  controls/controls.sites.vcf.gz (+ .tbi)                      sites-only, normalised, same region, from the control genomes
  snpeff/                                                      snpEff dataDir + config with database UU_Cfam_GSD_1.0_NCBI106
  fallback/*.vcf.gz                                            fallback copies of every step's output
Module versions below are the ones used in last year's course; check they still exist (module avail) and add a BCFtools module.
-->

# Variant Calling: Why Can't These Puppies Walk?

## The case

A Weimaraner breeder contacts you. In a litter of **11 puppies, 3 cannot walk normally**. From 3–4 weeks of age they move their hind legs together, hopping like rabbits instead of walking one leg after the other. Their front legs are also affected, but less. Both parents and the other 8 puppies are healthy.

Watch the video of the affected puppies (Supporting Information **Video S1** of the paper; link at the end of this page).

The three puppies were examined by a veterinary neurologist and, after euthanasia, by veterinary pathologists. The brain, spinal cord, nerves, muscles and every other organ looked **completely normal**, even under the microscope with special stains.

So the microscope found nothing. Is it genetic? You are the genomics team.

We will use the real whole-genome sequences of the three affected puppies (**WE0025, WE0028, WE0029**) from this study:

> Schwarz C, Bartenschlager F, *et al.*, Leeb T. *EFNB3 frameshift variant in Weimaraner dogs with a condition resembling a congenital mirror movement disorder.* Movement Disorders 2025;40:1883–1891. [doi:10.1002/mds.30243](https://doi.org/10.1002/mds.30243)

**Don't read the paper yet!** First find the answer yourself, as the authors had to.

#### Questions before you start (write the answers in your README)
1. Healthy parents, 3 of 11 puppies affected. Which mode of inheritance does this suggest? Which other mode can't you exclude yet, given that all three affected puppies are male?
2. If the disease is autosomal recessive, what **genotype** do you expect all three puppies to share at the causal variant? What genotype do you expect for the parents?
3. Would you expect a healthy dog of a different breed to carry the causal variant in the homozygous state?
4. Each dog genome carries ~4–5 million variants compared with the reference. How would you find the one that causes the disease?

Keep your answers. You will check them at the end.

#### The plan

```
reads ─► BAM ─► mark duplicates ─► HaplotypeCaller (one gVCF per dog)
                                          │
                       CombineGVCFs ◄─────┘
                            │
                      GenotypeGVCFs  ─► one VCF with all three dogs
                            │
                    hard filtering  ─► annotation (snpEff)
                            │
              THE FUNNEL: from thousands of variants to one
```

To keep run times short we only analyse a **10 Mb region of chromosome 5** (`NC_049226.1:27,800,000-37,800,000`). Why this region? The authors genotyped all 13 family members on a SNP chip. Linkage analysis and autozygosity mapping (looking for long stretches where all three affected puppies are homozygous for the same haplotype) left only **5 regions on 5 chromosomes, 43.2 Mb in total**, where the causal variant can be. Our region lies inside one of them.

All positions refer to the **UU_Cfam_GSD_1.0** dog reference genome (a German Shepherd, also called CanFam4). Chromosome 5 is named `NC_049226.1` in this assembly.

## Setup

```shell
cd ~/course
mkdir -p variantCalling_EFNB3/{markduplicates,gvcf,variants,funnel}
cd variantCalling_EFNB3
ln -s /data/courses/courseB/EFNB3/ref ref
ln -s /data/courses/courseB/EFNB3/bams bams
ls -l ref/ bams/
```
Questions:
1. What does `ln -s` do, and why is it better than copying the reference into your directory?
2. Which index files exist for the reference? Which program needs which index?

## Part 1: Mark duplicates

PCR can copy the same DNA fragment many times. These duplicate reads are not independent evidence. If one of them carries a PCR error, the copies make the error look like a real variant. We **mark** duplicates rather than delete them, so the variant caller can ignore them while the information stays in the file.

Two of the three dogs (WE0025, WE0029) were deposited in the European Nucleotide Archive as BAM files with duplicates already marked. You will mark duplicates for **WE0028** yourself.

```shell
cd markduplicates
code markduplicates.sh
```

```
#!/bin/bash
# Slurm options
#SBATCH --mail-user=<your.email@example.com>
#SBATCH --mail-type=fail,end
#SBATCH --job-name="markDup"
#SBATCH --chdir=.
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --partition=pcourseb
#SBATCH --output=markDup_%j.out
#SBATCH --error=markDup_%j.err

module load GATK/4.2.6.1-GCCcore-10.3.0-Java-11
module load SAMtools/1.13-GCC-10.3.0

gatk MarkDuplicates \
    -I ../bams/WE0028.sorted.bam \
    -O WE0028.dedup.bam \
    -M WE0028.dup_metrics.txt

samtools index WE0028.dedup.bam
```

```shell
sbatch markduplicates.sh
```
Questions:
1. Open `WE0028.dup_metrics.txt`. What percentage of reads are duplicates (`PERCENT_DUPLICATION`)? Is that a lot?
2. Compare `samtools flagstat` of the input and output BAM. What changed and what didn't?

## Part 2: Look before you call

Before running any variant caller, look at the reads. Open an interactive session and use `samtools tview` on one of the affected dogs:

```shell
srun --partition=pcourseb --time=1:00:00 --mem=4G --pty /bin/bash
module load SAMtools/1.13-GCC-10.3.0
samtools tview bams/WE0025.dedup.bam ref/UU_Cfam_GSD_1.0.chr5.fa
```
Press `g` and type `NC_049226.1:32799600`. Scroll right slowly with the arrow keys.

Questions:
1. Do you see a position where **all** reads disagree with the reference? Is it a substitution, an insertion or a deletion?
2. What does it mean if (almost) all reads show the change, compared with about half of them?

Don't tell your neighbours yet! Type `q` to quit and `exit` to leave the compute node.

## Part 3: One gVCF per dog (job array)

GATK's recommended workflow calls variants **per sample** in **gVCF mode** first. A gVCF records the variants and also **how confident we are that every other position is reference**. Later we combine the gVCFs and genotype all dogs together. This is how large projects (the authors used 1,489 control genomes) add new dogs without re-analysing everything.

Instead of writing three scripts, we use a **Slurm job array**: one script, run once for each dog. Slurm sets `$SLURM_ARRAY_TASK_ID` to 1, 2 or 3.

```shell
cd ../gvcf
printf "../bams/WE0025.dedup.bam\n../markduplicates/WE0028.dedup.bam\n../bams/WE0029.dedup.bam\n" > bams.list
cat bams.list
code haplotypeCaller.sh
```

```
#!/bin/bash
# Slurm options
#SBATCH --mail-user=<your.email@example.com>
#SBATCH --mail-type=fail,end
#SBATCH --job-name="hapCall"
#SBATCH --chdir=.
#SBATCH --time=3:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=pcourseb
#SBATCH --array=1-3
#SBATCH --output=hapCall_%A_%a.out
#SBATCH --error=hapCall_%A_%a.err

module load GATK/4.2.6.1-GCCcore-10.3.0-Java-11

# pick line number $SLURM_ARRAY_TASK_ID from bams.list
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bams.list)
SAMPLE=$(basename $BAM .dedup.bam)

gatk --java-options "-Xmx6g" HaplotypeCaller \
    -R ../ref/UU_Cfam_GSD_1.0.chr5.fa \
    -I $BAM \
    -L ../ref/EFNB3_region.interval_list \
    -ERC GVCF \
    -O ${SAMPLE}.g.vcf.gz
```

```shell
sbatch haplotypeCaller.sh
squeue -u $USER
```
Questions:
1. How many jobs do you see in `squeue`? What do `%A` and `%a` mean in the output file names?
2. What does `-L` do? What would happen without it?
3. Look inside a gVCF: `zcat WE0025.g.vcf.gz | grep -v "^#" | head`. What do the lines with `<NON_REF>` and an `END=` tag mean?

If the jobs take too long, copy the gVCFs: `cp /data/courses/courseB/EFNB3/gvcf/*.g.vcf.gz* .`

## Part 4: Joint genotyping

```shell
cd ../variants
code jointGenotype.sh
```

```
#!/bin/bash
# Slurm options
#SBATCH --mail-user=<your.email@example.com>
#SBATCH --mail-type=fail,end
#SBATCH --job-name="genotype"
#SBATCH --chdir=.
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --partition=pcourseb
#SBATCH --output=genotype_%j.out
#SBATCH --error=genotype_%j.err

module load GATK/4.2.6.1-GCCcore-10.3.0-Java-11

REF=../ref/UU_Cfam_GSD_1.0.chr5.fa

# merge the three per-dog gVCFs into one multi-sample gVCF
gatk CombineGVCFs -R $REF \
    -V ../gvcf/WE0025.g.vcf.gz \
    -V ../gvcf/WE0028.g.vcf.gz \
    -V ../gvcf/WE0029.g.vcf.gz \
    -O WE.combined.g.vcf.gz

# genotype all dogs together
gatk GenotypeGVCFs -R $REF \
    -V WE.combined.g.vcf.gz \
    -O WE.vcf.gz
```
Questions:
1. Which samples are in the VCF, and in which order? (`bcftools query -l WE.vcf.gz`). **Write down the order; you need it in Part 7.**
2. Why genotype the dogs together instead of using the three single-dog VCFs?

## Part 5: Hard filtering

Each variant has quality annotations: how many reads support it, whether they come from both DNA strands, how well the reads map, and so on. We flag variants that look like artefacts. SNPs and indels have different error profiles, so we filter them separately and then merge them again.

```
#!/bin/bash
# Slurm options
#SBATCH --mail-user=<your.email@example.com>
#SBATCH --mail-type=fail,end
#SBATCH --job-name="filter"
#SBATCH --chdir=.
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --partition=pcourseb
#SBATCH --output=filter_%j.out
#SBATCH --error=filter_%j.err

module load GATK/4.2.6.1-GCCcore-10.3.0-Java-11

gatk SelectVariants -V WE.vcf.gz -select-type SNP   -O WE.snps.vcf.gz
gatk SelectVariants -V WE.vcf.gz -select-type INDEL -O WE.indels.vcf.gz

gatk VariantFiltration -V WE.snps.vcf.gz -O WE.flt.snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

gatk VariantFiltration -V WE.indels.vcf.gz -O WE.flt.indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

gatk MergeVcfs -I WE.flt.snps.vcf.gz -I WE.flt.indels.vcf.gz -O WE.flt.vcf.gz
```
Save it as `filter.sh` and submit it with `sbatch filter.sh`.

Question: Look at a few lines with `zcat WE.flt.vcf.gz | grep -v "^#" | head`. Where do you see whether a variant passed?

## Part 6: What does each variant do? (snpEff)

A position number means nothing to a biologist. We want to know: **which gene, and what happens to the protein?** snpEff compares each variant to the gene annotation (NCBI annotation release 106 for UU_Cfam_GSD_1.0) and adds an `ANN` field with the predicted effect.

```
#!/bin/bash
# Slurm options
#SBATCH --mail-user=<your.email@example.com>
#SBATCH --mail-type=fail,end
#SBATCH --job-name="snpEff"
#SBATCH --chdir=.
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --partition=pcourseb
#SBATCH --output=snpEff_%j.out
#SBATCH --error=snpEff_%j.err

module load BCFtools  # check the exact module name with: module avail BCFtools

snpEff="apptainer exec -B $TMPDIR:/tmp -B /data/courses/courseB/EFNB3/snpeff /mnt/containers/apptainer/snpeff:5.2--hdfd78af_1 snpEff"

$snpEff ann -Xmx6g \
    -c /data/courses/courseB/EFNB3/snpeff/snpEff.config \
    -dataDir /data/courses/courseB/EFNB3/snpeff/data \
    -stats WE.snpEff.html \
    UU_Cfam_GSD_1.0_NCBI106 WE.flt.vcf.gz \
    | bcftools norm -f ../ref/UU_Cfam_GSD_1.0.chr5.fa -m -any -Oz -o WE.ann.vcf.gz

bcftools index -t WE.ann.vcf.gz
```
`bcftools norm` splits sites with several alternative alleles into one line each and left-aligns indels. This lets us compare our variants with the control dogs in Part 7.

Download `WE.snpEff.html` to your computer and open it in a browser.

Questions:
1. How many variants are in the region? How many are in exons?
2. Which effect types (missense, synonymous, frameshift, …) are most common? Does that make sense?

## Part 7: The funnel, from thousands of variants to one

This is the heart of the exercise. At each step, write the number of remaining variants in your README. Before each step, **predict** whether it will remove many or few variants.

```shell
cd ../funnel
ln -s ../variants/WE.ann.vcf.gz* .
module load BCFtools  # check the exact module name with: module avail BCFtools
```

**Step 1: All variants in the region**
```shell
bcftools view -H WE.ann.vcf.gz | wc -l
```

**Step 2: Only high-quality variants (PASS)**
```shell
bcftools view -f PASS -Oz -o step2.pass.vcf.gz WE.ann.vcf.gz
bcftools index -t step2.pass.vcf.gz
bcftools view -H step2.pass.vcf.gz | wc -l
```

**Step 3: Fits the recessive model**

All three affected puppies must be **homozygous for the alternative allele** (`"AA"` in bcftools). `GT[0]`, `GT[1]`, `GT[2]` are the three samples in the order you wrote down in Part 4.
```shell
bcftools view -i 'GT[0]="AA" && GT[1]="AA" && GT[2]="AA"' -Oz -o step3.recessive.vcf.gz step2.pass.vcf.gz
bcftools index -t step3.recessive.vcf.gz
bcftools view -H step3.recessive.vcf.gz | wc -l
```
Question: This step removes most variants but leaves many. Why do the three puppies share so many homozygous variants? (Hint: what is a dog breed, genetically? And what does the reference genome come from?)

**Step 4: Absent from healthy dogs**

A variant causing a rare, severe disease **cannot be common in healthy dogs**. The file `controls.sites.vcf.gz` lists every variant seen in the control dog genomes for this region. We keep only variants that are **not** in the controls (`-C` = complement):
```shell
CTRL=/data/courses/courseB/EFNB3/controls/controls.sites.vcf.gz
bcftools view -H $CTRL | wc -l
bcftools isec -C -w1 step3.recessive.vcf.gz $CTRL -Oz -o step4.private.vcf.gz
bcftools index -t step4.private.vcf.gz
bcftools view -H step4.private.vcf.gz | wc -l
```

**Step 5: Changes a protein**

snpEff classifies effects as `HIGH` (e.g. stop-gained, frameshift, splice site), `MODERATE` (e.g. missense), `LOW` (e.g. synonymous) and `MODIFIER` (intergenic, intronic, …).
```shell
bcftools view -H step4.private.vcf.gz | grep -E '\|(HIGH|MODERATE)\|' > step5.candidates.txt
wc -l step5.candidates.txt
cut -f 1-5 step5.candidates.txt
```

Fill in the funnel:

| Step | Filter | Variants left (your region) | Paper (whole genome) |
|---|---|---|---|
| 0 | SNP chip linkage ∩ autozygosity | 10 Mb region | 5 regions, 43.2 Mb |
| 1 | All variant calls | | |
| 2 | PASS | | |
| 3 | Homozygous ALT in all 3 cases | | |
| 4 | Absent from control dogs | | |
| 5 | HIGH or MODERATE effect | | **1** |

Questions:
1. Which step removed the most variants? Which step was the most *informative*?
2. Do the numbers match your predictions?
3. Step 4 used control dogs of **many different breeds**, not just Weimaraners. Why is that still useful? Could it ever remove the true causal variant?

## Part 8: Meet your candidate

Look at the annotation of your candidate(s). The `ANN` field holds one entry per transcript, separated by commas:
```shell
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\n' step4.private.vcf.gz \
  | grep -E '\|(HIGH|MODERATE)\|' \
  | tr ',' '\n' | cut -d'|' -f 1-4,7,10,11 | head
```
The fields are: allele | effect | impact | gene | transcript | cDNA change | protein change.

Questions:
1. Which gene is it? What kind of variant (SNP, insertion, deletion)? How many bases?
2. What is the predicted effect on the protein? Why does a change of this size have such a drastic effect? (Hint: how many bases is a codon?)
3. Compare with what you saw in `tview` in Part 2. Is it the same position?
4. Check the genotypes and read depths of all three dogs (replace `<POS>` with your candidate's position):
```shell
bcftools query -H -f '%POS\t%REF\t%ALT[\t%GT:%AD]\n' -r NC_049226.1:<POS> step4.private.vcf.gz
```
`AD` = allelic depth: reads supporting reference, alternative. Is the evidence convincing?

#### Look at it in IGV

Download the three BAM files for the region around the candidate (links on ILIAS) and load them in IGV together with the UU_Cfam_GSD_1.0 genome. Go to `NC_049226.1:32,799,580-32,799,680`.
1. Do all reads in all three dogs show the variant?
2. Is the region repetitive or low-complexity? Could this be an alignment artefact?
3. The gene track shows exons. In which exon is the variant?

## Part 9: Is it biologically plausible?

A variant that passes every filter is only a **candidate**. Now read about the gene:

1. **OMIA** (Online Mendelian Inheritance in Animals): search https://omia.org for your gene. Is it already known in dogs?
2. **Mouse knockouts**: search the gene at https://www.informatics.jax.org (MGI) and open "Phenotypes". What happens to mice without this gene? Compare with the puppy video.
3. **Function**: read the gene summary on https://www.ncbi.nlm.nih.gov/gene or GeneCards. Where in the body is the protein made, and when? What does a "midline barrier" in the spinal cord do?
4. **Evolution**: search NCBI Gene for this gene in chicken (*Gallus gallus*) or zebra finch. What do you find? Why might losing this gene be an **advantage** for a bird? (Think about how wings move compared with legs.)
5. **The protein**: the full protein has a signal peptide, a receptor-binding domain, a transmembrane domain and a short cytoplasmic tail. Using the protein change from Part 8, which domains does the mutant protein lose?

#### Rule out the known disease

Weimaraners already have a known recessive disease with a hopping gait: **spinal dysraphism**, caused by a frameshift in *NKX2-8*. Affected dogs also have visible malformations of the spinal cord.
1. Why is it important to exclude it?
2. The pathology was normal. Does that support or argue against spinal dysraphism?
<!-- INSTRUCTOR: once the NKX2-8 c.449delinsTT position on UU_Cfam_GSD_1.0 is confirmed, add a bcftools/tview check here and include that region in the course BAMs. -->

## Part 10: Report to the breeder

Write a short report (5–8 sentences, plain language, no jargon) to the breeder in your README. It must answer:

1. Is the disease genetic? What is the cause (gene, type of variant), in words a dog breeder understands?
2. How is it inherited? Why were the parents healthy?
3. What are the parents' genotypes, almost certainly? What fraction of their puppies would be expected to be affected if they are mated again?
4. The 8 healthy littermates: what is the chance that each of them is a carrier? (Careful: it is not 1/2!)
5. What do you recommend for breeding? Should carriers be excluded from breeding altogether?

#### For the brave: design the DNA test

The authors genotyped 134 Weimaraners with PCR and Sanger sequencing using these primers:
```
Forward: TGGACACACCTGGATTCTGA
Reverse: TCTCATAGTGGGGGCAGAAG
```
1. Find the forward primer in the reference (`samtools faidx ref/UU_Cfam_GSD_1.0.chr5.fa NC_049226.1:32799000-32800500 | grep -v ">" | tr -d '\n' | grep -ob TGGACACACCTGGATTCTGA`).
2. The reverse primer binds the other strand. Make its reverse complement with `echo TCTCATAGTGGGGGCAGAAG | rev | tr ACGT TGCA` and find it in the same way.
3. How long is the PCR product from a normal dog? From an affected dog? Could you tell them apart on an agarose gel? What about a carrier?

The authors found **3 affected, 13 carriers and 118 normal** dogs among the 134 Weimaraners. Many of them were relatives of the litter. Why does this make 13/134 an **overestimate** of the carrier frequency in the breed?

## Part 11: Check your predictions, then read the paper

Go back to your answers from the beginning. Which were right? Now read the paper (open access): https://doi.org/10.1002/mds.30243. Video S1 is in the Supporting Information.

#### Journal club: is it really a "mirror movement disorder"?

The authors named the disease after a human condition, congenital mirror movements (CMM). A group of human neurologists disagreed in a published comment:

> Trouillard O, *et al.* Comment on Schwarz *et al.* Movement Disorders 2025. https://doi.org/10.1002/mds.70001

The authors replied: https://doi.org/10.1002/mds.70000

Read the comment (it is one page) and discuss in your group:
1. What is the difference between a hopping gait and a mirror movement?
2. The commenters checked *EFNB3* in 59 human CMM patients and found nothing. Does that mean the dog finding is wrong?
3. Horses with a *DMRT3* variant can pace or tölt instead of trotting. How does that relate to the puppies?
4. Why does the name of a disease matter?

## Add your work to git

```shell
cd ~/course
code .gitignore
```
Make sure the following are ignored (add them if missing):
```
*.bam
*.bai
*.vcf.gz
*.vcf.gz.tbi
*.g.vcf.gz
*.html
```
```shell
git add variantCalling_EFNB3/*/*.sh variantCalling_EFNB3/gvcf/bams.list README.md .gitignore
git commit -m "Add EFNB3 variant calling workflow and case notes"
git push
```
Your README should now contain: your predictions, the filled funnel table, your candidate variant, and your report to the breeder.
