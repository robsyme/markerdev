#!/usr/bin/env nextflow

params.bams = "data/bams/*.bam"
params.reference = "data/reference/oviAri3.fa"

reference = file(params.reference)

String runClassifier(Path path) {
  match = path.getFileName().toString() =~ /Mo_Sheep_(?<key>[^_-]+)[_-]HiQ_.*.bam/
  match.matches()
  return match.group("key")
}

Channel
.fromPath(params.bams)
.map { [it, it + ".bai"] }
.set { bams }

process generateRegions {
  input:
  file reference
  
  output:
  file 'regions.txt' into regions

  """
samtools faidx ${reference}
cut -f1 ${reference}.fai > regions.txt
  """
}

regions
.splitText()
.map { it.trim() }
.spread(bams)
.groupTuple()
.set { freebayesInputs }

process freebayes {
  tag { region.trim() }
  cpus 2
  
  input:
  file reference
  set region, 'input.*.bam', 'input.*.bam.bai' from freebayesInputs

  output:
  file("out.gvcf") into gvcfs

  """
freebayes \
 --gvcf \
 --ploidy 2 \
 --region ${region} \
 --genotype-qualities \
 -f ${reference} \
 input.*.bam \
> out.gvcf
  """
}

process vcfConcat {
  publishDir "$baseDir/results/vcf", mode: 'link'

  input:
  file "split.*.vcf" from gvcfs.toList()

  output:
  set 'out.gvcf.gz', 'out.gvcf.gz.tbi' into gvcfFinal

  """
for vcf in *.vcf; do
  bgzip \$vcf
  tabix -p vcf \$vcf.gz
done

vcf-concat *.vcf.gz | vcf-sort | bgzip > out.gvcf.gz
tabix out.gvcf.gz
  """
}

process filterSites {
  publishDir "$baseDir/results/filtering", mode: 'link'
  cpus 4

  input:
  set 'in.gvcf.gz', 'in.gvcf.gz.tbi' from gvcfFinal

  output:
  file 'goodregions.bed' into goodregions

  """
zcat in.gvcf \
| bio-vcf --filter 'r.alt.first == "<*>"' --eval '[r.chrom, r.pos-1, r.info.end]' --seval 's.gq > 200' --num-threads ${task.cpus} \
| grep -v false \
| cut -f 1-3 \
| sort -k1,1 -k2,2n \
> goodregions.bed
  """
}
