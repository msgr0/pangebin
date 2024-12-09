#!/usr/bin/env nextflow

process NCBI_RETRIVE {
  beforeScript 'sleep 2'
  cache 'lenient'
  maxForks 1
  errorStrategy 'ignore'

  input:
  val(meta)

  output:
  tuple val(meta), path(fna_gz), emit: genome

  script:
  fna_gz = "${meta.id}.gen.fna.gz"
  """
  link=\$(python $projectDir/ncbi_link.py --input ${meta.id})
    rsync -vPp --chmod=u+w \$link ${fna_gz}
  """
}


process BLAST {
  conda 'bioconda::htslib blast'
  errorStrategy 'ignore'
  cache 'lenient'

  input:
  tuple val(meta), path(fasta), path(fna_gz)

  output: 
  tuple val(meta), path(blast_out), emit: blast_out
  tuple val(meta), path(blast_gt), emit: blast_gt

  script:
  fna = "${meta.id}.fna"
  blast_out = 
  blast_gt = 
  """
  bgzip -d -c ${genomic_fna_gz} > ${fna}

  makeblastdb -in ${fna} -dbtype nucl
  blastn -db ${fna} -query ${fasta} -out ${blast_out} -outfmt "6 qseqid sseqid pident length mismatch qstart qend bitscore"
  # BLAST TO GT, GRAPH, REFERENCE (if graph == reference, mode assembly, if graph != reference, mode pangenome, blast scores tsv)
  """
}

process BUILD_GT {
  cache 'lenient'

  input:
  tuple val(meta), path(blast_out), path(graph), path(reference)

  output:
  tuple val(meta), path(ground_truth), emit: gt

  script:
  gt_file = 
  """
  python3 $projectDir/bin/blast_utils.py --assembly ${fna} --csv ${blast_tsv} --fasta ${fasta_AB} --output ${ground_truth_csv}
  """
}

process STATS {
  input:
  tuple val(meta), path(graph), path(ground_truth)

  output:
  tuple val(meta), path(stats)

  script:
  stats = "${meta.id}.stats.txt"
  """
  python3 $projectDir/graph_stats.py --sample ${meta.id} --input ${graph} --groundt ${ground_truth} --output ${stats}
  """
}


process EVAL {
  input:
  tuple val(meta), path(prediction), path(ground_truth)
  val type

  output:
  tuple val(meta), path("*.pdf"), path("*.mix.csv"), emit: mix, optional: true

  script:
  """
  python3 $projectDir/evaluate_bins.py --sample ${meta.id} --gfa ${graph} --csv ${gt_csv} --bin ${bin_tsv} --type ${type} --output .
  
  """
}

// process COMPARE {
//   input:
//   tuple val(meta), path(scores_pangenome), path(scores_assembly)
//   // file ==
//   // prec:
//   // rec:
//   // f1:

//   output:
//   tuple val(meta), path(scores), path(graph), emit: results

//   script:
//   scores = 
//   // produce a comparison chart (maybe gnuplot for fast plotting???)
//   """
  
//   """

}

workflow EVALUATE {
  take:
  sample // meta, pangenome, assembly, prediction_pan, prediction_assembly
  
  main:
  pangenome = sample.map{meta, graph, graph_fasta, ref, ref_fasta, pred, ref_pred -> tuple(meta, graph, graph_fasta, pred)

  reference = sample.map{meta, graph, graph_fasta, ref, ref_fasta, pred, ref_pred -> tuple(meta, ref, ref_fasta, ref_pred)

  id = sample.map{meta, graph, graph_fasta, ref, ref_fasta, pred, ref_pred -> meta}

  NCBI_RETRIVE(id)
  EVALUATE_PANGENOME(pangenome, reference, NCBI_RETRIVE.out.genome)
  EVALUATE_ASSEMBLY(reference, NCBI_RETRIVE.out.genome)

  COMPARE(EVALUATE_PANGENOME.out.score.join(EVALUATE_ASSEMBLY.out.score))

  emit:
  pangenome_score = EVALUATE_PANGENOME.out.score
  assembly_score = EVALUATE_ASSEMBLY.out.score 
  comparative = COMPARE.out.results

}



workflow EVALUATE_ASSEMBLY {
  take:
  graph
  genome

}
workflow EVALUATE_PANGENOME {
  take:
  pangenome
  assembly
  genome

  main:
  BLAST(genome.join(pangenome.map{it, graph, fasta -> tuple(it, fasta)}))
  BUILD_GT(BLAST.out.blast_gt.join(
    assembly.map{it, graph, fasta -> tuple(it, fasta)}).join(
      pangenome.map{it, graph, fasta -> tuple(it, graph)}
    )
  )

  emit:
  ground_truth
  scores
}