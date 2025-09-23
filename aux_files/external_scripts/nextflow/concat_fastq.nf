nextflow.enable.dsl=2

if (!params.concat_samplesheet) {
  exit 1, 'Concatenating input samplesheet file not found!'
}

if (!params.s3_outdir) {
  exit 1, 'Missing s3 bucket path for merged fastq files!'
}

def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id = row.sample_id

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    return fastq_meta
}

process concat_fastq {
    publishDir path: "${params.s3_outdir}", mode: 'copy'
    container "${params.use_container}"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")

    output:
    tuple val(meta), path("*.fastq.gz")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]

    if (readList.size >= 2) {
        def read1 = []
        def read2 = []
        readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
        """
        cat ${read1.join(' ')} > ${prefix}_R1_001.fastq.gz
        cat ${read2.join(' ')} > ${prefix}_R2_001.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
        END_VERSIONS
        """
    }
}

workflow {
  Channel.fromPath(params.concat_samplesheet)
    .splitCsv( header:true, sep:',' )
    .map { create_fastq_channel(it) }
    .groupTuple(by: [0]) //group by sample_id and merge fastq filepaths?
    .map {
        meta, fastq -> [ meta, fastq.flatten() ]
    } |
    concat_fastq
}