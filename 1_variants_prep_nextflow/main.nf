nextflow.enable.dsl=2

import groovy.yaml.YamlSlurper

/*
helper utilities
*/

def stripExtensions(String name) {
    name.replaceAll(/(\.(txt|tsv|csv|gz))+$/, '')
}

def derivePrefix(String genelistPath) {
    def base = new File(genelistPath).getName()
    def prefix = stripExtensions(base)
    if (!prefix) {
        error "could not derive prefix from genelist name: ${base}"
    }
    return prefix
}

def deriveOutdir(String genelistPath, String prefix, String override) {
    if (override) {
        return override
    }

    def envDs = System.getenv("DATASET_ID")
    def envSample = System.getenv("SAMPLE_ID")
    if (envDs && envSample) {
        def root = System.getenv("ROOT_DIR") ?: System.getenv("PDC_TMP") ?: "/cfs/klemming/scratch/m/mmarandi"
        return "${root}/experiments/${envDs}/${envSample}/01_inputs/"
    }

    def marker = "/data/initial/"
    if (genelistPath.contains(marker)) {
        def parts = genelistPath.split(marker, 2)
        return "${parts[0]}${marker.replace('/initial/', '/intermediate/')}${prefix}/"
    }
    return "${projectDir}/results/${prefix}/"
}

def readRuntimeConf(String path) {
    if (!path) {
        return [:]
    }
    def f = file(path)
    if (!f.exists()) {
        error "runtime_conf not found: ${path}"
    }
    new YamlSlurper().parse(f)
}

def padLabel(Integer pad) {
    if (pad % 1000 == 0) {
        return "${(int) (pad / 1000)}kb"
    }
    return "${pad}bp"
}

/*
parameters and defaults
*/

def runtime = readRuntimeConf(params.runtime_conf)

if (!params.genelist) {
    error "--genelist is required"
}
if (!params.mane) {
    error "--mane is required"
}

def genelistPath = file(params.genelist)
def manePath = file(params.mane)

if (!genelistPath.exists()) {
    error "genelist not found: ${genelistPath}"
}
if (!manePath.exists()) {
    error "mane not found: ${manePath}"
}

def padTss = (params.pad_tss ?: runtime.pad_tss ?: 10_000) as Integer
def prefix = params.prefix ?: derivePrefix(genelistPath.toString())
def outdir = deriveOutdir(genelistPath.toString(), prefix, params.outdir ?: runtime.outdir)
def gnomadHt = params.gnomad_ht ?:
        runtime.gnomad_ht ?:
        "gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht/"
def sparkConf = params.spark_conf ?: runtime.spark_conf
def gcsConnectorJar = params.gcs_connector_jar ?: runtime.gcs_connector_jar
def hailHome = params.hail_home ?: runtime.hail_home
def tmpDir = params.tmp_dir ?: runtime.tmp_dir ?: "./tmp"
def sparkConfFile = file(sparkConf ?: "${projectDir}/conf/spark_local.json")
if (!sparkConfFile.exists()) {
    error "spark_conf not found: ${sparkConfFile}"
}

def bedLabel = "tss±${padLabel(padTss)}"
def bedFilename = "${prefix}_${bedLabel}.bed"
def qcFilename = "${prefix}_mapping_qc.tsv"
def missingFilename = "${prefix}_missing_ensg.txt"
def mergedFilename = "${prefix}_variants.tsv.gz"

log.info "prefix: ${prefix}"
log.info "outdir: ${outdir}"
log.info "pad_tss: ${padTss}"
log.info "bed filename: ${bedFilename}"
log.info "gnomad_ht: ${gnomadHt}"
log.info "spark_conf: ${sparkConfFile}"
log.info "tmp_dir: ${tmpDir}"

/*
channels
*/

Channel
    .of(genelistPath)
    .set { genelist_ch }

Channel
    .of(manePath)
    .set { mane_ch }

Channel
    .value(prefix)
    .set { prefix_ch }

Channel
    .value(outdir)
    .set { outdir_ch }

Channel
    .value(padTss)
    .set { pad_ch }

Channel
    .value(gnomadHt)
    .set { gnomad_ht_ch }

Channel
    .value(tmpDir)
    .set { tmp_dir_ch }

Channel
    .value(gcsConnectorJar)
    .set { gcs_connector_jar_ch }

Channel
    .value(hailHome)
    .set { hail_home_ch }

Channel
    .fromPath(sparkConfFile)
    .set { spark_conf_ch }

/*
processes
*/

process EXTRACT_ENSG {
    tag "$prefix"
    publishDir { outdir }, mode: 'copy', overwrite: true, pattern: 'ensg_list.txt'

    input:
    path genelist
    val prefix
    val outdir

    output:
    path "ensg_list.txt", emit: ensg_list

    script:
    """
    extract_ensg.py --genelist ${genelist} --output ensg_list.txt
    """
}

process BUILD_TSS_BED {
    tag "$prefix"
    publishDir { outdir }, mode: 'copy', overwrite: true

    input:
    path ensg_list
    path mane
    val pad_tss
    val prefix
    val outdir
    val bed_filename
    val qc_filename
    val missing_filename

    output:
    path "${bed_filename}", emit: bed
    path "${qc_filename}", emit: qc
    path "${missing_filename}", emit: missing

    script:
    """
    build_tss_bed.py \
      --ensg-list ${ensg_list} \
      --mane ${mane} \
      --pad-tss ${pad_tss} \
      --bed-out ${bed_filename} \
      --qc-out ${qc_filename} \
      --missing-out ${missing_filename}
    """
}

process HAIL_EXPORT_PER_GENE {
    tag "$prefix"
    publishDir { outdir }, mode: 'copy', overwrite: true, pattern: "per_gene/*"

    input:
    path bed
    path spark_conf
    val gnomad_ht
    val tmp_dir
    val gcs_connector_jar
    val hail_home
    val prefix
    val outdir

    output:
    path "per_gene", emit: per_gene_dir

    script:
    """
    spark_arg=""
    if [[ -n "${spark_conf}" && "${spark_conf}" != "null" ]]; then
      spark_arg="--spark-conf ${spark_conf}"
    fi
    gcs_arg=""
    if [[ -n "${gcs_connector_jar}" && "${gcs_connector_jar}" != "null" ]]; then
      gcs_arg="--gcs-connector-jar ${gcs_connector_jar}"
    fi
    hail_arg=""
    if [[ -n "${hail_home}" && "${hail_home}" != "null" ]]; then
      hail_arg="--hail-home ${hail_home}"
    fi

    hail_export_per_gene.py \
      --bed ${bed} \
      --gnomad-ht ${gnomad_ht} \
      --outdir per_gene \
      --tmp-dir ${tmp_dir} \
      ${spark_arg} \
      ${gcs_arg} \
      ${hail_arg}
    """
}

process MERGE_VARIANTS {
    tag "$prefix"
    publishDir { outdir }, mode: 'copy', overwrite: true

    input:
    path per_gene_dir
    val prefix
    val outdir
    val merged_filename

    output:
    path "${merged_filename}", emit: merged

    script:
    """
    merge_variants.py \
      --per-gene-dir ${per_gene_dir} \
      --output ${merged_filename}
    """
}

/*
workflow
*/

workflow {
    ensg = EXTRACT_ENSG(genelist_ch, prefix_ch, outdir_ch)
    bed = BUILD_TSS_BED(
        ensg.ensg_list,
        mane_ch,
        pad_ch,
        prefix_ch,
        outdir_ch,
        bedFilename,
        qcFilename,
        missingFilename
    )
    per_gene = HAIL_EXPORT_PER_GENE(
        bed.bed,
        spark_conf_ch,
        gnomad_ht_ch,
        tmp_dir_ch,
        gcs_connector_jar_ch,
        hail_home_ch,
        prefix_ch,
        outdir_ch
    )
    MERGE_VARIANTS(per_gene.per_gene_dir, prefix_ch, outdir_ch, mergedFilename)
}
