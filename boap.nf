#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* ------------------------------------------------------------------------
 * ONT WGS ASSEMBLY PIPELINE
 * ------------------------------------------------------------------------ */

def version = "1.0.0"

// Display version and exit if --version is called
if (workflow.commandLine.contains('--version')) {
    println "BAOP Version: $version"
    exit 0
}

// --- 1. CONFIGURATION ---
params.input          = "input_fastq"
params.reads          = "${params.input}/*.fastq.gz" 
params.threads        = 30
params.outdir         = "baop_output"
params.coverage       = 100
params.min_contig_len = 1000
params.gsize          = null
params.force_model    = false 

// --- 2. GSIZE PARSING ---
def manual_gsize_val = null
if (params.gsize) {
    def input_str = params.gsize.toString().toLowerCase()
    if (input_str.endsWith('mb')) {
        manual_gsize_val = (input_str.minus('mb') as BigDecimal) * 1000000
    } else if (input_str.endsWith('m')) {
        manual_gsize_val = (input_str.minus('m') as BigDecimal) * 1000000
    } else {
        manual_gsize_val = input_str as BigDecimal
    }
    manual_gsize_val = manual_gsize_val.toLong()
}


// --- 3. Processes ---

process BAM2FASTQ {
    tag "$sampleid"
    cpus 4

    input:  tuple val(sampleid), path(bam_file)
    output: tuple val(sampleid), path("${sampleid}.fastq.gz"), emit: fastq

    script:
    """
    samtools fastq -T '*' -@ ${task.cpus} ${bam_file} | gzip > ${sampleid}.fastq.gz
    """
}

process NANOQ {
    tag "$sampleid"
    cpus 1
    
    input:  tuple val(sampleid), path(reads)
    output: tuple val(sampleid), path("${sampleid}_q10.fastq.gz"), emit: filtered_reads 

    script:
    """
    nanoq --input ${reads} --min-len ${params.min_contig_len} --min-qual 10 --output-type g -o ${sampleid}_q10.fastq.gz
    """
}

process GENOME_SIZE {
    tag "$sampleid"
    cpus 4

    input:  tuple val(sampleid), path(reads)
    output: tuple val(sampleid), path(reads), path("gsize.txt"), emit: info

    script:
    """
    raven ${reads} --disable-checkpoints -p 0 -t ${task.cpus} | awk '!/^>/ {len += length(\$0)} END {print len}' > gsize.txt
    """
}

process GSIZE {
    tag "$sampleid"
    executor 'local' 
    
    input: 
        tuple val(sampleid), path(reads)
        val(manual_size)
    output: 
        tuple val(sampleid), path(reads), path("gsize.txt"), emit: info

    script:
    """
    echo "${manual_size}" > gsize.txt
    """
}

process FILTLONG {
    tag "$sampleid"
    cpus 1

    input:  tuple val(sampleid), path(reads), path(gsize_file)
    output: tuple val(sampleid), path("${sampleid}_${params.coverage}x.fastq.gz"), path(gsize_file), emit: reads

    script:
    """
    TARGET_BASES=\$(( \$(cat ${gsize_file}) * ${params.coverage} ))
    filtlong --target_bases \$TARGET_BASES --length_weight 0 --window_q_weight 0 ${reads} | gzip > ${sampleid}_${params.coverage}x.fastq.gz
    """
}

process FLYE {
    tag "$sampleid"
    cpus 8
    
    publishDir "${params.outdir}/flye/${sampleid}", mode: "copy", pattern: "*.{txt,log,gfa,gv}"

    input:  tuple val(sampleid), path(reads), path(gsize_file)
    output: tuple val(sampleid), path("${sampleid}_flye.fasta"), path(reads), emit: assembly
            path "assembly_info.txt", emit: info
            path "flye.log"
            path "assembly_graph*"

    script:
    """
    flye --nano-hq ${reads} --out-dir . --threads ${task.cpus} --genome-size \$(cat ${gsize_file}) --iterations 1

    # Add info (cov, circular)
    awk 'NR>1 {printf "%s\\tcov=%s circular=%s\\n", \$1, \$3, (\$4=="Y"?"Y":"N")}' assembly_info.txt > contig_info.txt
    
    # Duplicate contig headers and merge info  
    awk '/^>/ {print \$0 " " substr(\$0,2)} !/^>/ {print \$0}' assembly.fasta \
        | seqkit replace -p " (.+)\$" -r " {kv}" -k contig_info.txt \
        | seqkit sort --reverse --by-length \
        | seqkit replace -p ^ -r "contig_{nr} orig_" --nr-width 5 \
        | seqkit seq -m ${params.min_contig_len} -o ${sampleid}_flye.fasta
    """
}

process DNAAPLER {
    tag "$sampleid"
    
    publishDir "${params.outdir}/flye/${sampleid}", mode: "copy", pattern: "*_summary.tsv"

    input:  tuple val(sampleid), path(assembly), path(reads)
    output: tuple val(sampleid), path("${sampleid}_reo.fasta"), path(reads), emit: assembly_reo
            path "dnaapler/contigs_all_reorientation_summary.tsv"

    script:
    """
    grep ">" ${assembly} | grep "circular=N" | sed 's/^>//g' > noncircular_contigs.txt || true
    
    dnaapler all --input ${assembly} --output dnaapler --threads ${task.cpus} --prefix contigs --ignore noncircular_contigs.txt
    
    mv dnaapler/contigs_reoriented.fasta ${sampleid}_reo.fasta
    """
}

process MEDAKA2_POLISH {
    tag "$sampleid"

    cpus 4

    // define the model flag conditionally based on params.force_model
    def model_flag = params.force_model ? "-m ${params.force_model}" : ""

    publishDir "${params.outdir}/assembly", mode: "copy", pattern: "*.fasta"

    input:  tuple val(sampleid), path(assembly), path(reads)
    output: tuple val(sampleid), path("${sampleid}_polished.fastq"), path("${sampleid}.fasta"), emit: polished

    script:
    
    """
    medaka_consensus -i ${reads} -d ${assembly} -o . -t ${task.cpus} ${model_flag} -q
    
    grep ">" ${assembly} | sed 's/^>//g' | sed 's/orig_[^ ]* //g' > original_headers_clean.txt

    POLISHED_FILE=\$(ls *.fastq | grep -v "${sampleid}" | head -n1)
    seqkit sort -n \$POLISHED_FILE > sorted_consensus.fastq

    seqkit fx2tab -l -g -n -i -H sorted_consensus.fastq > new_lengths.txt
    
    paste -d' ' <(awk '{print \$1}' original_headers_clean.txt) \
                <(awk 'NR>1{print "len="\$2}' new_lengths.txt) \
                <(awk '{\$1=""; print}' original_headers_clean.txt | sed 's/^ *//' | sed 's/  */ /g') \
                | sed 's/ /\\t/' > merged_headers.txt

    seqkit replace -p \\\$ -r " contig_{nr}" --nr-width 5 sorted_consensus.fastq \
        | seqkit replace -p ' (.+)\$' -r ' {kv}' -k merged_headers.txt --keep-key \
        > ${sampleid}_polished.fastq

    seqkit fq2fa ${sampleid}_polished.fastq -o ${sampleid}.fasta
    """
}

process ALPAQA {
    tag "$sampleid"
    
    input:  tuple val(sampleid), path(fastq), path(fasta)
    output: path "${sampleid}_alpaqa.tsv", emit: stats
            path "${sampleid}_contig_info.txt", emit: info

    script:
    """
    alpaqa.py -i ${fastq} -o ${sampleid}_alpaqa.tsv -t 1

    awk -v sid="${sampleid}" '
    BEGIN { 
        OFS="\\t"; 
        print "Strain","Contig","Length","Cov","Circ","Rot","Rot_gene"; 
    }
    /^>/ {
        sub(/^>/, "", \$1);
        contig_id = \$1;
        
        len = cov = circular = rotated = rotated_gene = "NA";

        for (i = 2; i <= NF; i++) {
            # Note: orig_ tags are stripped in previous step, so we dont parse them here
            if (\$i ~ /^len=/) len = substr(\$i, 5);
            else if (\$i ~ /^cov=/) cov = substr(\$i, 5);
            else if (\$i ~ /^circular=/) circular = \$i;       
            else if (\$i ~ /^rotated=/) rotated = \$i;         
            else if (\$i ~ /^rotated_gene=/) rotated_gene = \$i;
        }

        print sid, contig_id, len, cov, circular, rotated, rotated_gene;
    }' ${fasta} > ${sampleid}_contig_info.txt
    """
}

// --- 4. WORFLOW ---

workflow {
    def input_path = file(params.input)
    
    if (input_path.isFile()) {
        raw_ch = Channel.fromPath(input_path)
                        .map { file -> tuple(file.simpleName, file) }
    } else {
        raw_ch = Channel.fromPath( "${params.input}/*.{fastq.gz,bam}" )
                        .map { file -> tuple(file.simpleName, file) }
    }

    raw_ch.branch {
        bam:   it[1].name.endsWith('.bam')
        fastq: true
    }.set { input_sources }

    BAM2FASTQ(input_sources.bam)
    reads_ch = input_sources.fastq.mix(BAM2FASTQ.out.fastq)

    // --- PIPELINE START ---
    
    NANOQ(reads_ch)
    
    if (manual_gsize_val) {
        GSIZE(NANOQ.out.filtered_reads, manual_gsize_val)
        filtlong_input = GSIZE.out.info
    } else {
        GENOME_SIZE(NANOQ.out.filtered_reads)
        filtlong_input = GENOME_SIZE.out.info
    }

    FILTLONG(filtlong_input)
    FLYE(FILTLONG.out.reads)
    DNAAPLER(FLYE.out.assembly)
    MEDAKA2_POLISH(DNAAPLER.out.assembly_reo)
    ALPAQA(MEDAKA2_POLISH.out.polished)
    ALPAQA.out.stats.collectFile(name: 'alpaqa_report.tsv', storeDir: "${params.outdir}/summary", keepHeader: true, skip: 1)
    ALPAQA.out.info.collectFile(name: 'contig_report.tsv', storeDir: "${params.outdir}/summary", keepHeader: true, skip: 1)
}