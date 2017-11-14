#!/usr/bin/env nextflow

//params.setup = 'setup.tab'
params.in = '/lustre/scratch/users/elin.axelsson/first_set/*_{1,2}.fastq'
params.fasta = "/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
params.dna_fasta     = "/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
params.gtf 	     = "/lustre/scratch/projects/berger_common/backup_berger_common/gtf/Arabidopsis_thaliana.TAIR10.35.gtf"
params.fragment_len  = '180'
params.fragment_sd   = '20'
params.bootstrap     = '100'
params.output        = "results/"
params.normtosize = '119146348'


fasta=file(params.fasta)
fasta_dna=file(params.dna_fasta)
gtf=file(params.gtf)


// paired or single, but single have to be names _1.fq

Channel
    .fromFilePairs( params.in, size: -1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.in}" }
    .set { read_files }
  
// create indeces if needed

process STARindex {
storeDir '/lustre/scratch/projects/berger_common/backup_berger_common/'

    input:
    file fasta_dna
    file gtf

    output: 
    file 'star' into star_index

    script:
    """
    STAR --runThreadN 4 --runMode genomeGenerate --genomeDir star --genomeFastaFiles ${fasta_dna} --sjdbGTFfile ${gtf} 
    """
}



process kallistoIndex {
storeDir '/lustre/scratch/projects/berger_common/backup_berger_common'
    input:
    file fasta

    output:
    file "tair10_transcripts.idx" into transcriptome_index

    script:
    """
    kallisto index -i tair10_transcripts.idx ${fasta} 
    """
}

// start 

process trimgalore {
 	publishDir "$params.output/$name/trimgalore", mode: 'copy' , pattern: '*.{txt,html,zip}'
	tag "data: $name"

	input:
	set val(name), file(reads) from read_files 

	output:
	set name, file("*.fq") into trimmed_reads
	file "*trimming_report.txt"
	file "*_fastqc.html"
	file "*_fastqc.zip"

	script:
	def single = reads instanceof Path
    	if( single ) {
	"""
	trim_galore --fastqc ${reads}
	"""
	}
	else {
	"""
	trim_galore --fastqc --paired ${reads}
	"""
	}
}

// need this channel several times

trimmed_reads.into {trimmed_reads_star; trimmed_reads_kallisto}


process star {
	publishDir "$params.output/$name/star_${name}", mode: 'copy'
	tag "star: $name"

    	input:
    	file index from star_index
    	set name, file(fq) from trimmed_reads_star
    
   	output:
    	set name, file("Aligned.out.sam") into samfiles
	file "*.out"
	file "*.tab"   

    	script:
    	"""
    	STAR --genomeDir $index --readFilesIn $fq --runThreadN 4 --quantMode GeneCounts  
    	"""
}

process sam2bam {
	publishDir "$params.output/$name/bam_bw", mode: 'copy'
	tag "bam: $name"

	input:
	set name, file(sam) from samfiles
	

	output:
	set name, file("Aligned_sorted_${name}.bam") into bamfiles
	
	script:
	"""
	samtools view -b ${sam} | samtools sort -o Aligned_sorted_${name}.bam - 
	"""
}

process bam2bw {
	publishDir "$params.output/$name/bam_bw", mode: 'copy'
        tag "bw: $name"

	input:
	set name, file(bam) from bamfiles


	output:
	file("${name}.bw")
	
	script:
	"""
	export TMPDIR=\$(pwd)
   	samtools index ${bam}
   	bamCoverage -b ${bam} -o ${name}.bw --normalizeTo1x ${params.normtosize} --binSize=10	
	"""
}






process quantKallisto {
publishDir "$params.output/$name", mode: 'copy'
tag "fq: $name"

    input:
    file index from transcriptome_index
    set name, file(fq) from trimmed_reads_kallisto

    output:
    file "kallisto_${name}" into kallisto_dirs 

    script:
    def single = fq instanceof Path
    if( single ) {
    """
    kallisto quant -i ${index} -o kallisto_${name} --single -l ${params.fragment_len} -s ${params.fragment_sd} -b ${params.bootstrap} ${fq}
    """ 
    }
    else {
    """
    kallisto quant -i ${index} -o kallisto_${name} -b ${params.bootstrap} ${fq}
    """
    }
}


workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}










	
	




