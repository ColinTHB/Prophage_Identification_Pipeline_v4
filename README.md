Code for pipeline for identification and extraction of prophage elements in genomic assemblies. It depends on several conda environments with particular tools loaded to function.

The pipeline identifies prophage sequences with virsorter2 and then annotates "prophage" sequences picked up by virsorter2 with pharokka.
These prophage sequences are clipped to their respective "approximate" genome size based on the annotation of hallmark genes.
Outputs are prophage sequences as gbk files and tsv with virsorter2 hit details and predicted prophage boundaries on the host genome.

The current version of the pipeline depends on the presence of conda environments with the installed software.

##################################################################################################################################################################################

These are as follows.
  
  virsorter2_env - virsorter2 (v2.2.4)
 
  pullseq_env - pullseq (v1.0.2)
 
  pharokka_env - pharokka (v1.7.2)
  
  bedtools_env - bedtools (v2.31.1) & biopython (v1.84)

##################################################################################################################################################################################

Usage is as follows

  1. rewrite headers of input fast to ensure compatibility with the script

./rewrite_fasta_headers-4-prophage-hunting-v1.1.sh genome.fna input-genome-fixed.fna

  2. Run the prophage hunting script

./prophage-hunting-v4.1.sh -i input-genome-fixed.fna 

  3. Extract prophage sequences from pharokka annotation

./prophage_extract_gbk+tsv_here-v5.1.sh <path_to_virsorter_output>

####################################################################################################################################################################################

The last script calls on several scripts from the processing_scripts_v2 folder - ensure it is placed in the same directory as /prophage_extract_gbk+tsv_here-v5.sh script.
