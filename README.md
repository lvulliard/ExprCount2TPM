# ExprCount2TPM
Convert to TPM the expression raw counts of all the genes for which the length is known.

## Input
The script takes 2 command-line arguments: the path to an htseq read count files and the path to an annotation file from ENSEMBL (including identifier matching htseq file, and genomic positions of start and end of the feature). See script for details of the expected format.

## Output
Text file with feature identifier and TPM expression (tab-separated).

## Disclaimer: script not maintained
This script is provided without any warranty and support.
