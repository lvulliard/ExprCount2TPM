#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
	unsigned int i;
	if (argc == 3){
		// Read command line arguments
		char* htseq_path = argv[1];
		char* ensembl_path = argv[2];
		
		// ----- Parsing gene annotations

		FILE* ensembl_file = fopen(ensembl_path, "r");

		// Buffer
		char ensembl[1024];

		// Read number of lines
		unsigned int ENSEMBL_FILE_LENGTH = 0;
		while(fgets(ensembl, 1023, ensembl_file)){
			ENSEMBL_FILE_LENGTH++;
		}

		// Store the start, end and name of each gene
		int g_start[ENSEMBL_FILE_LENGTH];
		int g_end[ENSEMBL_FILE_LENGTH];
		char* g_name[ENSEMBL_FILE_LENGTH];

		// Where to split the strings
		const char sep1[2] = "	";
		const char sep2[2] = "\"";

		// Restart from beginning of file
		rewind(ensembl_file);
		unsigned int i = 0;
		while(fgets(ensembl, 1024, ensembl_file)){
			char *token;

			token = strtok(ensembl, sep1);
			token = strtok(NULL, sep1);
			token = strtok(NULL, sep1);
			token = strtok(NULL, sep1);
			sscanf(token, "%d", &g_start[i]);

			token = strtok(NULL, sep1);
			sscanf(token, "%d", &g_end[i]);

			token = strtok(NULL, sep1);
			token = strtok(NULL, sep1);
			token = strtok(NULL, sep1);
			token = strtok(NULL, sep1);
			token = strtok(token, sep2);
			token = strtok(NULL, sep2);
			g_name[i] = strdup(token);
			
			i++;
		}

		fclose(ensembl_file);

		// ----- Parsing expression file

		FILE* htseq_file = fopen(htseq_path, "r");

		// Buffer
		char htseq[1024];

		// Read number of lines
		unsigned int HTSEQ_FILE_LENGTH = 0;
		while(fgets(htseq, 1023, htseq_file)){
			HTSEQ_FILE_LENGTH++;
		}

		// Store the start, end and name of each gene
		float h_expr[HTSEQ_FILE_LENGTH];
		char* h_name[HTSEQ_FILE_LENGTH]; 

		// Restart from beginning of file
		rewind(htseq_file);
		int j = 0;
		double scaling_factor = 0;
		while(fgets(htseq, 1024, htseq_file)){
			h_name[j] = strdup(strtok(htseq, sep1));
			sscanf(strtok(NULL, sep1), "%f", &h_expr[j]);

			// Find matching gene in annotations
			int match_name = 0;
			for(i = 0; i < ENSEMBL_FILE_LENGTH; i++){
				if (strcmp(h_name[j], g_name[i]) == 0){
					h_expr[j] /= 1 + abs(g_start[i] - g_end[i]);
					scaling_factor += h_expr[j];
					match_name = 1;
					break;
				}
			}
			if (!match_name){
				fprintf(stderr, "No annotation found for %s\n", h_name[j]);
				h_expr[j] = -1;
			}
			j++;
		}

		fclose(htseq_file);

		// Convert to kilobases / million reads
		scaling_factor /= 1000;

		// Output
		for (i = 0; i < HTSEQ_FILE_LENGTH; i++){
			if (h_expr[i] >= 0){
				printf("%s	%f\n", h_name[i], h_expr[i]/scaling_factor);
			}
		}

	}	
	else {
		fprintf( stderr, "Convert to TPM the raw counts of all the genes for which the length is known.\nUsage:\n  ExprCount2TPM <rawCountFile> <biomartFile>\n");
		return(-1);
	}
}
