#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_LINE 10000
#define MAX_LOCI 35100000
#define MAX_INDS 1000

double **allele_frequencies;
double af[MAX_LOCI];
int af_count = 0;

void read_af_file(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Allele frequency file open failed");
        exit(1);
    }

    char line[MAX_LINE];
    while (fgets(line, sizeof(line), fp)) {
        af[af_count++] = atof(line);
    }

    fclose(fp);
}

void process_gl_file(const char *filename) {
    FILE *in = fopen(filename, "r");
    if (!in) {
        perror("GL file open failed");
        exit(1);
    }

    char outname[1024];
    snprintf(outname, sizeof(outname), "cpntest_%s", filename);
    char *ext = strstr(outname, ".gl");
    if (ext) strcpy(ext, ".txt");

    FILE *out = fopen(outname, "w");
    if (!out) {
        perror("Output file open failed");
        exit(1);
    }

    char line[MAX_LINE];
    int line_af_index = 0;

    while (fgets(line, sizeof(line), in)) {
        char *pos = strstr(line, ":");
        if (!pos) {
            printf("failed to match %s", line);
            continue;
        }

        // Skip to genotype data
        char *data = strchr(pos, ' ');
        if (!data) continue;
        data++; // move past the space

        // Parse genotype likelihoods
        double p = af[line_af_index++];
        double prior[3];
        prior[0] = pow(1 - p, 2);
        prior[1] = 2 * p * (1 - p);
        prior[2] = pow(p, 2);

        char *token = strtok(data, " \t\n");
        while (token) {
            double gl[3], sum = 0.0, gmean = 0.0;

            for (int i = 0; i < 3 && token; i++) {
                double val = atof(token);
                gl[i] = pow(10.0, val / -10.0) * prior[i];
                sum += gl[i];
                token = strtok(NULL, " \t\n");
            }

            for (int i = 0; i < 3; i++) {
                gl[i] /= sum;
                gmean += i * gl[i];
            }

            fprintf(out, "%.5f ", gmean);
        }

        fprintf(out, "\n");
    }

    fclose(in);
    fclose(out);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "USAGE: %s af_file.txt file.gl\n", argv[0]);
        return 1;
    }

    read_af_file(argv[1]);
    process_gl_file(argv[2]);

    return 0;
}

