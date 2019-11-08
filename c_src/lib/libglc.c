#include "GLC.h"
#include "stdio.h"
#include "stdlib.h"

void savef(char filename[], double *data, int n, int m, int showm)
{
    FILE *outfile;
    outfile = fopen(filename, "w");
    for (int i = 0; i < n; i++)
    {
        for (int ij = 0; ij < showm; ij++)
        {
            fprintf(outfile, "%15.7e    ", *(data + m * i + ij));
        }
        fprintf(outfile, "\n");
    }
    fclose(outfile);
}
