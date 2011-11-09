/* 
 * WHAM -- Weighted Histogram Analysis Method
 *         2 dimensional implementation
 *
 * Reference 1: Computer Physics Communications, 91(1995) 275-282, Benoit Roux
 * Body of code references equation numbers from this paper.
 *
 * This code is nearly identical to the 1D implementation.  I thought about 
 * making it one big program, but this seemed simpler, especially since the 
 * 2D case requires more command line arguments.
 *
 * $Revision$
 * $Author: alan $
 * $Date: 2003/12/12 22:44:02 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>

#include "wham-2d.h"

// Boltzmann's constant in kcal/mol K
#define k_B 0.001982923700

#define COMMAND_LINE "Command line:  wham-2d Px[=0|pi|val] hist_min_x hist_max_x num_bins_x Py[=0|pi|val] hist_min_y hist_max_y num_bins_y tol temperature numpad metadatafile freefile use_mask\n"

double HIST_MAXx,HIST_MINx,BIN_WIDTHx;
double HIST_MAXy,HIST_MINy,BIN_WIDTHy;
double TOL;
double kT;
double **HISTOGRAM;
int  NUM_BINSx, NUM_BINSy;
int PERIODICx, PERIODICy;
double PERIODx, PERIODy;

int main(int argc, char *argv[])
{
/* moved to global
double kT; // temperature
*/
int i,j;
int first;
int have_energy;
char *freefile;
FILE *METAFILE, *FREEFILE;
struct hist_group *hist_group;
double coor[2];
double error;
double **free_ene;
double **prob, **final_prob;
double *final_f;
double sum;
int iteration;
int max_iteration = 100000;
int numpad;
int **mask;
int use_mask;

if (argc != 15)
    {
    printf( COMMAND_LINE );
    exit(-1);
    }

// Print the command line out into the output file
printf("#");
for (i=0; i<argc; i++)
    {
    printf(" %s", argv[i]);
    }
printf("\n");

PERIODICx = parse_periodic(argv[1], &PERIODx);
if (PERIODICx)
    {
    printf("#Turning on periodicity in x with period = %f\n", PERIODx);
    }
    
HIST_MINx = atof(argv[2]);
HIST_MAXx = atof(argv[3]);
NUM_BINSx = atoi(argv[4]);
BIN_WIDTHx = (HIST_MAXx - HIST_MINx) / (double) NUM_BINSx;

// Parse command line arguments

PERIODICy =parse_periodic(argv[5], &PERIODy);
if (PERIODICy)
    {
    printf("#Turning on periodicity in y with period = %f\n", PERIODy);
    }
HIST_MINy = atof(argv[6]);
HIST_MAXy = atof(argv[7]);
NUM_BINSy = atoi(argv[8]);
BIN_WIDTHy = (HIST_MAXy - HIST_MINy) / (double) NUM_BINSy;

TOL = atof(argv[9]);
kT = atof(argv[10]) * k_B;

numpad = atoi(argv[11]);

METAFILE = fopen(argv[12], "r");
if (METAFILE == (FILE *)NULL)
    {
    printf("couldn't open metadatafile %s: %s\n", argv[12], strerror(errno));
    exit(errno);
    }

i = strlen(argv[13]);
freefile = (char *) malloc(i * sizeof(char));
freefile = argv[13];
if (!freefile)
    {
    printf("couldn't allocate space for freefile name: %s\n", strerror(errno));
    exit(errno);
    }

use_mask = atoi(argv[14]);

HISTOGRAM = (double **) malloc(sizeof(double *) * NUM_BINSx);
if (!HISTOGRAM)
    {
    printf("couldn't allocate space for HISTOGRAM: %s\n", strerror(errno));
    exit(errno);
    }
for (i=0; i<NUM_BINSx; i++)
    {
    HISTOGRAM[i] = (double *) malloc(sizeof(double) * NUM_BINSy);
    if (!HISTOGRAM[i])
        {
        printf("couldn't allocate space for HISTOGRAM[%d]: %s\n", 
                i,strerror(errno));
        exit(errno);
        }
    }

prob = (double **) malloc(sizeof(double *) * NUM_BINSx);
if (!prob)
    {
    printf("couldn't allocate space for prob: %s\n", strerror(errno));
    exit(errno);
    }
for (i=0; i<NUM_BINSx; i++)
    {
    prob[i] = (double *) malloc(sizeof(double) * NUM_BINSy);
    if (!prob[i])
        {
        printf("couldn't allocate space for prob[%d]: %s\n", 
                i,strerror(errno));
        exit(errno);
        }
    }

final_prob = (double **) malloc(sizeof(double *) * NUM_BINSx);
if (!final_prob)
    {
    printf("couldn't allocate space for final_prob: %s\n", strerror(errno));
    exit(errno);
    }
for (i=0; i<NUM_BINSx; i++)
    {
    final_prob[i] = (double *) malloc(sizeof(double) * NUM_BINSy);
    if (!final_prob[i])
        {
        printf("couldn't allocate space for final_prob[%d]: %s\n", 
                i,strerror(errno));
        exit(errno);
        }
    }

free_ene = (double **) malloc(sizeof(double *) * NUM_BINSx);
if (!free_ene)
    {
    printf("couldn't allocate space for free_ene: %s\n", strerror(errno));
    exit(errno);
    }
for (i=0; i<NUM_BINSx; i++)
    {
    free_ene[i] = (double *) malloc(sizeof(double) * NUM_BINSy);
    if (!free_ene[i])
        {
        printf("couldn't allocate space for free_ene[%d]: %s\n", 
                i,strerror(errno));
        exit(errno);
        }
    }

// allocate the mask
if (use_mask)
    {
    // allocate memory to store the mask
    mask = (int **) malloc(sizeof(int *) * NUM_BINSx);
    if (!mask)
        {
        printf("couldn't allocate space for mask: %s\n", strerror(errno));
        exit(errno);
        }
    for (i=0; i<NUM_BINSx; i++)
        {
        // Using calloc to ensure that the array is initialized as all zeros
        //mask[i] = (int *) malloc(sizeof(int) * NUM_BINSy);
        mask[i] = (int *) calloc(NUM_BINSy, sizeof(int));
        if (!mask[i])
            {
            printf("couldn't allocate space for mask[%d]: %s\n", 
                    i,strerror(errno));
            exit(errno);
            }
        }
    }

i = get_numwindows(METAFILE);
printf("#Number of windows = %d\n", i);

hist_group = make_hist_group(i);
//printf("From hist_group: %d\n", hist_group->num_windows);

i = read_metadata(METAFILE, hist_group, use_mask, mask);
assert(i == hist_group->num_windows);

// allocate memory to store the final F values, for when we do MC
// bootstrap error analysis

final_f = (double *)malloc(sizeof(double)*hist_group->num_windows);
if (!final_f)
    {
    printf("couldn't allocate space for final_f: %s\n", strerror(errno));
    exit(errno);
    }


// Figure out if we have trajectories at different temperatures.
// Missing temperatures are set to -1 in read_metadata, and since we 
// require that either all trajectories specify a temperature or all 
// trajectories are assumed to be at the WHAM temperature, we only have to
// check one of them

if (hist_group->kT[0] > 0)
    {
    have_energy = 1;
    }
else
    {
    have_energy = 0;
    for (i=0; i<hist_group->num_windows;i++)
        {
        hist_group->kT[i] = kT;
        }
    }

free(HISTOGRAM);

// for each window, zero out the estimated perturbation due to the restraints
for (i=0; i< hist_group->num_windows; i++)
    {
    hist_group->F[i]=0.0; 
    hist_group->F_old[i]=0.0;
    }

// Do the actual WHAM stuff, iterate to self consistency
iteration = 0;
first = 1;
while (! is_converged(hist_group) || first )
    {
    first = 0;
    save_free(hist_group);
    wham_iteration(hist_group, prob, have_energy, use_mask, mask);
    // Dump out some info
    iteration++;
    if (iteration % 10 == 0)
        {
        error = average_diff(hist_group);
        printf("#Iteration %d:  %f\n", iteration, error);
        }
    
    // Periodically dump out the histogram and free energy
    if (iteration % 100 == 0)
        {
        calc_free(free_ene,prob,kT, use_mask, mask);
        for (i=0; i< NUM_BINSx; i++)
            {
            for (j=0; j< NUM_BINSy; j++)
                {
                calc_coor(i,j,coor);
                printf("%f\t%f\t%f\t%f\n", coor[0], coor[1], free_ene[i][j], 
                                        prob[i][j]);
                }
            }

        // Write the bias values to stdout
        printf("# Dumping simulation biases, in the metadata file order \n");
        printf("# Window  F (free energy units)\n");
        for (j=0; j<hist_group->num_windows;j++)
            {
            printf("# %d\t%f\n", j, hist_group->F[j]);
            }
        }
    // Cheesy bailout if we're going on too long
    if (iteration >= max_iteration) 
        {
        printf("Too many iterations: %d\n", iteration);
        break;
        }
    }


// We're done, write out the free energy and histogram

// Write the bias values to stdout
printf("# Dumping simulation biases, in the metadata file order \n");
printf("# Window  F (free energy units)\n");
for (j=0; j<hist_group->num_windows;j++)
    {
    printf("# %d\t%f\n", j, hist_group->F[j]);
    }

calc_free(free_ene, prob,kT, use_mask, mask);

sum = 0.0;
for(i=0; i< NUM_BINSx; i++)
    {
    for (j=0; j < NUM_BINSy; j++)
        {
        sum += prob[i][j];
        }
    }

for(i=0; i< NUM_BINSx; i++)
    {
    for (j=0; j < NUM_BINSy; j++)
        {
        prob[i][j] /= sum;
        final_prob[i][j] = prob[i][j];
        }
    }

FREEFILE = fopen(freefile, "w");
if (!FREEFILE)
    {
    printf("couldn't open %s: %s\n", freefile, strerror(errno));
    printf("dumping free energy and probability to stdout\n");
    for (i=0; i< NUM_BINSx; i++)
        {
        for (j=0; j< NUM_BINSy; j++)
            {
            calc_coor(i,j,coor);
            printf("%f\t%f\t%f\t%f\n", coor[0], coor[1], 
                                       free_ene[i][j], final_prob[i][j]);
            }
        }
    exit(errno);
    }
else
    {
    // TODO: Add header like in the 1D case
    fprintf(FREEFILE, "#X\t\tY\t\tFree\t\tPro\n");
    // leading padded values in x 
    for (i=-numpad; i<0; i++)
        {
        // leading padding values in y
        for (j=-numpad; j<0; j++)
            {
            calc_coor(i,j,coor); 
            fprintf(FREEFILE,"%f\t%f\t%f\t%f\n", coor[0], coor[1],
                    free_ene[NUM_BINSx+i][NUM_BINSy+j], 
                    final_prob[NUM_BINSx+i][NUM_BINSy+j]);
            }
        // center values in y
        for (j=0; j<NUM_BINSy; j++)
            {
            calc_coor(i,j,coor);
            fprintf(FREEFILE,"%f\t%f\t%f\t%f\n", coor[0], coor[1], 
                           free_ene[NUM_BINSx+i][j], 
                           final_prob[NUM_BINSx+i][j]);
            }
        // trailing padding values in y
        for (j=0; j<numpad; j++)
            {
            calc_coor(i,NUM_BINSy+j,coor); 
            fprintf(FREEFILE,"%f\t%f\t%f\t%f\n", coor[0], coor[1],
                           free_ene[NUM_BINSx+i][j], 
                           final_prob[NUM_BINSx+i][j]);
            }
        fprintf(FREEFILE, "\n");
        }
    // center values in x
    for (i=0; i< NUM_BINSx; i++)
        {
        // leading padding values in y
        for (j=-numpad; j<0; j++)
            {
            calc_coor(i,j,coor); 
            fprintf(FREEFILE,"%f\t%f\t%f\t%f\n", coor[0], coor[1],
                           free_ene[i][NUM_BINSy+j], 
                           final_prob[i][NUM_BINSy+j]);
            }
        // center values in y
        for (j=0; j<NUM_BINSy; j++)
            {
            calc_coor(i,j,coor);
            fprintf(FREEFILE,"%f\t%f\t%f\t%f\n", coor[0], coor[1], 
                             free_ene[i][j], 
                             final_prob[i][j]);
            }
        // trailing padding values in y
        for (j=0; j<numpad; j++)
            {
            calc_coor(i,NUM_BINSy+j,coor); 
            fprintf(FREEFILE,"%f\t%f\t%f\t%f\n", coor[0], coor[1],
                             free_ene[i][j],
                             final_prob[i][j]);
            }
        fprintf(FREEFILE, "\n");
        }
    // trailing padding values in x
    for (i=0; i<numpad; i++)
        {
        // leading padding values in y
        for (j=-numpad; j<0; j++)
            {
            calc_coor(NUM_BINSx+i,j,coor); 
            fprintf(FREEFILE,"%f\t%f\t%f\t%f\n", coor[0], coor[1],
                             free_ene[i][NUM_BINSy+j],
                             final_prob[i][NUM_BINSy+j]);
            }
        // center values in y
        for (j=0; j<NUM_BINSy; j++)
            {
            calc_coor(NUM_BINSx+i,j,coor);
            fprintf(FREEFILE,"%f\t%f\t%f\t%f\n", coor[0], coor[1], 
                             free_ene[i][j], 
                             final_prob[i][j]);
            }
        // trailing padding values in y
        for (j=0; j<numpad; j++)
            {
            calc_coor(NUM_BINSx+i,NUM_BINSy+j,coor); 
            fprintf(FREEFILE,"%f\t%f\t%f\t%f\n", coor[0], coor[1],
                             free_ene[i][j], 
                             final_prob[i][j]);
            }
        fprintf(FREEFILE, "\n");
        }
    }

exit(0);
}

/* ********************** end main program ************************** */

int parse_periodic(char *c, double *period)
/* Read a command line argument, passed as c.
 * Return 1 if the argument says we should have periodicity, 0 otherwise.
 * Set period to be the period, or 0 if there is no periodicity.
 */
{
int is_periodic=1;
int len;
int i;

if (toupper(c[0]) != 'P')
    {
    printf( COMMAND_LINE );
    printf("died here: %s \n", c);
    exit(-1);
    }
else
    {
    len = strlen(c);
    if (len == 2)
        {
        *period = DEGREES;  // 360
        }
    else
        {
        c= &(c[3]);
        if (c[0] == '0')  // turn off periodicity
            {
            is_periodic = 0;
            *period = 0.0;
            }
        else if (isalpha(c[0]))
            {
            for (i=0; i<len-1;i++)
                {
                c[i] = toupper(c[i]);
                }
            if (strncmp(c,"PI",2) == 0)
                {
                *period = RADIANS;  // 2 pi
                }
            else 
                {
                printf( COMMAND_LINE );
                exit(-1);
                }
            }
        else
            {
            *period = atof(c);
            }
        }

    }
return is_periodic;
}

/*******************************************************************/

/*
 *  Perform a single WHAM iteration
 */   
void wham_iteration(struct hist_group* hist_group, double **prob, 
                    int have_energy, int use_mask, int **mask)
{
int i,j,k;
double num, denom, bias, bf, coor[2];
// loop over bins of global histogram
for (i=0; i<NUM_BINSx; i++)
    {
    for (k=0;k<NUM_BINSy; k++)
        {
        if (use_mask && !(mask[i][k])) continue;
        calc_coor(i,k,coor);
        num = 0.0;
        denom = 0.0;
        /*
         *   use previous biases to estimate probability
         *   Equation 8 in Reference 1
         */
        //printf("Bin %d:\n", i);
        for (j=0; j<hist_group->num_windows;j++)
            {
            num += (double) get_histval( &(hist_group->hists[j]),i,k);
            bias = calc_bias(hist_group,j,coor);
            //bf = exp((hist_group->F_old[j] - bias) / kT);
            bf = exp((hist_group->F_old[j] - bias) / hist_group->kT[j]);
            if (have_energy)
                {
                denom += (double) hist_group->partition[j] * bf;
                }
            else
                {
                denom += (double) hist_group->hists[j].num_points * bf;
                }
            }
        prob[i][k] = num / denom;
        //printf("#%d\t%d\t%f\t%f\t%f\n", i,k,num, denom, prob[i][k]);
        /*
         *   use new probability to update bias estimate
         *   Equation 9 from Reference 1
         */
        for (j=0; j<hist_group->num_windows;j++)
            {
            bias = calc_bias(hist_group,j,coor);
            bf = exp(-bias/hist_group->kT[j]) * prob[i][k];
            hist_group->F[j] += bf;
            }
        }
    }
/*
 *   take natural log of Equation 9 from Reference 1
 *   because we store F rather than exp(-F/kT)
 */
for (j=0; j<hist_group->num_windows;j++)
    {
    hist_group->F[j] = -hist_group->kT[j] * log(hist_group->F[j]);
    }
// probably unnecessary, but couldn't hurt
for (j=0; j<hist_group->num_windows;j++)
    {
    hist_group->F[j] -= hist_group->F[0];
    }
}
