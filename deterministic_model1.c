/*

Deterministic model based on the first model of Ehlers and Bataillon (2007).
Code by Allan Crossman.

Note that, internally, the genotypes are referred to by their E&B names and numbers.


FUNDAMENTAL PARAMETERS:

-h <value>
	Value of parameter h, the probability (between 0 and 1) that an inconstant reproduces as a cosex.

-S <value>
	Selfing rate for cosexes (between 0 and 1).

-d <value>
	Inbreeding depression; a fixed penalty (between 0 and 1) applied to the products of selfing.

-V <value>
	Viability of YY individuals (between 0 and 1). A value of 0 corresponds to "ancient dioecy". A value of 1 corresponds to "recent dioecy".

--PSatF <value>
	Pollen production by the population required for females to have all ovules fertilised. A value of 1 corresponds to the theoretical pollen output if the entire population was male. Higher values mean more severe pollen limitation. 0 (default) means no pollen limitation at all.

--ppY <value>
	Only implemented in Model 1. Viability of Y pollen (between 0 and 1).


TO RUN A SINGLE SIMULATION:

Ordinarily the program considers a range of values of Q and F and outputs a graph. To override this and use just 1 value of each, use:

	--onerun -Q <value> -F <value>

Or:

	--onerun -K <value> -k <value>

This will then output (as text) the final genotype frequencies.

The relationship between K and Q is: K = (1 / Q) - 1
The relationship between k and F is: k = (1 / F) - 1

One can also use the arguments --pi and --omega if one prefers, where pi = K + 1, and omega = k + 1.


OUTPUT FORMAT:

--gnuplot
	Text output of female frequencies at equilibrium, suitable for Gnuplot to produce a 3D graph.

--oldformat
	Use K and k parameterisation when creating graph, rather than Q and F.

--oldformatlimit <value>
	What value of k and K to use as the maximum for an "oldformat" (K and k) graph (default 4, i.e. produces a graph with axes of 0-4).

--subdivisions <value>
	Size of the output graph in pixels. Also controls the amount of text output in "gnuplot" mode.


MISC:

--pgd
	Start the population in a state of pseudo-gynodioecy, and attempt to invade males into it (rather than the default of starting with dioecy and attempting to invade inconstants into it).

--iterations
	Iterations before assuming equilibrium has been reached (default 10000, is usually sufficient).

*/


#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MODEL 1

#define PGD 1
#define SSD 2
#define DIO 3
#define PAD 4
#define INC 5


int ** result;

// The following values are defaults that can be changed with command-line options.

float h	= 0.5;					// Probability that inconstant male is cosex
float S = 0.0;					// Selfing rate
float d = 0.0;					// Inbreeding depression
float V = 1.0;					// Fitness of YY individuals, relative to XY individuals
float Q = 1.0;					// Cosex production of pollen, relative to male production
float F = 1.0;					// Cosex production of ovules, relative to female production
float PSatF = 0;				// Pollen saturation point for female receivers
float ppY = 1.0;				// Viability of Y pollen

int pgd = 0;					// Start off in PGD and see if males invade? (Instead of starting with DIO)

int subdivisions = 201;			// Width and height of the output graphics file
int gnuplot = 0;				// Output text of female frequencies suitable for GNU plot in 3D mode

int endpoint = 10000;			// Number of iterations to run
float threshold = 0.01;			// What frequency of a genotype is considered enough to count it as surviving
int onerun = 0;					// Just running once with user-specified Q and F parameters?
int oldformat = 0;				// Old style graph of k and K = 0 to 4?
int oldformatlimit = 4;			// Axis size if drawing graph in old (E&B style) format



void parsecommandline (int argc, char * argv[])
{
	int n;

	for (n = 1; n < argc; n++)
	{
		if (
		(strcmp(argv[n], "-H") == 0 || strcmp(argv[n], "-h") == 0 || strcmp(argv[n], "--inconstanth") == 0)
			&&
		(n < argc - 1)
		) {
			h = atof(argv[n + 1]);
			continue;
		}
		
		if (
		(strcmp(argv[n], "-S") == 0 || strcmp(argv[n], "-s") == 0 || strcmp(argv[n], "--selfing") == 0)
			&&
		(n < argc - 1)
		) {
			S = atof(argv[n + 1]);
			continue;
		}
		
		if (
		(strcmp(argv[n], "-D") == 0 || strcmp(argv[n], "-d") == 0 || strcmp(argv[n], "--depression") == 0)
			&&
		(n < argc - 1)
		) {
			d = atof(argv[n + 1]);
			continue;
		}
		
		if (strcmp(argv[n], "--PSatF") == 0 && n < argc - 1)
		{
			PSatF = atof(argv[n + 1]);
			continue;
		}
		
		if (strcmp(argv[n], "--ppY") == 0 && n < argc - 1)
		{
			ppY = atof(argv[n + 1]);
			continue;
		}
		
		if ((strcmp(argv[n], "-V") == 0 || strcmp(argv[n], "-v") == 0) && n < argc - 1)
		{
			V = atof(argv[n + 1]);
			continue;
		}
		
		if (strcmp(argv[n], "--yypenalty") == 0 && n < argc - 1)
		{
			V = 1 - atof(argv[n + 1]);
			continue;
		}
		
		if (strcmp(argv[n], "--ancient") == 0 || strcmp(argv[n], "--ancientdioecy") == 0)
		{
			V = 0;
			continue;
		}
		
		if (strcmp(argv[n], "--recent") == 0 || strcmp(argv[n], "--recentdioecy") == 0)
		{
			V = 1;
			continue;
		}
		
		if ((strcmp(argv[n], "-Q") == 0 || strcmp(argv[n], "-q") == 0) && n < argc - 1)
		{
			Q = atof(argv[n + 1]);
			continue;
		}
		
		if ((strcmp(argv[n], "-F") == 0 || strcmp(argv[n], "-f") == 0) && n < argc - 1)
		{
			F = atof(argv[n + 1]);
			continue;
		}
		
		if ((strcmp(argv[n], "-K") == 0 || strcmp(argv[n], "--malek") == 0) && n < argc - 1)
		{
			Q = 1.0 / (1.0 + atof(argv[n + 1]));
			continue;
		}
		
		if (strcmp(argv[n], "--pi") == 0)
		{
			Q = 1.0 / atof(argv[n + 1]);
			continue;
		}
		
		if ((strcmp(argv[n], "-k") == 0 || strcmp(argv[n], "--femalek") == 0) && n < argc - 1)
		{
			F = 1.0 / (1.0 + atof(argv[n + 1]));
			continue;
		}
		
		if (strcmp(argv[n], "--omega") == 0)
		{
			F = 1.0 / atof(argv[n + 1]);
			continue;
		}
		
		if (strcmp(argv[n], "--threshold") == 0 && n < argc - 1)
		{
			threshold = atof(argv[n + 1]);
			continue;
		}
		
		if (strcmp(argv[n], "--subdivisions") == 0 && n < argc - 1)
		{
			subdivisions = atoi(argv[n + 1]);		// atoi!
			continue;
		}
		
		if (strcmp(argv[n], "--oldformatlimit") == 0 && n < argc - 1)
		{
			oldformatlimit = atoi(argv[n + 1]);		// atoi!
			continue;
		}
		
		if ((strcmp(argv[n], "--endpoint") == 0 || strcmp(argv[n], "--iterations") == 0) && n < argc - 1)
		{
			endpoint = atoi(argv[n + 1]);			// atoi!
			continue;
		}
		
		if (strcmp(argv[n], "--onerun") == 0)
		{
			onerun = 1;
			continue;
		}
		
		if (strcmp(argv[n], "--pgd") == 0)
		{
			pgd = 1;
			continue;
		}
		
		if (strcmp(argv[n], "--oldformat") == 0)
		{
			oldformat = 1;
			continue;
		}
		
		if (strcmp(argv[n], "--gnuplot") == 0)
		{
			gnuplot = 1;
			continue;
		}
		
		if (argv[n][0] == '-' && isdigit(argv[n][1]) == 0)
		{
			printf("Unrecognised option %s\n", argv[n]);
			exit(1);
		}
	}
	
	return;
}

void drawbmp (char * filename, int magnify)
{
	unsigned int headers[13];
	FILE * outfile;
	int extrabytes;
	int paddedsize;
	int x; int y; int n;
	int i; int j;
	int red = 0, green = 0, blue = 0;

	extrabytes = 4 - ((subdivisions * magnify * 3) % 4);	// How many bytes of padding to add to each
															// horizontal line - the size of which must
															// be a multiple of 4 bytes.
	if (extrabytes == 4) extrabytes = 0;

	paddedsize = ((subdivisions * magnify * 3) + extrabytes) * subdivisions * magnify;

	// Headers...
	// Note that the "BM" identifier in bytes 0 and 1 is NOT included in these "headers".

	headers[0]  = paddedsize + 54;			// bfSize (whole file size)
	headers[1]  = 0;						// bfReserved (both)
	headers[2]  = 54;						// bfOffbits
	headers[3]  = 40;						// biSize
	headers[4]  = subdivisions * magnify;	// biWidth
	headers[5]  = subdivisions * magnify;	// biHeight

	// Would have biPlanes and biBitCount in position 6, but they're shorts.
	// It's easier to write them out separately (see below) than pretend
	// they're a single int, especially with endian issues...

	headers[7]  = 0;						// biCompression
	headers[8]  = paddedsize;				// biSizeImage
	headers[9]  = 0;						// biXPelsPerMeter
	headers[10] = 0;						// biYPelsPerMeter
	headers[11] = 0;						// biClrUsed
	headers[12] = 0;						// biClrImportant
	
	outfile = fopen(filename, "wb");
	if (outfile == NULL)
	{
		printf("Failed to create output file!\n");
		exit(1);
	}
	
	//
	// Headers begin...
	// When printing ints and shorts, we write out 1 character at a time to avoid endian issues.
	//

	fprintf(outfile, "BM");

	for (n = 0; n <= 5; n++)
	{
		fprintf(outfile, "%c", headers[n] & 0x000000FF);
		fprintf(outfile, "%c", (headers[n] & 0x0000FF00) >> 8);
		fprintf(outfile, "%c", (headers[n] & 0x00FF0000) >> 16);
		fprintf(outfile, "%c", (headers[n] & (unsigned int) 0xFF000000) >> 24);
	}

	// These next 4 characters are for the biPlanes and biBitCount fields.

	fprintf(outfile, "%c", 1);
	fprintf(outfile, "%c", 0);
	fprintf(outfile, "%c", 24);
	fprintf(outfile, "%c", 0);

	for (n = 7; n <= 12; n++)
	{
		fprintf(outfile, "%c", headers[n] & 0x000000FF);
		fprintf(outfile, "%c", (headers[n] & 0x0000FF00) >> 8);
		fprintf(outfile, "%c", (headers[n] & 0x00FF0000) >> 16);
		fprintf(outfile, "%c", (headers[n] & (unsigned int) 0xFF000000) >> 24);
	}

	//
	// Headers done, now write the data...
	//
	
	for (y = 0; y < subdivisions; y++)		// BMP image format is written from bottom to top...
	{
		for (j = 0; j < magnify; j++)
		{
			for (x = 0; x < subdivisions; x++)
			{
				red = 0; green = 0; blue = 0;
			
				if (result[x][y] == PGD)
				{
					red = 255; green = 127; blue = 127;
				}
				if (result[x][y] == DIO)
				{
					red = 127; green = 0; blue = 255;
				}
				if (result[x][y] == SSD)
				{
					red = 255; green = 255; blue = 0;
				}
				if (result[x][y] == PAD)
				{
					red = 180; green = 180; blue = 255;
				}
				if (result[x][y] == INC)
				{
					red = 255; green = 255; blue = 255;
				}
				
				// Also, it's written in (b,g,r) format...

				for (i = 0; i < magnify; i++)
				{
					fprintf(outfile, "%c", blue);
					fprintf(outfile, "%c", green);
					fprintf(outfile, "%c", red);
				}
			}
			if (extrabytes)		// See above - BMP lines must be of lengths divisible by 4.
			{
				for (n = 1; n <= extrabytes; n++)
				{
					fprintf(outfile, "%c", 0);
				}
			}
		}
	}
	fclose(outfile);
	return;
}

int main (int argc, char * argv[])
{
	// Plant frequencies... (zeroed here just to suppress a -Wall warning)
	float f_AA = 0;				// AA
	float f_Aa = 0;				// Aa
	float f_Aas = 0;			// Aa*
	float f_aa = 0;				// aa
	float f_aas = 0;			// aa*
	float f_asas = 0;			// a*a*
	
	float next_f_AA;
	float next_f_Aa;
	float next_f_Aas;
	float next_f_aa;
	float next_f_aas;
	float next_f_asas;
	
	// Pollen frequencies...
	float p_A;					// A
	float p_a;					// a
	float p_as;					// a*
	
	// Egg frequencies...
	float e_A;					// A
	float e_a;					// a
	float e_as;					// a*
	
	float male = 0;
	float female = 0;
	float inconstant = 0;
	
	float PSatC;				// Pollen saturation point for cosex receivers
	
	float K;
	float k;
	
	float totalpollen;
	float totalplants;
	int n;
	int x;
	int y;
	
	char bmp_filename[1024];
	char txt_filename[1024];
	
	FILE * textfile = NULL;
	
	
	parsecommandline(argc, argv);
	
	result = malloc(subdivisions * sizeof(int*));
	if (result == NULL)
	{
		printf("Out of memory!\n");
		exit(1);
	}
	for (n = 0; n < subdivisions; n++)
	{
		result[n] = malloc(subdivisions * sizeof(int));
		if (result[n] == NULL)
		{
			printf("Out of memory!\n");
			exit(1);
		}
	}
	
	// Print all settings...
	
	printf("\nModel %d\n\n", MODEL);
	
	if (onerun)
	{
		printf("Q = %G (K = %G, pi = %G)\n", Q, (1 / Q) - 1, 1 / Q);
		printf("F = %G (k = %G, \"omega\" = %G)\n\n", F, (1 / F) - 1, 1 / F);
	}
	
	printf("h = %G\n", h);
	printf("Selfing rate = %G\n", S);
	printf("Inbreeding depression = %G\n", d);
	printf("YY viability = %G (YY penalty = %G)\n", V, 1 - V);
	printf("PSatF = %G\n", PSatF);
	printf("ppY = %G\n\n", ppY);
	
	printf("Iterations = %d\n\n", endpoint);
	
	if (onerun == 0)
	{
		printf("Warning: --onerun option not received, therefore program will\n");
		printf("use %d Q and F combinations and produce a graph. This may\n", subdivisions * subdivisions);
		printf("take a long time. If this was not your intention, terminate now.\n\n");
		
		printf("                       %d |\n", oldformat ? oldformatlimit : 1);
		printf("Output format:       %s   |\n", oldformat ? "k" : "F");
		printf("                       0 |\n");
		printf("                          -----\n");
		printf("                          0   %d\n", oldformat ? oldformatlimit : 1);
		printf("                            %s\n\n", oldformat ? "K" : "Q");
	}
	
	// Choose names for .bmp and .txt output files (if needed)...
	
	sprintf(bmp_filename, "model%d_start%s_V%G_S%G_d%G_h%G_PSatF%G_ppY%G.bmp", MODEL, pgd ? "PGD" : "DIO", V, S, d, h, PSatF, ppY);
	sprintf(txt_filename, "model%d_start%s_V%G_S%G_d%G_h%G_PSatF%G_ppY%G.txt", MODEL, pgd ? "PGD" : "DIO", V, S, d, h, PSatF, ppY);
	
	// Open .txt output file (if needed)...
	
	if (onerun == 0 && gnuplot)
	{
		textfile = fopen(txt_filename, "w");
	}
	
	
	for (y = 0; y < subdivisions; y++)
	{
		for (x = 0; x < subdivisions; x++)
		{
			if (onerun == 0)
			{
			
				// Here we map the X,Y coordinates of our output .bmp file onto Q and F parameters...
			
				if (oldformat == 0)
				{
					Q = (float) x / (subdivisions - 1);						// X axis: Q values 0 to 1
					F = (float) y / (subdivisions - 1);						// Y axis: F values 0 to 1
				} else {
					K = ((float) x / (subdivisions - 1)) * oldformatlimit;	// X axis: K values 0 to oldformatlimit
					k = ((float) y / (subdivisions - 1)) * oldformatlimit;	// Y axis: k values 0 to oldformatlimit
					
					// But Q and F are the parameters actually used by the code, so calculate them:
					Q = 1 / (1 + K);
					F = 1 / (1 + k);
				}
			}
	
			// Set start plant frequencies...
			
			if (pgd == 0)				// Start with DIOECY, try inconstant invasion
			{
				f_AA = 0.499;
				f_Aa = 0.499;
				f_Aas = 0.002;
				f_aa = 0;
				f_aas = 0;
				f_asas = 0;
			} else {					// Start with PSEUDO-GYNODIOECY, try male invasion
				f_AA = 0.499;
				f_Aa = 0.002;
				f_Aas = 0.499;
				f_aa = 0;
				f_aas = 0;
				f_asas = 0;
			}
	
			for (n = 0; n < endpoint; n++)
			{
				// Outcrossed pollen frequencies....................................................
				//
				// Here we sum up the 3 types of pollen (containing the 3 alleles) from the various
				// possible sources. We could do this in 3 equations (as in the paper) but it's simpler
				// to consider each source in turn and add to the totals.
				
				p_A = 0;
				p_a = 0;
				p_as = 0;
				
				// From AA pure females (genotype 1)
				;
				
				// From Aa pure males (genotype 2)
				p_A += f_Aa * 0.5;
				p_a += f_Aa * 0.5;
				
				// From Aa* inconstants (genotype 3) as cosexes
				p_A += f_Aas * 0.5 * h * Q;
				p_as += f_Aas * 0.5 * h * Q;
				
				// From Aa* inconstants (genotype 3) as males
				p_A += f_Aas * 0.5 * (1 - h);
				p_as += f_Aas * 0.5 * (1 - h);
				
				// From aa pure males (genotype 4)
				p_a += f_aa;
				
				// From aa* inconstants (genotype 5) as cosexes
				p_a += f_aas * 0.5 * h * Q;
				p_as += f_aas * 0.5 * h * Q;
				
				// From aa* inconstants (genotype 5) as males
				p_a += f_aas * 0.5 * (1 - h);
				p_as += f_aas * 0.5 * (1 - h);
				
				// From a*a* inconstants (genotype 6) as cosexes
				p_as += f_asas * h * Q;
				
				// From a*a* inconstants (genotype 6) as males
				p_as += f_asas * (1 - h);
				
				// Apply Y pollen viability penalty.................................................
				
				p_a *= ppY;
				p_as *= ppY;
				
				// Normalise pollen frequencies to add up to 1......................................
				
				totalpollen = p_A + p_a + p_as;
				if (totalpollen > 0)
				{
					p_A /= totalpollen;
					p_a /= totalpollen;
					p_as /= totalpollen;
				}
				
				// Outcrossed egg frequencies.......................................................
				
				// Calculate pollen required to fertilise a cosex's outcrossing ovules:
				PSatC = PSatF * F * (1 - S);
				
				e_A = 0;
				e_a = 0;
				e_as = 0;
				
				// From AA pure females (genotype 1)
				if (totalpollen >= PSatF)
				{
					e_A += f_AA;
				} else {
					e_A += f_AA * totalpollen / PSatF;
				}
				
				// From Aa pure males (genotype 2)
				;
				
				// From Aa* inconstants (genotype 3) as cosexes
				if (totalpollen >= PSatC)
				{
					e_A += f_Aas * h * 0.5 * (1 - S) * F;
					e_as +=	f_Aas * h * 0.5 * (1 - S) * F;
				} else {
					e_A += f_Aas * h * 0.5 * (1 - S) * F * totalpollen / PSatC;
					e_as +=	f_Aas * h * 0.5 * (1 - S) * F * totalpollen / PSatC;
				}
				
				// From aa pure males (genotype 4)
				;
				
				// From aa* inconstants (genotype 5) as cosexes
				if (totalpollen >= PSatC)
				{
					e_a += f_aas * h * 0.5 * (1 - S) * F;
					e_as += f_aas * h * 0.5 * (1 - S) * F;
				} else {
					e_a += f_aas * h * 0.5 * (1 - S) * F * totalpollen / PSatC;
					e_as += f_aas * h * 0.5 * (1 - S) * F * totalpollen / PSatC;
				}
				
				// From a*a* inconstants (genotype 6) as cosexes
				if (totalpollen >= PSatC)
				{
					e_as += f_asas * h * (1 - S) * F;
				} else {
					e_as += f_asas * h * (1 - S) * F * totalpollen / PSatC;
				}
				
				// WE CANNOT AND MUST NOT NORMALISE THE EGG FREQUENCIES, AS WE HAVEN'T
				// YET CONSIDERED THE SELFED EGGS. BUT WE DON'T NEED TO NORMALISE.
				
				// Plant frequencies from outcrossing...............................................
				
				next_f_AA = p_A * e_A;
				next_f_Aa = p_A * e_a + p_a * e_A;
				next_f_Aas = p_A * e_as + p_as * e_A;
				next_f_aa = p_a * e_a;
				next_f_aas = p_a * e_as + p_as * e_a;
				next_f_asas = p_as * e_as;
				
				// Additional plants from selfing...................................................
				
				// From AA pure females (genotype 1)
				;
				
				// From Aa pure males (genotype 2)
				;
				
				// From Aa* inconstants (genotype 3)
				next_f_AA += f_Aas * (0.5 / (1 + ppY)) * S * (1 - d) * h * F;
				next_f_Aas += f_Aas * 0.5 * S * (1 - d) * h * F;
				next_f_asas += f_Aas * (0.5 * ppY / (1 + ppY)) * S * (1 - d) * h * F;
				
				// Aa* is the only genotype where there is competition between X and Y pollen
				// during selfing and where the ppY factor therefore is relevant...
				
				// Old versions without ppY:
				// next_f_AA += f_Aas * 0.25 * S * (1 - d) * h * F;
				// next_f_Aas += f_Aas * 0.5 * S * (1 - d) * h * F;
				// next_f_asas += f_Aas * 0.25 * S * (1 - d) * h * F;
				
				// From aa pure males (genotype 4)
				;
				
				// From aa* inconstants (genotype 5)
				next_f_aa += f_aas * 0.25 * S * (1 - d) * h * F;
				next_f_aas += f_aas * 0.5 * S * (1 - d) * h * F;
				next_f_asas += f_aas * 0.25 * S * (1 - d) * h * F;
				
				// From a*a* inconstants (genotype 6)
				next_f_asas += f_asas * S * (1 - d) * h * F;
				
				// Apply YY penalty.................................................................
				
				next_f_aa *= V;
				next_f_aas *= V;
				next_f_asas *= V;

				// Copy.............................................................................
				
				f_AA = next_f_AA;
				f_Aa = next_f_Aa;
				f_Aas = next_f_Aas;
				f_aa = next_f_aa;
				f_aas = next_f_aas;
				f_asas = next_f_asas;
					
				// Normalise plant frequencies to add up to 1.......................................
			
				totalplants = f_AA + f_Aa + f_Aas + f_aa + f_aas + f_asas;
				if (totalplants > 0)
				{
					f_AA /= totalplants;
					f_Aa /= totalplants;
					f_Aas /= totalplants;
					f_aa /= totalplants;
					f_aas /= totalplants;
					f_asas /= totalplants;
				}
			}
	
			// Calculate and save results...
			
			female = f_AA;
			male = f_Aa + f_aa;
			inconstant = f_Aas + f_aas + f_asas;
			
			if (male > threshold && female > threshold && inconstant > threshold)
			{
				result[x][y] = SSD;
			} else if (male > threshold && female > threshold) {
				result[x][y] = DIO;
			} else if (female > threshold && inconstant > threshold) {
				result[x][y] = PGD;
			} else if (male > threshold && inconstant > threshold) {
				result[x][y] = PAD;
			} else if (inconstant > threshold) {
				result[x][y] = INC;
			} else {
				result[x][y] = 0;
			}
			
			if (onerun == 0 && gnuplot)
			{
				fprintf(textfile, "%f", female);
				if (x == subdivisions - 1)
				{
					fprintf(textfile, "\n");
				} else {
					fprintf(textfile, "\t");
				}
			}
			
			if (onerun) break;
		}
		if (onerun) break;
	}
	
	if (onerun == 0)
	{
		drawbmp(bmp_filename, 1);
		printf("Saved %s\n", bmp_filename);
	} else {
		printf("Females       Males         Inconstants\n");
		printf("%.6f      %.6f      %.6f\n\n", female, male, inconstant);
	
		printf("Genotype frequencies, as notated by E&B (2007), or C&C (2012):\n\n");
		
		printf("E&B:  AA (1)    Aa (2)    Aa* (3)   aa (4)    aa* (5)   a*a* (6)\n");
		printf("C&C:  mm (1)    Mm (2)    M*m (4)   MM (3)    M*M (5)   M*M* (6)\n");
		printf("      %.6f  %.6f  %.6f  %.6f  %.6f  %.6f\n\n", f_AA, f_Aa, f_Aas, f_aa, f_aas, f_asas);
		
		if (male > threshold && female > threshold && inconstant > threshold)
		{
			printf("Final state: SSD\n");
		} else if (male > threshold && female > threshold) {
			printf("Final state: DIO\n");
		} else if (female > threshold && inconstant > threshold) {
			printf("Final state: PGD\n");
		} else if (male > threshold && inconstant > threshold) {
			printf("Final state: PAD\n");
		} else if (inconstant > threshold) {
			printf("Final state: INC\n");
		} else {
			printf("Final state: ???\n");
		}
	}
	return 0;
}
