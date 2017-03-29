/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * main.c
 * Copyright (C) 2016 kkozlov <mackoel@gmail.com>
 *
 * multiscale_ssa is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * multiscale_ssa is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/*
 * To change to float do:
 * cat main.c | sed -e "s/double/float/g" -e "s/%lf/%f/g" | sed -e "s/g_rand_float/g_rand_double/g"| grep -v EXTENSION > mssa_float.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <glib.h>
#include <glib/gstdio.h>

#define TMAX 1000000
#define LOWBITS 0x330E
#define BYTESIZE 8

#define SECONDS_PER_DAY 86400
#define MAX_SYNCH_ITER 1e2
#define MOLECULES_PER_CONCENTRATION 1e1
#define SYNCH_STEPS_PER_SLOW 20
#define FAST_TIME_MAX 0.50
#define FAST_STEPS_FRAC 0.3

#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num() 0
	#define omp_get_num_threads() 1
#endif

GRand *grand;

int ssa_prepare(int *species, int *stoichiometric_matrix, double *tau, double *parameter, int number_of_species, int number_of_reactions)
{
	int i, j, reaction_number;
	double *propensity;
	double *probability;
	double prop_sum = 0;
	double random = g_rand_double (grand);
	propensity = g_new0(double, number_of_reactions);
	probability = g_new0(double, number_of_reactions);
	for (i = 0; i < number_of_reactions; i++) {
		double prop = 1;
		int zero = 1;
		for (j = 0; j < number_of_species; j++) {
			if (stoichiometric_matrix[i * number_of_species + j] < 0) {
				zero = 0;
				prop *= species[j] / fabs(stoichiometric_matrix[i * number_of_species + j]);
			}
		}
		propensity[i] = (zero == 1) ? 0 : (parameter[i] * prop);
		prop_sum += propensity[i];
	}
	if (prop_sum <= 0.00000000001) {
		g_warning("Sum of propensities is too small %g!", prop_sum);
		(*tau) = TMAX;
		g_free(propensity);
		g_free(probability);
		return(-1);
	}
	reaction_number = -1;
	for (i = 0; i < number_of_reactions; i++) {
		probability[i] = 0;
		for (j = 0; j <= i; j++) {
			probability[i] += propensity[j] / prop_sum;
		}
		if (random < probability[i]) {
			reaction_number = i;
			break;
		}
	}
	(*tau) = -log(g_rand_double (grand)) / prop_sum;
	g_free(propensity);
	g_free(probability);
	return (reaction_number);
}

int ssa_step(int *species, int *stoichiometric_matrix, double *tau, double *parameter, int number_of_species, int number_of_reactions, int *matrix_reaction)
{
	int j, reaction_number;
	int reaction_type;
	int addon;
	reaction_number = ssa_prepare(species, stoichiometric_matrix, tau, parameter, number_of_species, number_of_reactions);
	if (reaction_number < 0) return(reaction_number);
	reaction_type = matrix_reaction[reaction_number];
	for (j = 0; j < number_of_species; j++) {
		addon = stoichiometric_matrix[reaction_number * number_of_species + j];
		if (addon < 0 && reaction_type == 1) continue;
		species[j] += addon;
		if (species[j] < 0) {
			g_warning("Specie %d < 0", j);
			species[j] = 0;
		}
	}
	return(reaction_number);
}

void evaluate_promotor_state(int *species, int *M, int *species_promoters_indices, int number_of_species, int number_of_binding_sites, int *type_species)
{
	int i, j;
	for (i = 0; i < number_of_species; i++) {
		int k = species_promoters_indices[i];
		if (k > -1) {
			double sum = 0;
			int inhibitors, activators, free;
			inhibitors = 0;
			activators = 0;
			free = 0;
			for (j = 0; j < number_of_binding_sites; j++) {
				int site = M[k * number_of_binding_sites + j];
				sum += site;
				switch (site) {
					case -1:
						inhibitors++;
					break;
					case 0:
						free++;
					break;
					case 1:
						activators++;
					break;
				}
			}
			switch (type_species[i]) {
				case -1:
/*					species[i] = (sum < 0) ? 1:0;*/
					species[i] = (inhibitors > 0) ? 1:0;
				break;
				case 0:
/*					species[i] = (sum == 0) ? 1:0;*/
					species[i] = (free > 0) ? 1:0;
//					species[i] = (inhibitors == 0 && activators == 0) ? 1:0;
//					species[i] = (inhibitors == 0 && free > 0) ? 1:0;
//					species[i] = (activators == 0 && free > 0) ? 1:0;
					break;
				case 1:
/*					species[i] = (sum > 0) ? 1:0;*/
//					species[i] = (inhibitors < activators) ? 1:0;
//					species[i] = (inhibitors == 0 && 0 < activators) ? 1:0;
					species[i] = (inhibitors == 0 && free < activators) ? 1:0;
//					species[i] = (inhibitors < activators && free < activators) ? 1:0;
				break;
			}
		}
	}
}

typedef struct {
	double t_start;
	double t_end;
	int type; /* propagate (1), add bias (2), divide (0), mitosis (-2), gastrulation (-1) etc */
	int n_nucs;
	int has_data;
	int kounter;
	double *data_protein; /* a table of concentrations/molecule numbers */
	int *solution_protein;
	int *bound_protein;
	double *data_mrna; /* a table of concentrations/molecule numbers */
	int *solution_mrna;
	double *corr;
	double *chisq;
} MSSA_Timeclass;

typedef struct {
	int coordinate;
	int length;
	double energy;
	int tf_index;
	int blocked_fw;
	int blocked_bk;
	int repressed_fw;
	int repressed_bk;
	int status; /* bound = 1, free = 0, blocked = -1 */
} MSSA_site;

typedef struct {
	double *T; /* activation coefficient */
	double *protein_degradation; /* n_target_genes length */
	double *mrna_degradation; /* n_target_genes length */
	double *translation;
	double *transport_mrna;
	double *transport_protein;
	double range;
} MSSA_Parameters;

typedef struct {
	int n_target_genes;
	int n_external_genes;
	int *target_gene_index; /* in the tables */
	int *external_gene_index; /* in the tables */
	int n_tfs; /* number of table columns */
	int n_nucs; /* number of table rows */
	GList *tc_list; /* list of MSSA_Timeclass structures */
	int *n_sites; /* number of sites in the regulatory region for each target gene */
	int *m_sites; /* start number of sites in the regulatory region for each target gene */
	int sum_sites; /* total number of sites */
	MSSA_site ***allele_0; /* sites on the first allele */
	MSSA_site ***allele_1; /* sites on the second allele */
	MSSA_Parameters *parameters;
	int repeat;
	int *tf_index_allele_0;
	int *tf_index_allele_1;
	int *status_allele_0;
	int *status_allele_1;
	double *energy_allele_0;
	double *energy_allele_1;
	int *blocked_fw_allele_0;
	int *blocked_fw_allele_1;
	int *repressed_fw_allele_0;
	int *repressed_fw_allele_1;
	int *blocked_bk_allele_0;
	int *blocked_bk_allele_1;
	int *repressed_bk_allele_0;
	int *repressed_bk_allele_1;
	int fsteps;
	int ssteps;
	int number_of_reactions_fast_per_nuc;
	int number_of_reactions_per_nuc;
	int cl_group_size;
	int repeats;
} MSSA_Problem;

static char*data_file;
static char*out_file;
static char*log_file;
static char*action;
static gboolean verbose = FALSE;
static gboolean scores = FALSE;
static gboolean dryrun = FALSE;
static int repeat;
static int mol_per_conc = MOLECULES_PER_CONCENTRATION;
static double fast_time_max = FAST_TIME_MAX;
static double fast_steps_frac = FAST_STEPS_FRAC;
static int parallel = 0; /* 0 - OpenMP, 1 - OpenMP+OpenCL, 2 - OpenCL */
static int gpu_number = 0;

MSSA_Problem *mssa_read_problem(gchar*filename)
{
	MSSA_Problem *problem = g_new(MSSA_Problem, 1);
	problem->repeats = repeat;
	FILE*fp = g_fopen(filename, "r");
	fscanf(fp, "%*s\n");
	fscanf(fp, "%d", &(problem->n_target_genes));
	problem->target_gene_index = g_new0(int, problem->n_target_genes);
	for (int i = 0; i < problem->n_target_genes; i++) {
		fscanf(fp, "%d", &(problem->target_gene_index[i]));
	}
	fscanf(fp, "%d", &(problem->n_tfs));
	problem->n_external_genes = problem->n_tfs - problem->n_target_genes;
	problem->external_gene_index = g_new0(int, problem->n_external_genes);
	for (int i = 0; i < problem->n_external_genes; i++) {
		fscanf(fp, "%d", &(problem->external_gene_index[i]));
	}
	fscanf(fp, "%d", &(problem->n_nucs));
	int ntc;
	fscanf(fp, "%*s");
	fscanf(fp, "%d", &ntc);
	fscanf(fp, "%*s");
	problem->tc_list = NULL;
	int kounter = 0;
	for (int k = 0; k < ntc; k++) {
		MSSA_Timeclass *tc = g_new(MSSA_Timeclass, 1);
		tc->kounter = kounter++;
		fscanf(fp, "%d %lf %lf %d\n", &(tc->type), &(tc->t_start), &(tc->t_end), &(tc->n_nucs));
		fscanf(fp, "%*s\n");
		tc->data_protein = g_new0(double, tc->n_nucs * problem->n_tfs);
		tc->solution_protein = g_new0(int, tc->n_nucs * problem->n_tfs * problem->repeats);
		tc->bound_protein = g_new0(int, tc->n_nucs * problem->n_tfs * problem->repeats);
		if (tc->type > 0) {
			tc->has_data = 1;
			for (int i = 0; i < tc->n_nucs; i++) {
				for (int j = 0; j < problem->n_tfs; j++) {
					fscanf(fp, "%lf", &(tc->data_protein[i * problem->n_tfs + j]));
				}
			}
			fscanf(fp, "%*s");
		} else {
			tc->has_data = 0;
			tc->type = -tc->type;
		}
		tc->data_mrna = g_new0(double, tc->n_nucs * problem->n_tfs);
		tc->solution_mrna = g_new0(int, tc->n_nucs * problem->n_tfs * problem->repeats);
		tc->corr = g_new0(double, problem->n_target_genes * problem->repeats);
		tc->chisq = g_new0(double, problem->n_target_genes * problem->repeats);
		problem->tc_list = g_list_append(problem->tc_list, tc);
	}
	problem->sum_sites = 0;
	problem->n_sites = g_new0(int, problem->n_target_genes);
	problem->m_sites = g_new0(int, problem->n_target_genes);
	for (int i = 0; i < problem->n_target_genes; i++) {
		problem->m_sites[i] = problem->sum_sites;
		fscanf(fp, "%d", &(problem->n_sites[i]));
		problem->sum_sites += problem->n_sites[i];
	}
	fscanf(fp, "%*s");
/*
	problem->allele_0 = g_new0(MSSA_site**, problem->n_nucs);
	problem->allele_1 = g_new0(MSSA_site**, problem->n_nucs);
*/
	problem->allele_0 = g_new0(MSSA_site**, 1);
	problem->allele_1 = g_new0(MSSA_site**, 1);

	problem->allele_0[0] = g_new0(MSSA_site*, problem->n_target_genes);
	problem->allele_1[0] = g_new0(MSSA_site*, problem->n_target_genes);
	for (int i = 0; i < problem->n_target_genes; i++) {
		problem->allele_0[0][i] = g_new0(MSSA_site, problem->n_sites[i]);
		problem->allele_1[0][i] = g_new0(MSSA_site, problem->n_sites[i]);
		for (int j = 0; j < problem->n_sites[i]; j++) {
			fscanf(fp, "%d %*c %d %*f %d %lf", &(problem->allele_0[0][i][j].coordinate),
			       &(problem->allele_0[0][i][j].tf_index), &(problem->allele_0[0][i][j].length),
			       &(problem->allele_0[0][i][j].energy));
			problem->allele_0[0][i][j].status = 0;
			problem->allele_1[0][i][j] = problem->allele_0[0][i][j];
		}
		fscanf(fp, "%*s");
	}

	problem->parameters = g_new(MSSA_Parameters, 1);
	problem->parameters->T = g_new0(double, problem->n_target_genes * problem->n_tfs);
	problem->parameters->protein_degradation = g_new0(double, problem->n_target_genes);
	problem->parameters->mrna_degradation = g_new0(double, problem->n_target_genes);
	problem->parameters->translation = g_new0(double, problem->n_target_genes);
	problem->parameters->transport_mrna = g_new0(double, problem->n_target_genes);
	problem->parameters->transport_protein = g_new0(double, problem->n_target_genes);
	for (int i = 0; i < problem->n_target_genes; i++) {
		for (int j = 0; j < problem->n_tfs; j++) {
			fscanf(fp, "%lf", &(problem->parameters->T[i * problem->n_tfs + j]));
		}
		fscanf(fp, "%lf", &(problem->parameters->translation[i]));
		fscanf(fp, "%lf", &(problem->parameters->protein_degradation[i]));
		fscanf(fp, "%lf", &(problem->parameters->mrna_degradation[i]));
		fscanf(fp, "%lf", &(problem->parameters->transport_mrna[i]));
		fscanf(fp, "%lf", &(problem->parameters->transport_protein[i]));
	}
	fscanf(fp, "%lf", &(problem->parameters->range));
	problem->fsteps = (int)(problem->sum_sites * fast_steps_frac);
	problem->ssteps = 100;
	problem->number_of_reactions_fast_per_nuc = 4 * problem->n_tfs * problem->n_target_genes; /* bind/unbind each tf in each promotor */
	problem->number_of_reactions_per_nuc = 2 * problem->n_target_genes + /* transcription */
		problem->n_target_genes + /* translation */
		problem->n_target_genes + /* degradation mrna */
		problem->n_target_genes + /* degradation protein */
		problem->n_target_genes + /* left transport mrna */
		problem->n_target_genes + /* right transport mrna */
		problem->n_target_genes + /* left transport protein */
		problem->n_target_genes; /* right transport protein */
	return(problem);
}

void setup_problem(MSSA_Problem*problem)
{
/*
  allocate and initialize CPU memory
*/
	int m;
	problem->tf_index_allele_0 = g_new0(int, problem->sum_sites);
	problem->tf_index_allele_1 = g_new0(int, problem->sum_sites);
	problem->status_allele_0 = g_new0(int, problem->sum_sites * problem->n_nucs * problem->repeats);
	problem->status_allele_1 = g_new0(int, problem->sum_sites * problem->n_nucs * problem->repeats);
	problem->energy_allele_0 = g_new0(double, problem->sum_sites);
	problem->energy_allele_1 = g_new0(double, problem->sum_sites);
	problem->blocked_fw_allele_0 = g_new0(int, problem->sum_sites);
	problem->blocked_fw_allele_1 = g_new0(int, problem->sum_sites);
	problem->repressed_fw_allele_0 = g_new0(int, problem->sum_sites);
	problem->repressed_fw_allele_1 = g_new0(int, problem->sum_sites);
	problem->blocked_bk_allele_0 = g_new0(int, problem->sum_sites);
	problem->blocked_bk_allele_1 = g_new0(int, problem->sum_sites);
	problem->repressed_bk_allele_0 = g_new0(int, problem->sum_sites);
	problem->repressed_bk_allele_1 = g_new0(int, problem->sum_sites);

	for (int i = 0; i < problem->n_target_genes; i++) {
		m = problem->m_sites[i];
		for (int j = 0; j < problem->n_sites[i]; j++) {
			problem->tf_index_allele_0[m + j] = problem->allele_0[0][i][j].tf_index;
			problem->tf_index_allele_1[m + j] = problem->allele_1[0][i][j].tf_index;
			problem->energy_allele_0[m + j] = problem->allele_0[0][i][j].energy;
			problem->energy_allele_1[m + j] = problem->allele_1[0][i][j].energy;
			problem->blocked_fw_allele_0[m + j] = problem->allele_0[0][i][j].blocked_fw;
			problem->blocked_fw_allele_1[m + j] = problem->allele_1[0][i][j].blocked_fw;
			problem->repressed_fw_allele_0[m + j] = problem->allele_0[0][i][j].repressed_fw;
			problem->repressed_fw_allele_1[m + j] = problem->allele_1[0][i][j].repressed_fw;
			problem->blocked_bk_allele_0[m + j] = problem->allele_0[0][i][j].blocked_bk;
			problem->blocked_bk_allele_1[m + j] = problem->allele_1[0][i][j].blocked_bk;
			problem->repressed_bk_allele_0[m + j] = problem->allele_0[0][i][j].repressed_bk;
			problem->repressed_bk_allele_1[m + j] = problem->allele_1[0][i][j].repressed_bk;
		}
	}
}

/* Stoichiometric matrix for fast reactions
 * reaction	 allele_0[i] allele_1[i] tf[k]
 * bind_tf[j]
 * unbind_tf[j]
 */

/* Stoichiometric matrix for slow reactions
 * reaction
 * traslation
 * transcription
 * degradation
 */

void mssa_print_timeclass (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	for (int l = 0; l < problem->repeats; l++) {
		fprintf(stdout, "%d time %f %f protein:\n", problem->repeat, tc->t_start, tc->t_end);
		for (int i = 0; i < tc->n_nucs; i++) {
			for (int j = 0; j < problem->n_tfs; j++) {
				fprintf(stdout, "%d ", tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
			}
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "time %f %f bound:\n", tc->t_start, tc->t_end);
		for (int i = 0; i < tc->n_nucs; i++) {
			for (int j = 0; j < problem->n_tfs; j++) {
				fprintf(stdout, "%d ", tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
			}
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "time %f %f mrna:\n", tc->t_start, tc->t_end);
		for (int i = 0; i < tc->n_nucs; i++) {
			for (int j = 0; j < problem->n_tfs; j++) {
				fprintf(stdout, "%d ", tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
			}
			fprintf(stdout, "\n");
		}
	}
}

/* Use this func
 * to print the result of propagation,
 * hence end time is printed.
 * Free and bound proteins are summed.
 */

void mssa_out_timeclass (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (out_file == NULL) return;
	FILE*fp = fopen(out_file, "a");
	for (int l = 0; l < problem->repeats; l++) {
		for (int i = 0; i < tc->n_nucs; i++) {
			fprintf(fp, "%d %d %d %6.3f", l, tc->kounter, i, tc->t_end);
			for (int j = 0; j < problem->n_target_genes; j++) {
				int k = problem->target_gene_index[j];
				fprintf(fp, " %d", tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k] + tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k]);
			}
			for (int j = 0; j < problem->n_target_genes; j++) {
				int k = problem->target_gene_index[j];
				fprintf(fp, " %d", tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k]);
			}
			for (int j = 0; j < problem->n_target_genes; j++) {
				int k = problem->target_gene_index[j];
				fprintf(fp, " %d", tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k]);
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
}

void add_bias (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa %d add bias %d\n", problem->repeat, tc->kounter);
#pragma omp parallel for collapse(3) schedule(static) default(none) shared(tc, problem, mol_per_conc)
	for (int l = 0; l < problem->repeats; l++) {
		for (int i = 0; i < tc->n_nucs; i++) {
			for (int j = 0; j < problem->n_target_genes; j++) {
				int k = problem->target_gene_index[j];
				tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k] += (int)(tc->data_mrna[i * problem->n_tfs + k]);
				tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k] += (int)(tc->data_protein[i * problem->n_tfs + k] * mol_per_conc);
			}
		}
	}
	if (verbose) mssa_print_timeclass (tc, problem);
}

void score (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	int n = tc->n_nucs;
	if (scores) printf("multiscale_ssa %d %d chisq ", problem->repeat, tc->kounter);
	for (int l = 0; l < problem->repeats; l++) {
		for (int j = 0; j < problem->n_target_genes; j++) {
			if (scores) printf("g%d =", j);
			double chisq = 0;
			int k = problem->target_gene_index[j];
#pragma omp parallel for schedule(static) default(none) shared(tc, problem, l, k, n, mol_per_conc) reduction(+:chisq)
			for (int i = 0; i < n; i++) {
				double difference = tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k] + tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k] - tc->data_protein[i * problem->n_tfs + k] * mol_per_conc;
				chisq += difference * difference;
			}
			tc->chisq[j] += chisq;
			if (scores) printf("%15.6f ", chisq);
		}
		if (scores) printf("\n");
		if (scores) printf("multiscale_ssa %d %d corr ", l, tc->kounter);
		for (int j = 0; j < problem->n_target_genes; j++) {
			if (scores) printf("g%d =", j);
			int k = problem->target_gene_index[j];
			double sx, sy, mx, my, xx, yy, nx, ny, xy, sxy, cxy;
			sx = sy = mx = my = xx = yy = nx = ny = xy = sxy = cxy = 0;
#pragma omp parallel for schedule(static) default(none) shared(tc, problem, l, k, n) reduction(+:nx) reduction(+:ny) reduction(+:xx) reduction(+:xy) reduction(+:yy)
			for (int i = 0; i < n; i++) {
				double x = tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k] + tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k];
				double y = tc->data_protein[i * problem->n_tfs + k];
				nx += x;
				ny += y;
				xx += x * x;
				yy += y * y;
				xy += x * y;
			}
			sx = (double)n * xx - nx * nx;
			sy = (double)n * yy - ny * ny;
			sxy = (double)n * xy - nx * ny;
			cxy = (sx > 0 && sy > 0) ? sxy / sqrt(sx * sy) : 0;
			tc->corr[j] += (1 - cxy);
			if (scores) printf("%15.6f ", 1 - cxy);
		}
	}
	if (scores) printf("\n");
}

double corr(double *x, double *y, int n, double *mx_ptr, double *my_ptr, double *sx_ptr, double *sy_ptr, double *sxy_ptr)
{
	double sx, sy, mx, my, xx, yy, nx, ny, xy, sxy, cxy;
	sx = sy = mx = my = xx = yy = nx = ny = xy = sxy = cxy = 0;
#pragma omp parallel for schedule(static) default(none) shared(x, y, n) reduction(+:nx) reduction(+:ny) reduction(+:xx) reduction(+:xy) reduction(+:yy)
	for (int i = 0; i < n; i++) {
		nx += x[i];
		ny += y[i];
		xx += x[i] * x[i];
		yy += y[i] * y[i];
		xy += x[i] * y[i];
	}
	mx = nx / (double)n;
	my = ny / (double)n;
	sx = (xx - 2 * mx * nx + mx * mx) / (double)(n - 1);
	sy = (yy - 2 * my * ny + my * my) / (double)(n - 1);
	sxy = (xy - my * nx - mx * ny + mx * my) / (double)(n - 1);
	cxy = sxy / (sx * sy);
	*mx_ptr = mx;
	*my_ptr = my;
	*sx_ptr = sx;
	*sy_ptr = sy;
	*sxy_ptr = sxy;
	return cxy;
}

void print_scores (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	int n = tc->n_nucs;
	if (tc->type == 1 && tc->has_data == 1) {
		printf("multiscale_ssa %d %d corr ", problem->repeats, tc->kounter);
		for (int j = 0; j < problem->n_target_genes; j++) {
			printf("g%d = ", j);
			printf("%.9f ", tc->corr[j]);
		}
		printf("chisq ");
		for (int j = 0; j < problem->n_target_genes; j++) {
			printf("g%d = ", j);
			printf("%.9f ", tc->chisq[j]);
		}
		printf("\n");
	}
}

double drnd(int* seed)
{
    int const a = 16807;
    int const m = 2147483647;
    double const reciprocal_m = 1.0/m;
    *seed = ((long)((*seed) * a))%m;
    double result = (*seed) * reciprocal_m;
    result = (result < 0) ? -result : result;
    return(result);
}

#ifdef OPENCL
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#include "cl-helper.h"

static cl_context ctx;
static cl_command_queue queue;
static cl_int status;
static cl_mem solution_protein;
static cl_mem solution_mrna;
static cl_mem bound_protein;
static cl_kernel knl_propagate;

static cl_mem T;
static cl_mem translation;
static cl_mem transport_mrna;
static cl_mem transport_protein;
static cl_mem protein_degradation;
static cl_mem mrna_degradation;
static cl_mem n_sites;
static cl_mem m_sites;
static cl_mem target_gene_index;

static cl_mem tf_index_allele_0;
static cl_mem tf_index_allele_1;
static cl_mem status_allele_0;
static cl_mem status_allele_1;
static cl_mem energy_allele_0;
static cl_mem energy_allele_1;
static cl_mem blocked_fw_allele_0;
static cl_mem blocked_fw_allele_1;
static cl_mem repressed_fw_allele_0;
static cl_mem repressed_fw_allele_1;
static cl_mem blocked_bk_allele_0;
static cl_mem blocked_bk_allele_1;
static cl_mem repressed_bk_allele_0;
static cl_mem repressed_bk_allele_1;

void setup_device (MSSA_Problem *problem)
{
	problem->cl_group_size = problem->n_nucs;
/*
 allocate device memory
*/
	if (gpu_number == -1) {
		create_context_on(CHOOSE_INTERACTIVELY, CHOOSE_INTERACTIVELY, 0, &ctx, &queue, 0);
	} else {
		create_context_on_gpu(gpu_number, &ctx, &queue, 0);
	}
	if (verbose)
		print_device_info_from_queue(queue);

	solution_protein = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->n_nucs * problem->n_tfs * problem->repeats, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	solution_mrna = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->n_nucs * problem->n_tfs * problem->repeats, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	bound_protein = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->n_nucs * problem->n_tfs * problem->repeats, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");

	T = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(double) * problem->n_tfs * problem->n_target_genes, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	translation = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(double) * problem->n_target_genes, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	transport_mrna = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(double) * problem->n_target_genes, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	transport_protein = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(double) * problem->n_target_genes, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	protein_degradation = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(double) * problem->n_target_genes, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	mrna_degradation = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(double) * problem->n_target_genes, 0, &status);

	n_sites = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->n_target_genes, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	m_sites = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->n_target_genes, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	target_gene_index = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->n_target_genes, 0, &status);

	tf_index_allele_0 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->sum_sites, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	tf_index_allele_1 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->sum_sites, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	status_allele_0 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->sum_sites * problem->n_nucs * problem->repeats, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	status_allele_1 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->sum_sites * problem->n_nucs * problem->repeats, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	energy_allele_0 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(double) * problem->sum_sites, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	energy_allele_1 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(double) * problem->sum_sites, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	blocked_fw_allele_0 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->sum_sites, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	blocked_fw_allele_1 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->sum_sites, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	repressed_fw_allele_0 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->sum_sites, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	repressed_fw_allele_1 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->sum_sites, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	blocked_bk_allele_0 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->sum_sites, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	blocked_bk_allele_1 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->sum_sites, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	repressed_bk_allele_0 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->sum_sites, 0, &status);
	CHECK_CL_ERROR(status, "clCreateBuffer");
	repressed_bk_allele_1 = clCreateBuffer(ctx, CL_MEM_READ_WRITE,
	  sizeof(int) * problem->sum_sites, 0, &status);
/*
 transfer to device
*/
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, T, /*blocking*/ CL_TRUE, /*offset*/ 0,
		problem->n_tfs * problem->n_target_genes * sizeof(double), problem->parameters->T,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, translation, /*blocking*/ CL_TRUE, /*offset*/ 0,
		problem->n_target_genes * sizeof(double), problem->parameters->translation,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, transport_mrna, /*blocking*/ CL_TRUE, /*offset*/ 0,
		problem->n_target_genes * sizeof(double), problem->parameters->transport_mrna,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, transport_protein, /*blocking*/ CL_TRUE, /*offset*/ 0,
		problem->n_target_genes * sizeof(double), problem->parameters->transport_protein,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, protein_degradation, /*blocking*/ CL_TRUE, /*offset*/ 0,
		problem->n_target_genes * sizeof(double), problem->parameters->protein_degradation,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, mrna_degradation, /*blocking*/ CL_TRUE, /*offset*/ 0,
		problem->n_target_genes * sizeof(double), problem->parameters->mrna_degradation,
		0, NULL, NULL));

	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, n_sites, /*blocking*/ CL_TRUE, /*offset*/ 0,
		problem->n_target_genes * sizeof(int), problem->n_sites,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, m_sites, /*blocking*/ CL_TRUE, /*offset*/ 0,
		problem->n_target_genes * sizeof(int), problem->m_sites,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, target_gene_index, /*blocking*/ CL_TRUE, /*offset*/ 0,
		problem->n_target_genes * sizeof(int), problem->target_gene_index,
		0, NULL, NULL));

	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, energy_allele_0, /*blocking*/ CL_TRUE, /*offset*/ 0,
		problem->sum_sites * sizeof(double), problem->energy_allele_0,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, energy_allele_1, /*blocking*/ CL_TRUE, /*offset*/ 0,
		problem->sum_sites * sizeof(double), problem->energy_allele_1,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, tf_index_allele_0, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->sum_sites, problem->tf_index_allele_0,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, tf_index_allele_1, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->sum_sites, problem->tf_index_allele_1,
		0, NULL, NULL));

	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, blocked_fw_allele_0, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->sum_sites, problem->blocked_fw_allele_0,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, blocked_fw_allele_1, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->sum_sites, problem->blocked_fw_allele_1,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, repressed_fw_allele_0, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->sum_sites, problem->repressed_fw_allele_0,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, repressed_fw_allele_1, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->sum_sites, problem->repressed_fw_allele_1,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, blocked_bk_allele_0, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->sum_sites, problem->blocked_bk_allele_0,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, blocked_bk_allele_1, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->sum_sites, problem->blocked_bk_allele_1,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, repressed_bk_allele_0, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->sum_sites, problem->repressed_bk_allele_0,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, repressed_bk_allele_1, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->sum_sites, problem->repressed_bk_allele_1,
		0, NULL, NULL));

/*
 run code on device
*/
	CALL_CL_GUARDED(clFinish, (queue));
	GError*gerror = NULL;
	GString*krn = g_string_new("\n");
	g_string_append_printf(krn, "#pragma OPENCL EXTENSION cl_khr_fp64: enable\n");
	g_string_append_printf(krn, "\n");
/*langdon_2009_CIGPU.pdf*/
	g_string_append_printf(krn, "double drnd_device(int* seed)\n"); // 1 <= *seed < m
	g_string_append_printf(krn, "{\n");
	g_string_append_printf(krn, "    int const a = 16807;\n"); //ie 7**5
	g_string_append_printf(krn, "    int const m = 2147483647;\n"); //ie 2**31-1
	g_string_append_printf(krn, "    double const reciprocal_m = 1.0 / ((double)m);\n");
	g_string_append_printf(krn, "    *seed = ((long)((*seed) * a)) %% m;\n");
	g_string_append_printf(krn, "    double result = ((double)(*seed)) * reciprocal_m;\n");
    g_string_append_printf(krn, "    result = (result < 0) ? -result : result;\n");
	g_string_append_printf(krn, "    return(result);\n");
	g_string_append_printf(krn, "}\n");
	g_string_append_printf(krn, "\n");

	g_string_append_printf(krn, "__kernel void mssa_propagate(\n");
	g_string_append_printf(krn, "	int seed_rng,\n");
	g_string_append_printf(krn, "	__global int *solution_mrna,\n");
	g_string_append_printf(krn, "	__global const int *target_gene_index,\n");
	g_string_append_printf(krn, "	__global const double *protein_degradation,\n");
	g_string_append_printf(krn, "	__global const double *mrna_degradation,\n");
	g_string_append_printf(krn, "	__global const double *translation,\n");
	g_string_append_printf(krn, "	__global const double *transport_mrna,\n");
	g_string_append_printf(krn, "	__global const double *transport_protein,\n");
	g_string_append_printf(krn, "	int interphase,\n");
	g_string_append_printf(krn, "	__global int *solution_protein,\n");
	g_string_append_printf(krn, "	__global int *bound_protein,\n");
	g_string_append_printf(krn, "	__global const int *tf_index_allele_0,\n");
	g_string_append_printf(krn, "	__global const int *tf_index_allele_1,\n");
	g_string_append_printf(krn, "	__global int *status_allele_0,\n");
	g_string_append_printf(krn, "	__global int *status_allele_1,\n");
	g_string_append_printf(krn, "	__global const double *energy_allele_0,\n");
	g_string_append_printf(krn, "	__global const double *energy_allele_1,\n");
	g_string_append_printf(krn, "	__global const int *blocked_fw_allele_0,\n");
	g_string_append_printf(krn, "	__global const int *blocked_fw_allele_1,\n");
	g_string_append_printf(krn, "	__global const int *repressed_fw_allele_0,\n");
	g_string_append_printf(krn, "	__global const int *repressed_fw_allele_1,\n");
	g_string_append_printf(krn, "	__global const int *blocked_bk_allele_0,\n");
	g_string_append_printf(krn, "	__global const int *blocked_bk_allele_1,\n");
	g_string_append_printf(krn, "	__global const int *repressed_bk_allele_0,\n");
	g_string_append_printf(krn, "	__global const int *repressed_bk_allele_1,\n");
	g_string_append_printf(krn, "	__global const int *n_sites,\n");
	g_string_append_printf(krn, "	__global const int *m_sites,\n");
	g_string_append_printf(krn, "	__global const double *T,\n");
	g_string_append_printf(krn, "	int n_nucs,\n");
	g_string_append_printf(krn, "	int n_reps,\n");
	g_string_append_printf(krn, "	double t_start_slow,\n");
	g_string_append_printf(krn, "	double t_stop_slow)\n");
	g_string_append_printf(krn, "{\n");
  	g_string_append_printf(krn, "	int nuc_id = get_local_id(0);\n");
  	g_string_append_printf(krn, "	int nuc_size = get_local_size(0);\n");
  	g_string_append_printf(krn, "	if (nuc_id < 0 || nuc_id >= n_nucs) return;\n");
  	g_string_append_printf(krn, "	int rep_id = get_group_id(0);\n");
  	g_string_append_printf(krn, "	int rep_size =  get_num_groups(0);\n");
  	g_string_append_printf(krn, "	if (rep_id < 0 || rep_id >= n_reps) return;\n");
	g_string_append_printf(krn, "	__local double propensity[%d];\n", problem->number_of_reactions_per_nuc * problem->n_nucs);
	g_string_append_printf(krn, "	__local double propensity_fast[%d];\n", problem->number_of_reactions_fast_per_nuc * problem->n_nucs);
	g_string_append_printf(krn, "	__local int site_tab[%d];\n", problem->number_of_reactions_fast_per_nuc * problem->n_nucs);
	g_string_append_printf(krn, "	__local double prop_sum_loc[%d];\n", problem->cl_group_size); /* should be group size */
	g_string_append_printf(krn, "	__local double tau_slow_loc[%d];\n", problem->cl_group_size); /* should be group size */
	g_string_append_printf(krn, "	int seed = seed_rng + rep_id * n_nucs + nuc_id;\n");
	g_string_append_printf(krn, "	for (int ap = nuc_id; ap < n_nucs; ap += nuc_size) {\n");
	g_string_append_printf(krn, "		int reaction_number;\n");
	g_string_append_printf(krn, "		for (int i = 0; i < %d; i++) {\n", problem->n_target_genes);
	g_string_append_printf(krn, "			for (int j = 0; j < %d; j++) {\n", problem->n_tfs);
	g_string_append_printf(krn, "				reaction_number = i * %d + j;\n", problem->n_tfs);
	g_string_append_printf(krn, "				propensity_fast[ap * %d + reaction_number] = 0;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "				reaction_number = %d * %d + i * %d + j;\n", problem->n_tfs, problem->n_target_genes, problem->n_tfs);
	g_string_append_printf(krn, "				propensity_fast[ap * %d + reaction_number] = 0;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "				reaction_number = 2 * %d * %d + i * %d + j;\n", problem->n_tfs, problem->n_target_genes, problem->n_tfs);
	g_string_append_printf(krn, "				propensity_fast[ap * %d + reaction_number] = 0;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "				reaction_number = 3 * %d * %d + i * %d + j;\n", problem->n_tfs, problem->n_target_genes, problem->n_tfs);
	g_string_append_printf(krn, "				propensity_fast[ap * %d + reaction_number] = 0;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "			}\n");
	g_string_append_printf(krn, "			reaction_number = ap * %d + 0 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "			propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "			reaction_number = ap * %d + 1 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "			propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "			reaction_number = ap * %d + 2 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "			propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "			reaction_number = ap * %d + 3 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "			propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "			reaction_number = ap * %d + 4 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "			propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "			reaction_number = ap * %d + 5 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "			propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "			reaction_number = ap * %d + 6 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "			propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "			reaction_number = ap * %d + 7 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "			propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "			reaction_number = ap * %d + 8 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "			propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "		}\n");
	g_string_append_printf(krn, "	}\n");
	g_string_append_printf(krn, "	barrier(CLK_LOCAL_MEM_FENCE);\n"); /* write all statuses, solution and bound */
	g_string_append_printf(krn, "	for (int rep = rep_id; rep < n_reps; rep += rep_size) {\n");
  	g_string_append_printf(krn, "		double random;\n");
	g_string_append_printf(krn, "		double t_slow;\n");
	g_string_append_printf(krn, "		double tau_slow;\n");
  	g_string_append_printf(krn, "		t_slow = t_start_slow;\n");
	g_string_append_printf(krn, "		for (int t = 0; t < 100000; t++) {\n");
//	g_string_append_printf(krn, "		while (t_slow < t_stop_slow) {\n");
  	g_string_append_printf(krn, "			if (interphase == 1) { /* interphase */\n");
/* start of parallel region */
	g_string_append_printf(krn, "				for (int ap = nuc_id; ap < n_nucs; ap += nuc_size) {\n");
	g_string_append_printf(krn, "					int reaction_number, found;\n");
	g_string_append_printf(krn, "					double t_start_fast;\n");
	g_string_append_printf(krn, "					double t_stop_fast;\n");
	g_string_append_printf(krn, "					double t_fast;\n");
	g_string_append_printf(krn, "					t_start_fast = %f;\n", 0);
	g_string_append_printf(krn, "					t_stop_fast = %f;\n", fast_time_max);
	g_string_append_printf(krn, "					t_fast = t_start_fast;\n");
	g_string_append_printf(krn, "					for (int l = 0; l < %d; l++) {\n", problem->fsteps);
	g_string_append_printf(krn, "						double tau_fast;\n");
	g_string_append_printf(krn, "						int promoter_number, tf;\n");
	g_string_append_printf(krn, "						int site_number;\n");
	g_string_append_printf(krn, "						int reaction_type;\n");
	g_string_append_printf(krn, "						int i, j, k, s, allele;\n");
	g_string_append_printf(krn, "						double aggregate;\n");
	g_string_append_printf(krn, "						double prop_sum = 0;\n");
				/* Binding */
	g_string_append_printf(krn, "						for (i = 0; i < %d; i++) {\n", problem->n_target_genes);
	g_string_append_printf(krn, "							double prop = 0;\n");
	g_string_append_printf(krn, "							s = (int)(drnd_device(&seed) * n_sites[i]);\n");
	g_string_append_printf(krn, "							j = tf_index_allele_0[m_sites[i] + s];\n");
	g_string_append_printf(krn, "							for (k = 0; k < n_sites[i]; k++) {\n");
	g_string_append_printf(krn, "								if (status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[i] + s] == 0) {\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "									prop = solution_protein[rep * n_nucs * %d + ap * %d + j] * energy_allele_0[m_sites[i] + s];\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									reaction_number = i * %d + j;\n", problem->n_tfs);
	g_string_append_printf(krn, "									propensity_fast[ap * %d + reaction_number] += prop;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "									site_tab[ap * %d + reaction_number] = s;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "									prop_sum += prop;\n");
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "								if (status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[i] + s] == 0) {\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "									prop = solution_protein[rep * n_nucs * %d + ap * %d + j] * energy_allele_1[m_sites[i] + s];\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									reaction_number = %d * %d + i * %d + j;\n", problem->n_tfs, problem->n_target_genes, problem->n_tfs);
	g_string_append_printf(krn, "									propensity_fast[ap * %d + reaction_number] += prop;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "									site_tab[ap * %d + reaction_number] = s;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "									prop_sum += prop;\n");
	g_string_append_printf(krn, "								}\n");
				/* Unbinding */
	g_string_append_printf(krn, "								if (status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[i] + s] == 1) {\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "									prop = bound_protein[rep * n_nucs * %d + ap * %d + j] * energy_allele_0[m_sites[i] + s];\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									reaction_number = 2 * %d * %d + i * %d + j;\n", problem->n_tfs, problem->n_target_genes, problem->n_tfs);
	g_string_append_printf(krn, "									propensity_fast[ap * %d + reaction_number] += prop;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "									site_tab[ap * %d + reaction_number] = s;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "									prop_sum += prop;\n");
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "								if (status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[i] + s] == 1) {\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "									prop = bound_protein[rep * n_nucs * %d + ap * %d + j] * energy_allele_1[m_sites[i] + s];\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									reaction_number = 3 * %d * %d + i * %d + j;\n", problem->n_tfs, problem->n_target_genes, problem->n_tfs);
	g_string_append_printf(krn, "									propensity_fast[ap * %d + reaction_number] += prop;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "									site_tab[ap * %d + reaction_number] = s;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "									prop_sum += prop;\n");
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "								s++;\n");
	g_string_append_printf(krn, "								if (s > n_sites[i] - 1) {\n");
	g_string_append_printf(krn, "									s = 0;\n");
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						}\n");
	g_string_append_printf(krn, "						found = -1;\n");
	g_string_append_printf(krn, "						aggregate = 0;\n");
	g_string_append_printf(krn, "						random = prop_sum * drnd_device(&seed);\n");
	g_string_append_printf(krn, "						for (i = 0; i < %d; i++) {\n", problem->n_target_genes);
				/* Binding */
	g_string_append_printf(krn, "							for (j = 0; j < %d; j++) {\n", problem->n_tfs);
	g_string_append_printf(krn, "								reaction_number = i * %d + j;\n", problem->n_tfs);
	g_string_append_printf(krn, "								aggregate += propensity_fast[ap * %d + reaction_number];\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "								if (random < aggregate) {\n");
	g_string_append_printf(krn, "									found = reaction_number;\n");
	g_string_append_printf(krn, "									promoter_number = i;\n");
	g_string_append_printf(krn, "									tf = j;\n");
	g_string_append_printf(krn, "									reaction_type = 1;\n");
	g_string_append_printf(krn, "									allele = 0;\n");
	g_string_append_printf(krn, "									break;\n");
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "								reaction_number = %d * %d + i * %d + j;\n", problem->n_tfs, problem->n_target_genes, problem->n_tfs);
	g_string_append_printf(krn, "								aggregate += propensity_fast[ap * %d + reaction_number];\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "								if (random < aggregate) {\n");
	g_string_append_printf(krn, "									found = reaction_number;\n");
	g_string_append_printf(krn, "									promoter_number = i;\n");
	g_string_append_printf(krn, "									tf = j;\n");
	g_string_append_printf(krn, "									reaction_type = 1;\n");
	g_string_append_printf(krn, "									allele = 1;\n");
	g_string_append_printf(krn, "									break;\n");
	g_string_append_printf(krn, "								}\n");
				/* Unbinding */
	g_string_append_printf(krn, "								reaction_number = 2 * %d * %d + i * %d + j;\n", problem->n_tfs, problem->n_target_genes, problem->n_tfs);
	g_string_append_printf(krn, "								aggregate += propensity_fast[ap * %d + reaction_number];\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "								if (random < aggregate) {\n");
	g_string_append_printf(krn, "									found = reaction_number;\n");
	g_string_append_printf(krn, "									promoter_number = i;\n");
	g_string_append_printf(krn, "									tf = j;\n");
	g_string_append_printf(krn, "									reaction_type = 0;\n");
	g_string_append_printf(krn, "									allele = 0;\n");
	g_string_append_printf(krn, "									break;\n");
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "								reaction_number = 3 * %d * %d + i * %d + j;\n", problem->n_tfs, problem->n_target_genes, problem->n_tfs);
	g_string_append_printf(krn, "								aggregate += propensity_fast[ap * %d + reaction_number];\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "								if (random < aggregate) {\n");
	g_string_append_printf(krn, "									found = reaction_number;\n");
	g_string_append_printf(krn, "									promoter_number = i;\n");
	g_string_append_printf(krn, "									tf = j;\n");
	g_string_append_printf(krn, "									reaction_type = 0;\n");
	g_string_append_printf(krn, "									allele = 1;\n");
	g_string_append_printf(krn, "									break;\n");
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "							if (found != -1) break;\n");
	g_string_append_printf(krn, "						}\n");
	g_string_append_printf(krn, "						tau_fast = -log(drnd_device(&seed)) / prop_sum;\n");
	g_string_append_printf(krn, "						site_number = -1;\n");
	g_string_append_printf(krn, "						if (found != -1) {\n");
	g_string_append_printf(krn, "							int kount_fw;\n");
	g_string_append_printf(krn, "							int kount_bk;\n");
	g_string_append_printf(krn, "							int block;\n");
	g_string_append_printf(krn, "							if (allele == 0) {\n");
	g_string_append_printf(krn, "								if (reaction_type == 1) { /* Binding */\n");
	g_string_append_printf(krn, "									s = site_tab[ap * %d + found];\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "									site_number = s;\n");
	g_string_append_printf(krn, "									status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s] = 1;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "									solution_protein[rep * n_nucs * %d + ap * %d + tf] += -1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									bound_protein[rep * n_nucs * %d + ap * %d + tf] -= -1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									if (solution_protein[rep * n_nucs * %d + ap * %d + tf] < 0) {\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "										solution_protein[rep * n_nucs * %d + ap * %d + tf] = 0;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									}\n");
	g_string_append_printf(krn, "									if (T[promoter_number * %d + tf] < 0) {\n", problem->n_tfs);
															if (problem->parameters->range > 0) {
	g_string_append_printf(krn, "										kount_fw = repressed_fw_allele_0[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										kount_bk = repressed_bk_allele_0[m_sites[promoter_number] + s];\n");
															} else {
	g_string_append_printf(krn, "										kount_fw = blocked_fw_allele_0[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										kount_bk = blocked_bk_allele_0[m_sites[promoter_number] + s];\n");
															}
	g_string_append_printf(krn, "										block = -1;\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_fw; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_bk; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "									} else {\n");
	g_string_append_printf(krn, "										kount_fw = blocked_fw_allele_0[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										kount_bk = blocked_bk_allele_0[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										block = -1;\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_fw; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_bk; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "									}\n");
	g_string_append_printf(krn, "								} else { /* Unbinding */\n");
	g_string_append_printf(krn, "									s = site_tab[ap * %d + found];\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "									site_number = s;\n");
	g_string_append_printf(krn, "									status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s] = 0;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "									solution_protein[rep * n_nucs * %d + ap * %d + tf] += 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									bound_protein[rep * n_nucs * %d + ap * %d + tf] -= 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									if (bound_protein[rep * n_nucs * %d + ap * %d + tf] < 0) {\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "										bound_protein[rep * n_nucs * %d + ap * %d + tf] = 0;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									}\n");
	g_string_append_printf(krn, "									if (T[promoter_number * %d + tf] < 0) {\n", problem->n_tfs);
															if (problem->parameters->range > 0) {
	g_string_append_printf(krn, "										kount_fw = repressed_fw_allele_0[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										kount_bk = repressed_bk_allele_0[m_sites[promoter_number] + s];\n");
															} else {
	g_string_append_printf(krn, "										kount_fw = blocked_fw_allele_0[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										kount_bk = blocked_bk_allele_0[m_sites[promoter_number] + s];\n");
															}
	g_string_append_printf(krn, "										block = 0;\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_fw; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_bk; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "									} else {\n");
	g_string_append_printf(krn, "										kount_fw = blocked_fw_allele_0[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										kount_bk = blocked_bk_allele_0[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										block = 0;\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_fw; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_bk; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "									}\n");
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "							if (allele == 1) {\n");
	g_string_append_printf(krn, "								if (reaction_type == 1) { /* Binding */\n");
	g_string_append_printf(krn, "									s = site_tab[ap * %d + found];\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "									site_number = s;\n");
	g_string_append_printf(krn, "									status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s] = 1;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "									solution_protein[rep * n_nucs * %d + ap * %d + tf] += -1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									bound_protein[rep * n_nucs * %d + ap * %d + tf] -= -1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									if (solution_protein[rep * n_nucs * %d + ap * %d + tf] < 0) {\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "										solution_protein[rep * n_nucs * %d + ap * %d + tf] = 0;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									}\n");
	g_string_append_printf(krn, "									if (T[promoter_number * %d + tf] < 0) {\n", problem->n_tfs);
															if (problem->parameters->range > 0) {
	g_string_append_printf(krn, "										kount_fw = repressed_fw_allele_1[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										kount_bk = repressed_bk_allele_1[m_sites[promoter_number] + s];\n");
															} else {
	g_string_append_printf(krn, "										kount_fw = blocked_fw_allele_1[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										kount_bk = blocked_bk_allele_1[m_sites[promoter_number] + s];\n");
															}
	g_string_append_printf(krn, "										block = -1;\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_fw; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_bk; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "									} else {\n");
	g_string_append_printf(krn, "										kount_fw = blocked_fw_allele_1[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										kount_bk = blocked_bk_allele_1[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										block = -1;\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_fw; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_bk; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "									}\n");
	g_string_append_printf(krn, "								} else { /* Unbinding */\n");
	g_string_append_printf(krn, "									s = site_tab[ap * %d + found];\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "									site_number = s;\n");
	g_string_append_printf(krn, "									status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s] = 0;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "									solution_protein[rep * n_nucs * %d + ap * %d + tf] += 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									bound_protein[rep * n_nucs * %d + ap * %d + tf] -= 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									if (bound_protein[rep * n_nucs * %d + ap * %d + tf] < 0) {\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "										bound_protein[rep * n_nucs * %d + ap * %d + tf] = 0;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									}\n");
	g_string_append_printf(krn, "									if (T[promoter_number * %d + tf] < 0) {\n", problem->n_tfs);
															if (problem->parameters->range > 0) {
	g_string_append_printf(krn, "										kount_fw = repressed_fw_allele_1[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										kount_bk = repressed_bk_allele_1[m_sites[promoter_number] + s];\n");
															} else {
	g_string_append_printf(krn, "										kount_fw = blocked_fw_allele_1[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										kount_bk = blocked_bk_allele_1[m_sites[promoter_number] + s];\n");
															}
	g_string_append_printf(krn, "										block = 0;\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_fw; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_bk; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "									} else {\n");
	g_string_append_printf(krn, "										kount_fw = blocked_fw_allele_1[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										kount_bk = blocked_bk_allele_1[m_sites[promoter_number] + s];\n");
	g_string_append_printf(krn, "										block = 0;\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_fw; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s + 1 + l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "										for (int l = 0; l < kount_bk; l++) {\n");
	g_string_append_printf(krn, "											if (status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] != 1)\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "												status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[promoter_number] + s - 1 - l] = block;\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "										}\n");
	g_string_append_printf(krn, "									}\n");
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						}\n");
	g_string_append_printf(krn, "						for (i = 0; i < %d; i++) {\n", problem->n_target_genes);
	g_string_append_printf(krn, "							for (j = 0; j < %d; j++) {\n", problem->n_tfs);
	g_string_append_printf(krn, "								reaction_number = i * %d + j;\n", problem->n_tfs);
	g_string_append_printf(krn, "								propensity_fast[ap * %d + reaction_number] = 0;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "								reaction_number = %d * %d + i * %d + j;\n", problem->n_tfs, problem->n_target_genes, problem->n_tfs);
	g_string_append_printf(krn, "								propensity_fast[ap * %d + reaction_number] = 0;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "								reaction_number = 2 * %d * %d + i * %d + j;\n", problem->n_tfs, problem->n_target_genes, problem->n_tfs);
	g_string_append_printf(krn, "								propensity_fast[ap * %d + reaction_number] = 0;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "								reaction_number = 3 * %d * %d + i * %d + j;\n", problem->n_tfs, problem->n_target_genes, problem->n_tfs);
	g_string_append_printf(krn, "								propensity_fast[ap * %d + reaction_number] = 0;\n", problem->number_of_reactions_fast_per_nuc);
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						}\n");
	g_string_append_printf(krn, "						t_fast += tau_fast;\n");

	g_string_append_printf(krn, "						barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n"); /* write all statuses, solution and bound */

	g_string_append_printf(krn, "						if (t_fast > t_stop_fast) break;\n");
	g_string_append_printf(krn, "					}\n");
	g_string_append_printf(krn, "				}\n"); /* endfor nuclei */
  	g_string_append_printf(krn, "			}\n"); /* endif interphase */

  		g_string_append_printf(krn, "			if (nuc_id == 0) {\n");

	g_string_append_printf(krn, "			for (int ap = 0; ap < n_nucs; ap++) {\n");
	g_string_append_printf(krn, "				prop_sum_loc[ap] = 0;\n");
	g_string_append_printf(krn, "			}\n");
//	g_string_append_printf(krn, "			barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n"); /* write all statuses, solution and bound */
	g_string_append_printf(krn, "			for (int ap = 0; ap < n_nucs; ap++) {\n");
//	g_string_append_printf(krn, "			for (int ap = nuc_id; ap < n_nucs; ap += nuc_size) {\n");
	g_string_append_printf(krn, "				for (int i = 0; i < %d; i++) {\n", problem->n_target_genes);
	g_string_append_printf(krn, "					if (interphase == 1) {\n");
/* transcription */
	g_string_append_printf(krn, "						int reaction_number = ap * %d + i;\n", problem->number_of_reactions_per_nuc);
	g_string_append_printf(krn, "						double prop = 0; /* Product of T of bound bs */\n");
	g_string_append_printf(krn, "						int zero = 1; /* flag */\n");
	g_string_append_printf(krn, "						int nact = 0;\n");
	g_string_append_printf(krn, "						int nrep = 0;\n");
/* Allele_0 */
	g_string_append_printf(krn, "						for (int k = 0; k < n_sites[i]; k++) {\n");
	g_string_append_printf(krn, "							if (status_allele_0[rep * n_nucs * %d + ap * %d + m_sites[i] + k] == 1) {\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "								prop += T[i * %d + tf_index_allele_0[m_sites[i] + k]];\n", problem->n_tfs);
	g_string_append_printf(krn, "								zero = 0;\n");
	g_string_append_printf(krn, "								nact += (T[i * %d + tf_index_allele_0[m_sites[i] + k]] > 0) ? 1:0;\n", problem->n_tfs);
	g_string_append_printf(krn, "								nrep += (T[i * %d + tf_index_allele_0[m_sites[i] + k]] < 0) ? 1:0;\n", problem->n_tfs);
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						}\n");
	g_string_append_printf(krn, "						propensity[reaction_number] = (nact > 0) ? exp(prop) : 0;\n");
	g_string_append_printf(krn, "						prop_sum_loc[nuc_id] += propensity[reaction_number];\n");
/* Allele_1 */
	g_string_append_printf(krn, "						reaction_number = ap * %d + %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "						prop = 0; /* Product of T of bound bs */\n");
	g_string_append_printf(krn, "						zero = 1; /* flag */\n");
	g_string_append_printf(krn, "						nact = 0;\n");
	g_string_append_printf(krn, "						nrep = 0;\n");
	g_string_append_printf(krn, "						for (int k = 0; k < n_sites[i]; k++) {\n");
	g_string_append_printf(krn, "							if (status_allele_1[rep * n_nucs * %d + ap * %d + m_sites[i] + k] == 1) {\n", problem->sum_sites, problem->sum_sites);
	g_string_append_printf(krn, "								prop += T[i * %d + tf_index_allele_1[m_sites[i] + k]];\n", problem->n_tfs);
	g_string_append_printf(krn, "								zero = 0;\n");
	g_string_append_printf(krn, "								nact += (T[i * %d + tf_index_allele_1[m_sites[i] + k]] > 0) ? 1:0;\n", problem->n_tfs);
	g_string_append_printf(krn, "								nrep += (T[i * %d + tf_index_allele_1[m_sites[i] + k]] < 0) ? 1:0;\n", problem->n_tfs);
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						}\n");
	g_string_append_printf(krn, "						propensity[reaction_number] = (nact > 0) ? exp(prop) : 0;\n");
	g_string_append_printf(krn, "						prop_sum_loc[nuc_id] += propensity[reaction_number];\n");
/* translation */
	g_string_append_printf(krn, "						reaction_number = ap * %d + 2 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "						int k = target_gene_index[i];\n");
	g_string_append_printf(krn, "						propensity[reaction_number] = solution_mrna[rep * n_nucs * %d + ap * %d + k] * translation[i];\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "						prop_sum_loc[nuc_id] += propensity[reaction_number];\n");
	g_string_append_printf(krn, "					}\n"); /* end if interphase */
	g_string_append_printf(krn, "					int k = target_gene_index[i];\n");
	g_string_append_printf(krn, "					int reaction_number = ap * %d + 3 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
/* mrna degradation */
	g_string_append_printf(krn, "					propensity[reaction_number] = solution_mrna[rep * n_nucs * %d + ap * %d + k] * mrna_degradation[i];\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "					prop_sum_loc[nuc_id] += propensity[reaction_number];\n");
	g_string_append_printf(krn, "					reaction_number = ap * %d + 4 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
/* protein degradation */
	g_string_append_printf(krn, "					propensity[reaction_number] = solution_protein[rep * n_nucs * %d + ap * %d + k] * protein_degradation[i];\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "					prop_sum_loc[nuc_id] += propensity[reaction_number];\n");
	g_string_append_printf(krn, "					double q;\n");
	g_string_append_printf(krn, "					if (ap > 0) {\n");
	g_string_append_printf(krn, "						reaction_number = ap * %d + 5 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "						q = (solution_mrna[rep * n_nucs * %d + (ap - 1) * %d + k] - solution_mrna[rep * n_nucs * %d + ap * %d + k]) * transport_mrna[i];\n", problem->n_tfs, problem->n_tfs, problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "						propensity[reaction_number] = (q > 0) ? q : 0;\n");
	g_string_append_printf(krn, "						prop_sum_loc[nuc_id] += propensity[reaction_number];\n");
/* left protein */
	g_string_append_printf(krn, "						reaction_number = ap * %d + 7 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "						q = (solution_protein[rep * n_nucs * %d + (ap - 1) * %d + k] - solution_protein[rep * n_nucs * %d + ap * %d + k]) * transport_protein[i];\n", problem->n_tfs, problem->n_tfs, problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "						propensity[reaction_number] = (q > 0) ? q : 0;\n");
	g_string_append_printf(krn, "						prop_sum_loc[nuc_id] += propensity[reaction_number];\n");
	g_string_append_printf(krn, "					}\n");
	g_string_append_printf(krn, "					if (ap < n_nucs - 1) {\n");
/* right mrna*/
	g_string_append_printf(krn, "						reaction_number = ap * %d + 6 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "						q = (solution_mrna[rep * n_nucs * %d + (ap + 1) * %d + k] - solution_mrna[rep * n_nucs * %d + ap * %d + k]) * transport_mrna[i];\n", problem->n_tfs, problem->n_tfs, problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "						propensity[reaction_number] = (q > 0) ? q : 0;\n");
	g_string_append_printf(krn, "						prop_sum_loc[nuc_id] += propensity[reaction_number];\n");
/* right protein */
	g_string_append_printf(krn, "						reaction_number = ap * %d + 8 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "						q = (solution_protein[rep * n_nucs * %d + (ap + 1) * %d + k] - solution_protein[rep * n_nucs * %d + ap * %d + k]) * transport_protein[i];\n", problem->n_tfs, problem->n_tfs, problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "						propensity[reaction_number] = (q > 0) ? q : 0;\n");
	g_string_append_printf(krn, "						prop_sum_loc[nuc_id] += propensity[reaction_number];\n");
	g_string_append_printf(krn, "					}\n");

//	g_string_append_printf(krn, "					barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n");

	g_string_append_printf(krn, "				}\n");
	g_string_append_printf(krn, "			}\n");
//	g_string_append_printf(krn, "			barrier(CLK_LOCAL_MEM_FENCE);\n");
//	g_string_append_printf(krn, "			double prop_sum_all = 0;\n");
//	g_string_append_printf(krn, "			for (int ap = 0; ap < %d; ap++) {\n", problem->cl_group_size);
//	g_string_append_printf(krn, "					prop_sum_all += prop_sum_loc[ap];\n");
//	g_string_append_printf(krn, "			}\n");
/*
	g_string_append_printf(krn, "			if (prop_sum_all <= 0.00000000001) {\n");
	g_string_append_printf(krn, "				tau_slow = %f;\n", TMAX);
	g_string_append_printf(krn, "				break;\n");
	g_string_append_printf(krn, "			}\n");
*/
//	g_string_append_printf(krn, "			if (nuc_id == 0) {\n");
	g_string_append_printf(krn, "			double prop_sum_all = 0;\n");
	g_string_append_printf(krn, "			for (int ap = 0; ap < n_nucs; ap++) {\n");
	g_string_append_printf(krn, "					prop_sum_all += prop_sum_loc[ap];\n");
	g_string_append_printf(krn, "			}\n");

	g_string_append_printf(krn, "				int reaction_index = -1; /* local reaction type */\n");
	g_string_append_printf(krn, "				int reaction_target = -1; /* local reaction target */\n");
	g_string_append_printf(krn, "				int reaction_nuc = -1; /* local reaction ap */\n");
	g_string_append_printf(krn, "				double aggregate = 0;\n");
	g_string_append_printf(krn, "				random = prop_sum_all * drnd_device(&seed);\n");
	g_string_append_printf(krn, "				for (int ap = 0; ap < n_nucs; ap++) {\n");
	g_string_append_printf(krn, "					for (int i = 0; i < %d; i++) {\n", problem->n_target_genes);
	g_string_append_printf(krn, "						if (interphase == 1) {\n");
/* transcription 1 */
	g_string_append_printf(krn, "							int reaction_number = ap * %d + i;\n", problem->number_of_reactions_per_nuc);
	g_string_append_printf(krn, "							aggregate += propensity[reaction_number];\n");
	g_string_append_printf(krn, "							if (random < aggregate) {\n");
	g_string_append_printf(krn, "								reaction_nuc = ap;\n");
	g_string_append_printf(krn, "								reaction_target = i;\n");
	g_string_append_printf(krn, "								reaction_index = 1;\n");
	g_string_append_printf(krn, "								break;\n");
	g_string_append_printf(krn, "							}\n");
/* transcription 2 */
	g_string_append_printf(krn, "							reaction_number = ap * %d + %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "							aggregate += propensity[reaction_number];\n");
	g_string_append_printf(krn, "							if (random < aggregate) {\n");
	g_string_append_printf(krn, "								reaction_nuc = ap;\n");
	g_string_append_printf(krn, "								reaction_target = i;\n");
	g_string_append_printf(krn, "								reaction_index = 2;\n");
	g_string_append_printf(krn, "								break;\n");
	g_string_append_printf(krn, "							}\n");
/* translation */
	g_string_append_printf(krn, "							reaction_number = ap * %d + 2 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "							aggregate += propensity[reaction_number];\n");
	g_string_append_printf(krn, "							if (random < aggregate) {\n");
	g_string_append_printf(krn, "								reaction_nuc = ap;\n");
	g_string_append_printf(krn, "								reaction_target = i;\n");
	g_string_append_printf(krn, "								reaction_index = 3;\n");
	g_string_append_printf(krn, "								break;\n");
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						}\n");
/* mrna degradation */
	g_string_append_printf(krn, "						int reaction_number = ap * %d + 3 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "						aggregate += propensity[reaction_number];\n");
	g_string_append_printf(krn, "						if (random < aggregate) {\n");
	g_string_append_printf(krn, "							reaction_nuc = ap;\n");
	g_string_append_printf(krn, "							reaction_target = i;\n");
	g_string_append_printf(krn, "							reaction_index = 4;\n");
	g_string_append_printf(krn, "							break;\n");
	g_string_append_printf(krn, "						}\n");
/* protein degradation */
	g_string_append_printf(krn, "						reaction_number = ap * %d + 4 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "						aggregate += propensity[reaction_number];\n");
	g_string_append_printf(krn, "						if (random < aggregate) {\n");
	g_string_append_printf(krn, "							reaction_nuc = ap;\n");
	g_string_append_printf(krn, "							reaction_target = i;\n");
	g_string_append_printf(krn, "							reaction_index = 5;\n");
	g_string_append_printf(krn, "							break;\n");
	g_string_append_printf(krn, "						}\n");
	g_string_append_printf(krn, "						if (ap > 0) {\n");
/* left mrna*/
	g_string_append_printf(krn, "							reaction_number = ap * %d + 5 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "							aggregate += propensity[reaction_number];\n");
	g_string_append_printf(krn, "							if (random < aggregate) {\n");
	g_string_append_printf(krn, "								reaction_nuc = ap;\n");
	g_string_append_printf(krn, "								reaction_target = i;\n");
	g_string_append_printf(krn, "								reaction_index = 6;\n");
	g_string_append_printf(krn, "								break;\n");
	g_string_append_printf(krn, "							}\n");
/* left protein */
	g_string_append_printf(krn, "							reaction_number = ap * %d + 7 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "							aggregate += propensity[reaction_number];\n");
	g_string_append_printf(krn, "							if (random < aggregate) {\n");
	g_string_append_printf(krn, "								reaction_nuc = ap;\n");
	g_string_append_printf(krn, "								reaction_target = i;\n");
	g_string_append_printf(krn, "								reaction_index = 8;\n");
	g_string_append_printf(krn, "								break;\n");
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						}\n");
	g_string_append_printf(krn, "						if (ap < n_nucs - 1) {\n");
/* right mrna*/
	g_string_append_printf(krn, "							reaction_number = ap * %d + 6 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "							aggregate += propensity[reaction_number];\n");
	g_string_append_printf(krn, "							if (random < aggregate) {\n");
	g_string_append_printf(krn, "								reaction_nuc = ap;\n");
	g_string_append_printf(krn, "								reaction_target = i;\n");
	g_string_append_printf(krn, "								reaction_index = 7;\n");
	g_string_append_printf(krn, "								break;\n");
	g_string_append_printf(krn, "							}\n");
/* right protein */
	g_string_append_printf(krn, "							reaction_number = ap * %d + 8 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "							aggregate += propensity[reaction_number];\n");
	g_string_append_printf(krn, "							if (random < aggregate) {\n");
	g_string_append_printf(krn, "								reaction_nuc = ap;\n");
	g_string_append_printf(krn, "								reaction_target = i;\n");
	g_string_append_printf(krn, "								reaction_index = 9;\n");
	g_string_append_printf(krn, "								break;\n");
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						}\n");

//		g_string_append_printf(krn, "					barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n");

	g_string_append_printf(krn, "					}\n");
	g_string_append_printf(krn, "					if (reaction_index != -1) {\n");
	g_string_append_printf(krn, "						break;\n");
	g_string_append_printf(krn, "					}\n");
	g_string_append_printf(krn, "				}\n");
/*
	g_string_append_printf(krn, "				if (reaction_index == -1) {\n");
	g_string_append_printf(krn, "								reaction_nuc = nuc_id;\n");
	g_string_append_printf(krn, "								reaction_target = 0;\n");
	g_string_append_printf(krn, "								reaction_index = 1;\n");
	g_string_append_printf(krn, "				}\n");
*/
//	g_string_append_printf(krn, "				if (reaction_index != -1) {\n");
	g_string_append_printf(krn, "					int k = target_gene_index[reaction_target];\n");
	g_string_append_printf(krn, "					switch (reaction_index) {\n");
	g_string_append_printf(krn, "						case 1: /* transcription 1 */\n");
	g_string_append_printf(krn, "							solution_mrna[rep * n_nucs * %d + reaction_nuc * %d + k] += 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "						break;\n");
	g_string_append_printf(krn, "						case 2: /* transcription 2 */\n");
	g_string_append_printf(krn, "							solution_mrna[rep * n_nucs * %d + reaction_nuc * %d + k] += 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "						break;\n");
	g_string_append_printf(krn, "						case 3: /* translation */\n");
	g_string_append_printf(krn, "							solution_protein[rep * n_nucs * %d + reaction_nuc * %d + k] += 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "						break;\n");
	g_string_append_printf(krn, "						case 4: /* mrna degradation */\n");
	g_string_append_printf(krn, "							solution_mrna[rep * n_nucs * %d + reaction_nuc * %d + k] -= 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "							if (solution_mrna[rep * n_nucs * %d + reaction_nuc * %d + k] < 0) {\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								solution_mrna[rep * n_nucs * %d + reaction_nuc * %d + k] = 0;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						break;\n");
	g_string_append_printf(krn, "						case 5: /* protein degradation */\n");
	g_string_append_printf(krn, "							solution_protein[rep * n_nucs * %d + reaction_nuc * %d + k] -= 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "							if (solution_protein[rep * n_nucs * %d + reaction_nuc * %d + k] < 0) {\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								solution_protein[rep * n_nucs * %d + reaction_nuc * %d + k] = 0;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						break;\n");
	g_string_append_printf(krn, "						case 6: /* left transport mrna */\n");
	g_string_append_printf(krn, "							if ((reaction_nuc - 1) > -1) {\n");
	g_string_append_printf(krn, "								solution_mrna[rep * n_nucs * %d + (reaction_nuc - 1) * %d + k] -= 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								solution_mrna[rep * n_nucs * %d + reaction_nuc * %d + k] += 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								if (solution_mrna[rep * n_nucs * %d + (reaction_nuc - 1) * %d + k] < 0) {\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									solution_mrna[rep * n_nucs * %d + (reaction_nuc - 1) * %d + k] = 0;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "							break;\n");
	g_string_append_printf(krn, "						case 7: /* right transport mrna */\n");
	g_string_append_printf(krn, "							if ((reaction_nuc + 1) < n_nucs) {\n");
	g_string_append_printf(krn, "								solution_mrna[rep * n_nucs * %d + (reaction_nuc + 1) * %d + k] -= 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								solution_mrna[rep * n_nucs * %d + reaction_nuc * %d + k] += 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								if (solution_mrna[rep * n_nucs * %d + (reaction_nuc + 1) * %d + k] < 0) {\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									solution_mrna[rep * n_nucs * %d + (reaction_nuc + 1) * %d + k] = 0;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						break;\n");
	g_string_append_printf(krn, "						case 8: /* left transport protein */\n");
	g_string_append_printf(krn, "							if ((reaction_nuc - 1) > -1) {\n");
	g_string_append_printf(krn, "								solution_protein[rep * n_nucs * %d + (reaction_nuc - 1) * %d + k] -= 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								solution_protein[rep * n_nucs * %d + reaction_nuc * %d + k] += 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								if (solution_protein[rep * n_nucs * %d + (reaction_nuc - 1) * %d + k] < 0) {\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									solution_protein[rep * n_nucs * %d + (reaction_nuc - 1) * %d + k] = 0;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						break;\n");
	g_string_append_printf(krn, "						case 9: /* right transport protein */\n");
	g_string_append_printf(krn, "							if ((reaction_nuc + 1) < n_nucs) {\n");
	g_string_append_printf(krn, "								solution_protein[rep * n_nucs * %d + (reaction_nuc + 1) * %d + k] -= 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								solution_protein[rep * n_nucs * %d + reaction_nuc * %d + k] += 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								if (solution_protein[rep * n_nucs * %d + (reaction_nuc + 1) * %d + k] < 0) {\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "									solution_protein[rep * n_nucs * %d + (reaction_nuc + 1) * %d + k] = 0;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "								}\n");
	g_string_append_printf(krn, "							}\n");
	g_string_append_printf(krn, "						break;\n");
	g_string_append_printf(krn, "						case -1: /* fail */\n");
	g_string_append_printf(krn, "						default: /* fail */\n");
	g_string_append_printf(krn, "							solution_mrna[rep * n_nucs * %d + 0 * %d + k] = 1111;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "						break;\n");
	g_string_append_printf(krn, "					}\n");
//	g_string_append_printf(krn, "				}\n");
	g_string_append_printf(krn, "				tau_slow = -log(drnd_device(&seed)) / prop_sum_all;\n");
	g_string_append_printf(krn, "				for (int ap = 0; ap < %d; ap++) {\n", problem->cl_group_size);
	g_string_append_printf(krn, "					tau_slow_loc[ap] = tau_slow;\n");
	g_string_append_printf(krn, "				}\n");

//	g_string_append_printf(krn, "			barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n");
	g_string_append_printf(krn, "				for (int ap = 0; ap < n_nucs; ap++) {\n");
	g_string_append_printf(krn, "					for (int i = 0; i < %d; i++) {\n", problem->n_target_genes);
		g_string_append_printf(krn, "					int k = target_gene_index[reaction_target];\n");
		g_string_append_printf(krn, "							solution_mrna[rep * n_nucs * %d + ap * %d + k] = nuc_id + 1;\n", problem->n_tfs, problem->n_tfs);
	g_string_append_printf(krn, "			}\n");
	g_string_append_printf(krn, "			}\n");

	g_string_append_printf(krn, "			}\n"); /* if nuc_id == 0 */
	g_string_append_printf(krn, "			for (int ap = nuc_id; ap < n_nucs; ap += nuc_size) {\n");
	g_string_append_printf(krn, "				for (int i = 0; i < %d; i++) {\n", problem->n_target_genes);
	g_string_append_printf(krn, "					int reaction_number = ap * %d + 0 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "					propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "					reaction_number = ap * %d + 1 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "					propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "					reaction_number = ap * %d + 2 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "					propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "					reaction_number = ap * %d + 3 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "					propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "					reaction_number = ap * %d + 4 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "					propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "					reaction_number = ap * %d + 5 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "					propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "					reaction_number = ap * %d + 6 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "					propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "					reaction_number = ap * %d + 7 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "					propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "					reaction_number = ap * %d + 8 * %d + i;\n", problem->number_of_reactions_per_nuc, problem->n_target_genes);
	g_string_append_printf(krn, "					propensity[reaction_number] = 0;\n");
	g_string_append_printf(krn, "				}\n");
	g_string_append_printf(krn, "			}\n");
//	g_string_append_printf(krn, "			barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n");
	g_string_append_printf(krn, "			tau_slow = tau_slow_loc[nuc_id];\n");
/*	g_string_append_printf(krn, "			iter_kounter++;\n");*/
	g_string_append_printf(krn, "			t_slow += tau_slow;\n");
	g_string_append_printf(krn, "			if (t_slow > t_stop_slow) break;\n");
	g_string_append_printf(krn, "		}\n");
	g_string_append_printf(krn, "	}\n"); /* end for repeats */
	g_string_append_printf(krn, "}\n");/* end function */

	if ( !g_file_set_contents("mssa_propagate.cl", krn->str, -1, &gerror) ) {
		g_message("Set file contents failed: mssa_propagate.cl");
		if (gerror) {
			g_warning(gerror->message);
			g_error_free(gerror);
			gerror = NULL;
		}
	}
	knl_propagate = kernel_from_string(ctx, krn->str, "mssa_propagate", NULL);
	g_string_free(krn, TRUE);
}

void propagate_with_transport_5 (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa %d propagate %d\n", problem->repeat, tc->kounter);
	cl_double t_start_slow;
	cl_double t_stop_slow;
	cl_int interphase = 1;//(tc->type == 3) ? 0 : 1;
	cl_int n_nucs = tc->n_nucs;
	cl_int n_reps = problem->repeats;
/* Set simulation time */
	t_start_slow = tc->t_start;
	t_stop_slow = tc->t_end;
/* Simulate */
	cl_int seed_rng = g_rand_int_range(grand, 1, problem->repeats * problem->n_nucs + 1);
	gint64 time1 = g_get_real_time();
	gint64 time2;
/*
 transfer to device
*/
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, solution_protein, /*blocking*/ CL_TRUE, /*offset*/ 0,
		tc->n_nucs * problem->n_tfs * problem->repeats * sizeof(int), tc->solution_protein,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, solution_mrna, /*blocking*/ CL_TRUE, /*offset*/ 0,
		tc->n_nucs * problem->n_tfs * problem->repeats * sizeof(int), tc->solution_mrna,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, bound_protein, /*blocking*/ CL_TRUE, /*offset*/ 0,
		tc->n_nucs * problem->n_tfs * problem->repeats * sizeof(int), tc->bound_protein,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, status_allele_1, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->repeats * problem->sum_sites * problem->n_nucs, problem->status_allele_1,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueWriteBuffer, (
		queue, status_allele_0, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->repeats * problem->sum_sites * problem->n_nucs, problem->status_allele_0,
		0, NULL, NULL));
/*
 run code on device
*/
	CALL_CL_GUARDED(clFinish, (queue));

	SET_32_KERNEL_ARGS(knl_propagate,
		seed_rng, // 1
		solution_mrna, // 2
		target_gene_index, // 3
		protein_degradation, // 4
		mrna_degradation, // 5
		translation, // 6
		transport_mrna, // 7
		transport_protein, // 8
		interphase, // 9
		solution_protein, // 10
		bound_protein, // 11
		tf_index_allele_0, // 12
		tf_index_allele_1, // 13
		status_allele_0, // 14
		status_allele_1, // 15
		energy_allele_0, // 16
		energy_allele_1, // 17
		blocked_fw_allele_0, // 18
		blocked_fw_allele_1, // 19
		repressed_fw_allele_0, // 20
		repressed_fw_allele_1, // 21
		blocked_bk_allele_0, // 22
		blocked_bk_allele_1, // 23
		repressed_bk_allele_0, // 24
		repressed_bk_allele_1, // 25
		n_sites, // 26
		m_sites, // 27
		T, // 28
		n_nucs, // 29
		n_reps, // 30
		t_start_slow, // 31
		t_stop_slow); // 32

	/* the number of work-items that make up a work-group (also referred to as the size of the work-group)
	* the number of work-items that make up a work-group (also referred to as the size of the work-group)
	* The total number of work-items in the work-group must be less than or equal to the CL_DEVICE_MAX_WORK_GROUP_SIZE
	* the number of work-items specified in local_work_size[0],... local_work_size[work_dim - 1] must be less than or equal to the corresponding values specified by CL_DEVICE_MAX_WORK_ITEM_SIZES[0],.... CL_DEVICE_MAX_WORK_ITEM_SIZES[work_dim - 1] */

	size_t ldim[] = { problem->cl_group_size };
//	size_t ldim[] = { 1 };
	/* the number of global work-items: 2^(CL_DEVICE_ADDRESS_BITS = 32) - 1 */

	size_t gdim[] = { ((problem->cl_group_size * problem->repeats + ldim[0] - 1)/ldim[0]) * ldim[0] };
//	size_t gdim[] = { 1 };

	CALL_CL_GUARDED(clEnqueueNDRangeKernel,
		(queue, knl_propagate,
		 /* dimensions: 1 or 2 or 3 */ 1, NULL, gdim, ldim,
		 0, NULL, NULL));
	CALL_CL_GUARDED(clFinish, (queue));
/*
 transfer back & check
*/
	CALL_CL_GUARDED(clEnqueueReadBuffer, (
		queue, solution_protein, /*blocking*/ CL_TRUE, /*offset*/ 0,
		tc->n_nucs * problem->n_tfs * problem->repeats * sizeof(int), tc->solution_protein,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueReadBuffer, (
		queue, solution_mrna, /*blocking*/ CL_TRUE, /*offset*/ 0,
		tc->n_nucs * problem->n_tfs * problem->repeats * sizeof(int), tc->solution_mrna,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueReadBuffer, (
		queue, bound_protein, /*blocking*/ CL_TRUE, /*offset*/ 0,
		tc->n_nucs * problem->n_tfs * problem->repeats * sizeof(int), tc->bound_protein,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueReadBuffer, (
		queue, status_allele_1, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->repeats * problem->sum_sites * problem->n_nucs, problem->status_allele_1,
		0, NULL, NULL));
	CALL_CL_GUARDED(clEnqueueReadBuffer, (
		queue, status_allele_0, /*blocking*/ CL_TRUE, /*offset*/ 0,
		sizeof(int) * problem->repeats * problem->sum_sites * problem->n_nucs, problem->status_allele_0,
		0, NULL, NULL));
	CALL_CL_GUARDED(clFinish, (queue));
	time2 = g_get_real_time();
	double elapsed = (time2 - time1)/G_USEC_PER_SEC;
	printf("tc %d: %f s\n", tc->kounter, elapsed);
	if (verbose) mssa_print_timeclass (tc, problem);
}

void clear_device(MSSA_Problem *problem)
{
/*
 * clean up
*/
	CALL_CL_GUARDED(clReleaseKernel, (knl_propagate ));
	CALL_CL_GUARDED(clReleaseCommandQueue, (queue));
	CALL_CL_GUARDED(clReleaseContext, (ctx));

	g_free(problem->tf_index_allele_0);
	g_free(problem->tf_index_allele_1);
	g_free(problem->status_allele_0);
	g_free(problem->status_allele_1);
	g_free(problem->energy_allele_0);
	g_free(problem->energy_allele_1);
	g_free(problem->blocked_fw_allele_0);
	g_free(problem->blocked_fw_allele_1);
	g_free(problem->repressed_fw_allele_0);
	g_free(problem->repressed_fw_allele_1);
	g_free(problem->blocked_bk_allele_0);
	g_free(problem->blocked_bk_allele_1);
	g_free(problem->repressed_bk_allele_0);
	g_free(problem->repressed_bk_allele_1);

}

#endif


#ifdef DOCL
#define __kernel
#define __local
#define __global
#define CLK_LOCAL_MEM_FENCE 1
#define CLK_GLOBAL_MEM_FENCE 1
void barrier(int arg) {};
int get_local_id(int arg) {return 0;};
int get_local_size(int arg) {return 1;};
int get_global_id(int arg) {return 0;};
int get_global_size(int arg) {return 1;};
int get_group_id(int arg) {return 0;};
int get_num_groups(int arg) {return 1;};

#include "mssa_propagate.cl"

void propagate_with_transport_6 (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa %d propagate %d\n", problem->repeat, tc->kounter);
	int number_of_reactions;
	int slow_only = (tc->type == 3) ? 1 : 0;
	int interphase = 1 - slow_only;
	int seed_rng = 1;
	gint64 time1 = g_get_real_time();
	gint64 time2;
			mssa_propagate(
				seed_rng, // 1
				tc->solution_mrna, // 2
				problem->target_gene_index, // 3
				problem->parameters->protein_degradation, // 4
				problem->parameters->mrna_degradation, // 5
				problem->parameters->translation, // 6
				problem->parameters->transport_mrna, // 7
				problem->parameters->transport_protein, // 8
				interphase, // 9
				tc->solution_protein, // 10
				tc->bound_protein, // 11
				problem->tf_index_allele_0, // 12
				problem->tf_index_allele_1, // 13
				problem->status_allele_0, // 14
				problem->status_allele_1, // 15
				problem->energy_allele_0, // 16
				problem->energy_allele_1, // 17
				problem->blocked_fw_allele_0, // 18
				problem->blocked_fw_allele_1, // 19
				problem->repressed_fw_allele_0, // 20
				problem->repressed_fw_allele_1, // 21
				problem->blocked_bk_allele_0, // 22
				problem->blocked_bk_allele_1, // 23
				problem->repressed_bk_allele_0, // 24
				problem->repressed_bk_allele_1, // 25
				problem->n_sites, // 26
				problem->m_sites, // 27
				problem->parameters->T, // 28
				tc->n_nucs, // 29
				problem->repeats,
				tc->t_start, // 30
				tc->t_end); // 31
	time2 = g_get_real_time();
	double elapsed = (time2 - time1)/G_USEC_PER_SEC;
	printf("tc %d: %f s\n", tc->kounter, elapsed);
	if (verbose) mssa_print_timeclass (tc, problem);
}
#endif

void propagate_with_transport_4 (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa %d propagate %d\n", problem->repeat, tc->kounter);
	double t_start_slow;
	double t_stop_slow;
	int interphase = (tc->type == 3) ? 0 : 1;
	double *propensity_fast;
	double *site_tab;
	if (interphase == 1) { /* interphase */
		propensity_fast = g_new0(double, problem->number_of_reactions_fast_per_nuc * tc->n_nucs * problem->repeats);
		site_tab = g_new0(double, problem->number_of_reactions_fast_per_nuc * tc->n_nucs * problem->repeats);
	}
/* Set simulation time */
	t_start_slow = tc->t_start;
	t_stop_slow = tc->t_end;
/* Simulate */
	int iter_kounter = 0;
	gint64 time1 = g_get_real_time();
	gint64 time2;
	int seedrng = 1;
	double *propensity = g_new0(double, problem->number_of_reactions_per_nuc * tc->n_nucs * problem->repeats);
	int *seedrn = g_new0(int, tc->n_nucs * problem->repeats);

#pragma omp parallel for collapse(2) schedule(static) default(none) shared(problem, tc, seedrn, seedrng)
	for (int rep = 0; rep < problem->repeats; rep++) {
		for (int ap = 0; ap < tc->n_nucs; ap++) {
			seedrn[rep * tc->n_nucs + ap] = seedrng + rep * tc->n_nucs + ap;
		}
	}
/*
 * omp_set_nested(1);
 * OMP_NESTED=1
#pragma omp parallel for schedule(static) default(none) shared(problem, tc, propensity, propensity_fast, site_tab, t_start_slow, t_stop_slow, interphase, seedrn, fast_time_max)
*/
	for (int rep = 0; rep < problem->repeats; rep++) {
		double tau_slow, t_slow;
	  	t_slow = t_start_slow;
	  	while (t_slow < t_stop_slow) {
	  		if (interphase == 1) {
#pragma omp parallel for schedule(static) default(none) shared(problem, tc, propensity_fast, site_tab, seedrn, fast_time_max, rep)
				for (int ap = 0; ap < tc->n_nucs; ap++) {
					double t_start_fast;
					double t_stop_fast;
					double t_fast;
					t_start_fast = 0;
					t_stop_fast = fast_time_max;
					t_fast = t_start_fast;
					for (int l = 0; l < problem->fsteps; l++) {
						double tau_fast;
						int promoter_number, tf;
						int site_number;
						int reaction_type;
						int reaction_number, found;
						int i, j, k, s, sl, allele;
						double aggregate;
						double prop_sum = 0;
						double random;
					/* Binding */
						for (i = 0; i < problem->n_target_genes; i++) {
							s = (int)(drnd(&(seedrn[rep * tc->n_nucs + ap])) * problem->n_sites[i]);
							double prop = 0;
							j = problem->tf_index_allele_0[problem->m_sites[i] + s];
							for (k = 0; k < problem->n_sites[i]; k++) {
								if (problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[i] + s] == 0) {
									prop = tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + j] * problem->energy_allele_0[problem->m_sites[i] + s];
									reaction_number = i * problem->n_tfs + j;
									propensity_fast[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number] += prop;
									site_tab[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number] = s;
									prop_sum += prop;
								}
								if (problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[i] + s] == 0) {
									prop = tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + j] * problem->energy_allele_1[problem->m_sites[i] + s];
									reaction_number = problem->n_tfs * problem->n_target_genes + i * problem->n_tfs + j;
									propensity_fast[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number] += prop;
									site_tab[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number] = s;
									prop_sum += prop;
								}
					/* Unbinding */
								if (problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[i] + s] == 1) {
									prop = tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + j] * problem->energy_allele_0[problem->m_sites[i] + s];
									reaction_number = 2 * problem->n_tfs * problem->n_target_genes + i * problem->n_tfs + j;
									propensity_fast[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number] += prop;
									site_tab[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number] = s;
									prop_sum += prop;
								}
								if (problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[i] + s] == 1) {
									prop = tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + j] * problem->energy_allele_1[problem->m_sites[i] + s];
									reaction_number = 3 * problem->n_tfs * problem->n_target_genes + i * problem->n_tfs + j;
									propensity_fast[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number] += prop;
									site_tab[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number] = s;
									prop_sum += prop;
								}
								s++;
								if (s > problem->n_sites[i] - 1) {
									s = 0;
								}
							}
						}
						if (prop_sum <= 0.00000000001) {
							tau_fast = TMAX;
							promoter_number = -1;
							tf = -1;
							reaction_type = -1;
							site_number = -1;
							break;
						}
						found = -1;
						aggregate = 0;
						random = prop_sum * drnd(&(seedrn[rep * tc->n_nucs + ap]));
						for (i = 0; i < problem->n_target_genes; i++) {
							for (j = 0; j < problem->n_tfs; j++) {
					/* Binding */
								reaction_number = i * problem->n_tfs + j;
								aggregate += propensity_fast[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number];
								if (random < aggregate) {
									found = reaction_number;
									promoter_number = i;
									tf = j;
									reaction_type = 1;
									allele = 0;
									break;
								}
								reaction_number = problem->n_tfs * problem->n_target_genes + i * problem->n_tfs + j;
								aggregate += propensity_fast[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number];
								if (random < aggregate) {
									found = reaction_number;
									promoter_number = i;
									tf = j;
									reaction_type = 1;
									allele = 1;
									break;
								}
				/* Unbinding */
								reaction_number = 2 * problem->n_tfs * problem->n_target_genes + i * problem->n_tfs + j;
								aggregate += propensity_fast[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number];
								if (random < aggregate) {
									found = reaction_number;
									promoter_number = i;
									tf = j;
									reaction_type = 0;
									allele = 0;
									break;
								}
								reaction_number = 3 * problem->n_tfs * problem->n_target_genes + i * problem->n_tfs + j;
								aggregate += propensity_fast[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number];
								if (random < aggregate) {
									found = reaction_number;
									promoter_number = i;
									tf = j;
									reaction_type = 0;
									allele = 1;
									break;
								}
							}
							if (found != -1) break;
						}
						tau_fast = -log(drnd(&(seedrn[rep * tc->n_nucs + ap]))) / prop_sum;
						site_number = -1;
						if (found != -1) {
							int kount_fw;
							int kount_bk;
							int block;
							if (allele == 0) {
								if (reaction_type == 1) { /* Binding */
									s = site_tab[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + found];
									site_number = s;
									problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s] = 1;
									tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] += -1;
									tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] -= -1;
									if (tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] < 0) {
	 									g_warning("fast bind reaction n %d t %d: ap %d tf %d target %d < 0", found, reaction_type, ap, tf, promoter_number);
										tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] = 0;
									}
									if (problem->parameters->T[promoter_number * problem->n_tfs + tf] < 0) {
										if (problem->parameters->range > 0) {
											kount_fw = problem->repressed_fw_allele_0[problem->m_sites[promoter_number] + s];
											kount_bk = problem->repressed_bk_allele_0[problem->m_sites[promoter_number] + s];
										} else {
											kount_fw = problem->blocked_fw_allele_0[problem->m_sites[promoter_number] + s];
											kount_bk = problem->blocked_bk_allele_0[problem->m_sites[promoter_number] + s];
										}
										block = -1;
										for (int l = 0; l < kount_fw; l++) {
											if (problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] != 1)
												problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] = block;
										}
										for (int l = 0; l < kount_bk; l++) {
											if (problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] != 1)
												problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] = block;
										}
									} else {
										kount_fw = problem->blocked_fw_allele_0[problem->m_sites[promoter_number] + s];
										kount_bk = problem->blocked_bk_allele_0[problem->m_sites[promoter_number] + s];
										block = -1;
										for (int l = 0; l < kount_fw; l++) {
											if (problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] != 1)
												problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] = block;
										}
										for (int l = 0; l < kount_bk; l++) {
											if (problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] != 1)
												problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] = block;
										}
									}
								} else { /* Unbinding */
									s = site_tab[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + found];
									site_number = s;
									problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s] = 0;
									tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] += 1;
									tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] -= 1;
									if (tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] < 0) {
	 									g_warning("fast unbind reaction n %d t %d: ap %d tf %d target %d < 0", found, reaction_type, ap, tf, promoter_number);
										tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] = 0;
									}
									if (problem->parameters->T[promoter_number * problem->n_tfs + tf] < 0) {
										if (problem->parameters->range > 0) {
											kount_fw = problem->repressed_fw_allele_0[problem->m_sites[promoter_number] + s];
											kount_bk = problem->repressed_bk_allele_0[problem->m_sites[promoter_number] + s];
										} else {
											kount_fw = problem->blocked_fw_allele_0[problem->m_sites[promoter_number] + s];
											kount_bk = problem->blocked_bk_allele_0[problem->m_sites[promoter_number] + s];
										}
										block = 0;
										for (int l = 0; l < kount_fw; l++) {
											if (problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] != 1)
												problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] = block;
										}
										for (int l = 0; l < kount_bk; l++) {
											if (problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] != 1)
												problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] = block;
										}
									} else {
										kount_fw = problem->blocked_fw_allele_0[problem->m_sites[promoter_number] + s];
										kount_bk = problem->blocked_bk_allele_0[problem->m_sites[promoter_number] + s];
										block = 0;
										for (int l = 0; l < kount_fw; l++) {
											if (problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] != 1)
												problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] = block;
										}
										for (int l = 0; l < kount_bk; l++) {
											if (problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] != 1)
												problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] = block;
										}
									}
								}
							}
							if (allele == 1) {
								if (reaction_type == 1) { /* Binding */
									s = site_tab[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + found];
									site_number = s;
									problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s] = 1;
									tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] += -1;
									tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] -= -1;
									if (tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] < 0) {
	 									g_warning("fast bind reaction n %d t %d: ap %d tf %d target %d < 0", found, reaction_type, ap, tf, promoter_number);
										tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] = 0;
									}
									if (problem->parameters->T[promoter_number * problem->n_tfs + tf] < 0) {
										if (problem->parameters->range > 0) {
											kount_fw = problem->repressed_fw_allele_1[problem->m_sites[promoter_number] + s];
											kount_bk = problem->repressed_bk_allele_1[problem->m_sites[promoter_number] + s];
										} else {
											kount_fw = problem->blocked_fw_allele_1[problem->m_sites[promoter_number] + s];
											kount_bk = problem->blocked_bk_allele_1[problem->m_sites[promoter_number] + s];
										}
										block = -1;
										for (int l = 0; l < kount_fw; l++) {
											if (problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] != 1)
												problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] = block;
										}
										for (int l = 0; l < kount_bk; l++) {
											if (problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] != 1)
												problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] = block;
										}
									} else {
										kount_fw = problem->blocked_fw_allele_1[problem->m_sites[promoter_number] + s];
										kount_bk = problem->blocked_bk_allele_1[problem->m_sites[promoter_number] + s];
										block = -1;
										for (int l = 0; l < kount_fw; l++) {
											if (problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] != 1)
												problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] = block;
										}
										for (int l = 0; l < kount_bk; l++) {
											if (problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] != 1)
												problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] = block;
										}
									}
								} else { /* Unbinding */
									s = site_tab[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + found];
									site_number = s;
									problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s] = 0;
									tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] += 1;
									tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] -= 1;
									if (tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] < 0) {
										g_warning("fast unbind reaction n %d t %d: ap %d tf %d target %d < 0", found, reaction_type, ap, tf, promoter_number);
										tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] = 0;
									}
									if (problem->parameters->T[promoter_number * problem->n_tfs + tf] < 0) {
										if (problem->parameters->range > 0) {
											kount_fw = problem->repressed_fw_allele_1[problem->m_sites[promoter_number] + s];
											kount_bk = problem->repressed_bk_allele_1[problem->m_sites[promoter_number] + s];
										} else {
											kount_fw = problem->blocked_fw_allele_1[problem->m_sites[promoter_number] + s];
											kount_bk = problem->blocked_bk_allele_1[problem->m_sites[promoter_number] + s];
										}
										block = 0;
										for (int l = 0; l < kount_fw; l++) {
											if (problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] != 1)
												problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] = block;
										}
										for (int l = 0; l < kount_bk; l++) {
											if (problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] != 1)
												problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] = block;
										}
									} else {
										kount_fw = problem->blocked_fw_allele_1[problem->m_sites[promoter_number] + s];
										kount_bk = problem->blocked_bk_allele_1[problem->m_sites[promoter_number] + s];
										block = 0;
										for (int l = 0; l < kount_fw; l++) {
											if (problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] != 1)
												problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s + 1 + l] = block;
										}
										for (int l = 0; l < kount_bk; l++) {
											if (problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] != 1)
												problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[promoter_number] + s - 1 - l] = block;
										}
									}
								}
							}
						}
						for (i = 0; i < problem->n_target_genes; i++) {
							for (j = 0; j < problem->n_tfs; j++) {
								reaction_number = i * problem->n_tfs + j;
								propensity_fast[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number] = 0;
								reaction_number = problem->n_tfs * problem->n_target_genes + i * problem->n_tfs + j;
								propensity_fast[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number] = 0;
								reaction_number = 2 * problem->n_tfs * problem->n_target_genes + i * problem->n_tfs + j;
								propensity_fast[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number] = 0;
								reaction_number = 3 * problem->n_tfs * problem->n_target_genes + i * problem->n_tfs + j;
								propensity_fast[rep * problem->number_of_reactions_fast_per_nuc * tc->n_nucs + ap * problem->number_of_reactions_fast_per_nuc + reaction_number] = 0;
							}
						}
						t_fast += tau_fast;
						if (t_fast > t_stop_fast) break;
					}
				} /* endfor nuclei */
	  		} /* endif interphase */

			int reaction_index = -1; /* local reaction type */
			int reaction_target = -1; /* local reaction target */
			int reaction_nuc = -1; /* local reaction ap */
			int reaction_number, i, k, found;
			double aggregate;
			double prop_sum_all = 0;
			double random;

#pragma omp parallel for collapse(2) schedule(static) default(none) shared(problem, tc, propensity, rep, interphase) reduction(+:prop_sum_all)
			for (int ap = 0; ap < tc->n_nucs; ap++) {
				for (int i = 0; i < problem->n_target_genes; i++) {
					if (interphase == 1) {
/* transcription */
						int reaction_number = ap * problem->number_of_reactions_per_nuc + i;
						double prop = 0; /* Product of T of bound bs */
						int zero = 1; /* flag */
						int nact = 0;
						int nrep = 0;
/* Allele_0 */
						for (int k = 0; k < problem->n_sites[i]; k++) {
							if (problem->status_allele_0[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[i] + k] == 1) {
								prop += problem->parameters->T[i * problem->n_tfs + problem->tf_index_allele_0[problem->m_sites[i] + k]];
								zero = 0;
								nact += (problem->parameters->T[i * problem->n_tfs + problem->tf_index_allele_0[problem->m_sites[i] + k]] > 0) ? 1:0;
								nrep += (problem->parameters->T[i * problem->n_tfs + problem->tf_index_allele_0[problem->m_sites[i] + k]] < 0) ? 1:0;
							}
						}
						propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = (nact > 0) ? exp(prop) : 0;
						prop_sum_all += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
/* Allele_1 */
						reaction_number = ap * problem->number_of_reactions_per_nuc + problem->n_target_genes + i;
						prop = 0; /* Product of T of bound bs */
						zero = 1; /* flag */
						nact = 0;
						nrep = 0;
						for (int k = 0; k < problem->n_sites[i]; k++) {
							if (problem->status_allele_1[rep * tc->n_nucs * problem->sum_sites + ap * problem->sum_sites + problem->m_sites[i] + k] == 1) {
								prop += problem->parameters->T[i * problem->n_tfs + problem->tf_index_allele_1[problem->m_sites[i] + k]];
								zero = 0;
								nact += (problem->parameters->T[i * problem->n_tfs + problem->tf_index_allele_1[problem->m_sites[i] + k]] > 0) ? 1:0;
								nrep += (problem->parameters->T[i * problem->n_tfs + problem->tf_index_allele_1[problem->m_sites[i] + k]] < 0) ? 1:0;
							}
						}
						propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = (nact > 0) ? exp(prop) : 0;
						prop_sum_all += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
/* translation */
						reaction_number = ap * problem->number_of_reactions_per_nuc + 2 * problem->n_target_genes + i;
						int k = problem->target_gene_index[i];
						propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + k] * problem->parameters->translation[i];
						prop_sum_all += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
					}
					int k = problem->target_gene_index[i];
					double q;
					int reaction_number = ap * problem->number_of_reactions_per_nuc + 3 * problem->n_target_genes + i;
/* mrna degradation */
					propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + k] * problem->parameters->mrna_degradation[i];
					prop_sum_all += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
					reaction_number = ap * problem->number_of_reactions_per_nuc + 4 * problem->n_target_genes + i;
/* protein degradation */
					propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + k] * problem->parameters->protein_degradation[i];
					prop_sum_all += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
					if (ap > 0) {
						reaction_number = ap * problem->number_of_reactions_per_nuc + 5 * problem->n_target_genes + i;
						q = (tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + (ap - 1) * problem->n_tfs + k] - tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + k]) * problem->parameters->transport_mrna[i];
						propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = (q > 0) ? q : 0;
						prop_sum_all += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
		/* left protein */
						reaction_number = ap * problem->number_of_reactions_per_nuc + 7 * problem->n_target_genes + i;
						q = (tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + (ap - 1) * problem->n_tfs + k] - tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + k]) * problem->parameters->transport_protein[i];
						propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = (q > 0) ? q : 0;
						prop_sum_all += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
					}
					if (ap < tc->n_nucs - 1) {
		/* right mrna*/
						reaction_number = ap * problem->number_of_reactions_per_nuc + 6 * problem->n_target_genes + i;
						q = (tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + (ap + 1) * problem->n_tfs + k] - tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + k]) * problem->parameters->transport_mrna[i];
						propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = (q > 0) ? q : 0;
						prop_sum_all += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
		/* right protein */
						reaction_number = ap * problem->number_of_reactions_per_nuc + 8 * problem->n_target_genes + i;
						q = (tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + (ap + 1) * problem->n_tfs + k] - tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + k]) * problem->parameters->transport_protein[i];
						propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = (q > 0) ? q : 0;
						prop_sum_all += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
					}
				}
			}
			if (prop_sum_all <= 0.00000000001) {
				g_warning("slow Sum of propensities is too small %g!", prop_sum_all);
				tau_slow = TMAX;
				break;
			}
			found = 0;
			aggregate = 0;
			random = prop_sum_all * drnd(&(seedrn[rep * tc->n_nucs]));
			for (int ap = 0; ap < tc->n_nucs; ap++) {
				for (i = 0; i < problem->n_target_genes; i++) {
					if (interphase == 1) {
		/* transcription 1 */
						reaction_number = ap * problem->number_of_reactions_per_nuc + i;
						aggregate += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
						if (random < aggregate) {
							reaction_nuc = ap;
							reaction_target = i;
							reaction_index = 1;
							found = 1;
							break;
						}
		/* transcription 2 */
						reaction_number = ap * problem->number_of_reactions_per_nuc + problem->n_target_genes + i;
						aggregate += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
						if (random < aggregate) {
							reaction_nuc = ap;
							reaction_target = i;
							reaction_index = 2;
							found = 1;
							break;
						}
		/* translation */
						reaction_number = ap * problem->number_of_reactions_per_nuc + 2 * problem->n_target_genes + i;
						aggregate += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
						if (random < aggregate) {
							reaction_nuc = ap;
							reaction_target = i;
							reaction_index = 3;
							found = 1;
							break;
						}
					}
		/* mrna degradation */
					reaction_number = ap * problem->number_of_reactions_per_nuc + 3 * problem->n_target_genes + i;
					aggregate += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
					if (random < aggregate) {
						reaction_nuc = ap;
						reaction_target = i;
						reaction_index = 4;
						found = 1;
						break;
					}
		/* protein degradation */
					reaction_number = ap * problem->number_of_reactions_per_nuc + 4 * problem->n_target_genes + i;
					aggregate += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
					if (random < aggregate) {
						reaction_nuc = ap;
						reaction_target = i;
						reaction_index = 5;
						found = 1;
						break;
					}
					if (ap > 0) {
		/* left mrna*/
						reaction_number = ap * problem->number_of_reactions_per_nuc + 5 * problem->n_target_genes + i;
						aggregate += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
						if (random < aggregate) {
							reaction_nuc = ap;
							reaction_target = i;
							reaction_index = 6;
							found = 1;
							break;
						}
		/* left protein */
						reaction_number = ap * problem->number_of_reactions_per_nuc + 7 * problem->n_target_genes + i;
						aggregate += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
						if (random < aggregate) {
							reaction_nuc = ap;
							reaction_target = i;
							reaction_index = 8;
							found = 1;
							break;
						}
					}
					if (ap < tc->n_nucs - 1) {
		/* right mrna*/
						reaction_number = ap * problem->number_of_reactions_per_nuc + 6 * problem->n_target_genes + i;
						aggregate += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
						if (random < aggregate) {
							reaction_nuc = ap;
							reaction_target = i;
							reaction_index = 7;
							found = 1;
							break;
						}
		/* right protein */
						reaction_number = ap * problem->number_of_reactions_per_nuc + 8 * problem->n_target_genes + i;
						aggregate += propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number];
						if (random < aggregate) {
							reaction_nuc = ap;
							reaction_target = i;
							reaction_index = 9;
							found = 1;
							break;
						}
					}
				}
				if (found == 1) {
					break;
				}
			}
			tau_slow = -log(drnd(&(seedrn[rep * tc->n_nucs + reaction_nuc]))) / prop_sum_all;
			switch (reaction_index) {
				case 1: /* transcription 1 */
					k = problem->target_gene_index[reaction_target];
					tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k]++;
				break;
				case 2: /* transcription 2 */
					k = problem->target_gene_index[reaction_target];
					tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k]++;
				break;
				case 3: /* translation */
					k = problem->target_gene_index[reaction_target];
					tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k]++;
				break;
				case 4: /* mrna degradation */
					k = problem->target_gene_index[reaction_target];
					tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k]--;
					if (tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k] < 0) {
						g_warning("slow mrna: ap %d tf %d < 0", reaction_nuc, k);
						tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k] = 0;
					}
				break;
				case 5: /* protein degradation */
					k = problem->target_gene_index[reaction_target];
					tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k]--;
					if (tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k] < 0) {
						g_warning("slow prot: ap %d tf %d < 0", reaction_nuc, k);
						tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k] = 0;
					}
				break;
				case 6: /* left transport mrna */
					k = problem->target_gene_index[reaction_target];
					if ((reaction_nuc - 1) > -1) {
						tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + (reaction_nuc - 1) * problem->n_tfs + k]--;
						tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k]++;
						if (tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + (reaction_nuc - 1) * problem->n_tfs + k] < 0) {
							g_warning("slow prot: ap %d tf %d < 0", reaction_nuc, k);
							tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + (reaction_nuc - 1) * problem->n_tfs + k] = 0;
						}
					} else {
						g_warning("slow prot: %d ap %d bad %d", reaction_index, reaction_nuc, tc->n_nucs);
					}
				break;
				case 7: /* right transport mrna */
					k = problem->target_gene_index[reaction_target];
					if ((reaction_nuc + 1) < tc->n_nucs) {
						tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + (reaction_nuc + 1) * problem->n_tfs + k]--;
						tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k]++;
						if (tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + (reaction_nuc + 1) * problem->n_tfs + k] < 0) {
							g_warning("slow prot: ap %d tf %d < 0", reaction_nuc, k);
							tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + (reaction_nuc + 1) * problem->n_tfs + k] = 0;
						}
					} else {
						g_warning("slow prot: %d ap %d bad %d", reaction_index, reaction_nuc, tc->n_nucs);
					}
				break;
				case 8: /* left transport protein */
					k = problem->target_gene_index[reaction_target];
					if ((reaction_nuc - 1) > -1) {
						tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + (reaction_nuc - 1) * problem->n_tfs + k]--;
						tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k]++;
						if (tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + (reaction_nuc - 1) * problem->n_tfs + k] < 0) {
							g_warning("slow prot: ap %d tf %d < 0", reaction_nuc, k);
							tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + (reaction_nuc - 1) * problem->n_tfs + k] = 0;
						}
					} else {
						g_warning("slow prot: %d ap %d bad %d", reaction_index, reaction_nuc, tc->n_nucs);
					}
				break;
				case 9: /* right transport protein */
					k = problem->target_gene_index[reaction_target];
					if ((reaction_nuc + 1) < tc->n_nucs) {
						tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + (reaction_nuc + 1) * problem->n_tfs + k]--;
						tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k]++;
						if (tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + (reaction_nuc + 1) * problem->n_tfs + k] < 0) {
							g_warning("slow prot: ap %d tf %d < 0", reaction_nuc, k);
							tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + (reaction_nuc + 1) * problem->n_tfs + k] = 0;
						}
					} else {
						g_warning("slow prot: %d ap %d bad %d", reaction_index, reaction_nuc, tc->n_nucs);
					}
				break;
			}
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(problem, tc, propensity, rep)
			for (int ap = 0; ap < tc->n_nucs; ap++) {
				for (i = 0; i < problem->n_target_genes; i++) {
					int reaction_number = ap * problem->number_of_reactions_per_nuc + 0 * problem->n_target_genes + i;
					propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = 0;
					reaction_number = ap * problem->number_of_reactions_per_nuc + 1 * problem->n_target_genes + i;
					propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = 0;
					reaction_number = ap * problem->number_of_reactions_per_nuc + 2 * problem->n_target_genes + i;
					propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = 0;
					reaction_number = ap * problem->number_of_reactions_per_nuc + 3 * problem->n_target_genes + i;
					propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = 0;
					reaction_number = ap * problem->number_of_reactions_per_nuc + 4 * problem->n_target_genes + i;
					propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = 0;
					reaction_number = ap * problem->number_of_reactions_per_nuc + 5 * problem->n_target_genes + i;
					propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = 0;
					reaction_number = ap * problem->number_of_reactions_per_nuc + 6 * problem->n_target_genes + i;
					propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = 0;
					reaction_number = ap * problem->number_of_reactions_per_nuc + 7 * problem->n_target_genes + i;
					propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = 0;
					reaction_number = ap * problem->number_of_reactions_per_nuc + 8 * problem->n_target_genes + i;
					propensity[rep * problem->number_of_reactions_per_nuc * tc->n_nucs + reaction_number] = 0;
				}
			}
	/*		iter_kounter++;*/
			t_slow += tau_slow;
		}
	}
	g_free(propensity);
	g_free(seedrn);
	if (interphase == 1) { /* interphase */
		g_free(propensity_fast);
		g_free(site_tab);
	}
	time2 = g_get_real_time();
	double elapsed = (time2 - time1)/G_USEC_PER_SEC;
	printf("tc %d: %f s\n", tc->kounter, elapsed);
	if (verbose) mssa_print_timeclass (tc, problem);
}

/* External inputs
 * are added
 */

void inject (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa inject %d\n", tc->kounter);
#pragma omp parallel for collapse(3) schedule(static) default(none) shared(problem, tc, mol_per_conc)
	for (int l = 0; l < problem->repeats; l++) {
		for (int i = 0; i < tc->n_nucs; i++) {
			for (int j = 0; j < problem->n_external_genes; j++) {
				int k = problem->external_gene_index[j];
				int conc = (int)(tc->data_protein[i * problem->n_tfs + k] * mol_per_conc) - tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k];
				tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k] = MAX(conc, 0);
			}
		}
	}
	if (verbose) mssa_print_timeclass (tc, problem);
}

void zero_structure (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa zero_structure %d\n", tc->kounter);
#pragma omp parallel for collapse(3) schedule(static) default(none) shared(problem, tc)
	for (int l = 0; l < problem->repeats; l++) {
		for (int i = 0; i < tc->n_nucs; i++) {
			for (int j = 0; j < problem->n_tfs; j++) {
				tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j] = 0;
				tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j] = 0;
				tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j] = 0;
			}
		}
	}
	if (verbose) mssa_print_timeclass (tc, problem);
}

void unbound (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa unbound %d\n", tc->kounter);
#pragma omp parallel for collapse(3) schedule(static) default(none) shared(problem, tc)
	for (int l = 0; l < problem->repeats; l++) {
		for (int i = 0; i < tc->n_nucs; i++) {
			for (int j = 0; j < problem->n_tfs; j++) {
				tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j] += tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j];
				tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j] = 0;
			}
		}
	}
#pragma omp parallel for collapse(3) schedule(static) default(none) shared(problem)
	for (int l = 0; l < problem->repeats; l++) {
		for(int i = 0; i < problem->n_nucs; i++) {
			for (int j = 0; j < problem->sum_sites; j++) {
				problem->status_allele_0[l * problem->n_nucs * problem->sum_sites + i * problem->sum_sites + j] = 0;
				problem->status_allele_1[l * problem->n_nucs * problem->sum_sites + i * problem->sum_sites + j] = 0;
			}
		}
	}
	if (verbose) mssa_print_timeclass (tc, problem);
}

/* Nucleaar division
 * All protein is assumed to be already unbound!
 */

void divide (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa divide %d\n", tc->kounter);
	MSSA_Timeclass *tc_prev = (MSSA_Timeclass *) g_list_nth_data (problem->tc_list, (tc->kounter - 1));
	if (!tc_prev) return;
	if (tc->n_nucs > tc_prev->n_nucs * 2) {
#pragma omp parallel for schedule(static) default(none) shared(problem, tc, tc_prev)
		for (int l = 0; l < problem->repeats; l++) {
			int i = 0;
			for (int j = 0; j < problem->n_tfs; j++) {
				tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + 2 * i * problem->n_tfs + j] = tc_prev->solution_mrna[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j] / 2;
				tc->solution_protein[l * tc->n_nucs * problem->n_tfs + 2 * i * problem->n_tfs + j] = tc_prev->solution_protein[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j] / 2;
			}
			for (int i = 1; i < tc_prev->n_nucs - 1; i++) {
				for (int j = 0; j < problem->n_tfs; j++) {
					tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + 2 * i * problem->n_tfs + j] = tc_prev->solution_mrna[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j] / 2;
					tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + (2 * i + 1) * problem->n_tfs + j] = tc_prev->solution_mrna[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j] / 2;
					tc->solution_protein[l * tc->n_nucs * problem->n_tfs + 2 * i * problem->n_tfs + j] = tc_prev->solution_protein[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j] / 2;
					tc->solution_protein[l * tc->n_nucs * problem->n_tfs + (2 * i + 1) * problem->n_tfs + j] = tc_prev->solution_protein[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j] / 2;
				}
			}
			i = tc_prev->n_nucs - 1;
			for (int j = 0; j < problem->n_tfs; j++) {
				tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + (2 * i + 1) * problem->n_tfs + j] = tc_prev->solution_mrna[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j] / 2;
				tc->solution_protein[l * tc->n_nucs * problem->n_tfs + (2 * i + 1) * problem->n_tfs + j] = tc_prev->solution_protein[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j] / 2;
			}
		}
	} else {
#pragma omp parallel for schedule(static) default(none) shared(problem, tc, tc_prev)
		for (int l = 0; l < problem->repeats; l++) {
			for (int i = 0; i < tc_prev->n_nucs; i++) {
				for (int j = 0; j < problem->n_tfs; j++) {
					tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + 2 * i * problem->n_tfs + j] = tc_prev->solution_mrna[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j] / 2;
					tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + (2 * i + 1) * problem->n_tfs + j] = tc_prev->solution_mrna[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j] / 2;
					tc->solution_protein[l * tc->n_nucs * problem->n_tfs + 2 * i * problem->n_tfs + j] = tc_prev->solution_protein[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j] / 2;
					tc->solution_protein[l * tc->n_nucs * problem->n_tfs + (2 * i + 1) * problem->n_tfs + j] = tc_prev->solution_protein[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j] / 2;
				}
			}
		}
	}
	if (verbose) mssa_print_timeclass (tc, problem);
}

void konnect (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa konnect %d\n", tc->kounter);
	MSSA_Timeclass *tc_prev = (MSSA_Timeclass *) g_list_nth_data (problem->tc_list, (tc->kounter - 1));
	if (!tc_prev) return;
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(problem, tc, tc_prev)
	for (int l = 0; l < problem->repeats; l++) {
		for (int i = 0; i < tc_prev->n_nucs; i++) {
			for (int j = 0; j < problem->n_tfs; j++) {
				tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j] = tc_prev->solution_mrna[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j];
				tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j] = tc_prev->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j];
				tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j] = tc_prev->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j];
			}
		}
	}
	if (verbose) mssa_print_timeclass (tc, problem);
}

void mssa_mark_site_overlap (MSSA_Problem *problem, int range)
{
	int k = 0;
#pragma omp parallel for schedule(static) default(none) shared(range, problem, k)

		for (int i = 0; i < problem->n_target_genes; i++) {
			for (int j = 0; j < problem->n_sites[i]; j++) {
				int kount_fw = 0;
				int kount_bk = 0;
				for (int l = j + 1; l < problem->n_sites[i]; l++) {
					if (problem->allele_0[k][i][l].coordinate > problem->allele_0[k][i][j].coordinate + range)
						break;
					kount_fw++;
				}
				for (int l = j - 1; l > -1; l--) {
					if (problem->allele_0[k][i][l].coordinate < problem->allele_0[k][i][j].coordinate - range)
						break;
					kount_bk++;
				}
				problem->allele_0[k][i][j].repressed_fw = problem->allele_1[k][i][j].repressed_fw = kount_fw;
				problem->allele_0[k][i][j].repressed_bk = problem->allele_1[k][i][j].repressed_bk = kount_bk;
				kount_fw = 0;
				kount_bk = 0;
				for (int l = j + 1; l < problem->n_sites[i]; l++) {
					if (problem->allele_0[k][i][l].coordinate > problem->allele_0[k][i][j].coordinate + problem->allele_0[k][i][j].length)
						break;
					kount_fw++;
				}
				for (int l = j - 1; l > -1; l--) {
					if (problem->allele_0[k][i][l].coordinate < problem->allele_0[k][i][j].coordinate - problem->allele_0[k][i][j].length)
						break;
					kount_bk++;
				}
				problem->allele_0[k][i][j].blocked_fw = problem->allele_1[k][i][j].blocked_fw = kount_fw;
				problem->allele_0[k][i][j].blocked_bk = problem->allele_1[k][i][j].blocked_bk = kount_bk;
			}
		}
}

void run_parallel(MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	switch(parallel) {
		case 0:
			propagate_with_transport_4 (tc, problem);
		break;
		case 2:
#ifdef OPENCL
			propagate_with_transport_5 (tc, problem);
#else
			g_error("This type %d (OpenCL) is not compiled", parallel);
#endif
		break;
		case 3:
#ifdef DOCL
			propagate_with_transport_6 (tc, problem);
#else
			g_error("This type %d (OpenCL kernel as plain C) is not compiled", parallel);
#endif
		break;
		default:
			g_error("This type %d is not compiled", parallel);
		break;

	}
}

void integrate (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
/* Debug tag */
	if (verbose) printf("multiscale_ssa %d integrate %d\n", problem->repeat, tc->kounter);
/* Print initial condition */
	if (tc->kounter == 0 && verbose) mssa_print_timeclass (tc, problem);
/* Copy result from previous TC */
	if (tc->type > 0) konnect (tc, problem);
/* Add bias or set the initial cond */
	if (tc->type == 2) {
		inject (tc, problem);
		add_bias (tc, problem);
		if (!dryrun) {
			run_parallel (tc, problem);
		}
	}
/* Run the model */
	if (tc->type == 1) {
		if (tc->has_data == 1) {
			inject (tc, problem);
			score (tc, problem);
		}
		if (!dryrun) {
			run_parallel (tc, problem);
		}
	}
/* Nuclear division, All protein is to be unbound already */
	if (tc->type == 0) divide (tc, problem);
/* Run the model in the mitosis mode - without fast reactions,
 * translation and transcription.
 * Firstly, unbound all proteins!
 */
	if (tc->type == 3) {
		unbound (tc, problem);
		if (!dryrun) {
			run_parallel (tc, problem);
		}
	}
	mssa_out_timeclass (tc, problem);
}

/*
 * Standard gettext macros.
*/
#ifdef ENABLE_NLS
#  include <libintl.h>
#  undef _
#  define _(String) dgettext (PACKAGE, String)
#  ifdef gettext_noop
#    define N_(String) gettext_noop (String)
#  else
#    define N_(String) (String)
#  endif
#else
#  define textdomain(String) (String)
#  define gettext(String) (String)
#  define dgettext(Domain,Message) (Message)
#  define dcgettext(Domain,Message,Type) (Message)
#  define bindtextdomain(Domain,Directory) (Domain)
#  define _(String) (String)
#  define N_(String) (String)
#endif

static gboolean
option_version_cb (const gchar *option_name,
                   const gchar *value,
                   gpointer     data,
                   GError     **error)
{
  g_print ("%s %s\n", _("Multiscale SSA, version "), VERSION);

  exit (0);
  return FALSE;
}

static GOptionEntry entries[] =
{
	{ "dryrun", 'd', 0, G_OPTION_ARG_NONE, &dryrun, N_("Calculations are not performed"), N_("DRYRUN") },
	{ "logfile", 'l', 0, G_OPTION_ARG_STRING, &log_file, N_("File name for progress"), N_("FILENAME") },
	{ "outfile", 'o', 0, G_OPTION_ARG_STRING, &out_file, N_("File name for concentrations"), N_("FILENAME") },
	{ "action", 'a', 0, G_OPTION_ARG_STRING, &action, N_("What to do"), N_("OPERATION") },
	{ "gpunumber", 'g', 0, G_OPTION_ARG_INT, &gpu_number, N_("Selected gpu"), N_("GPU NUMBER") },
	{ "parallel", 'p', 0, G_OPTION_ARG_INT, &parallel, N_("Type of parallelilism"), N_("PARALLEL") },
	{ "repeat", 'r', 0, G_OPTION_ARG_INT, &repeat, N_("Number of repeats"), N_("REPEAT") },
	{ "verbose", 'e', 0, G_OPTION_ARG_NONE, &scores, N_("scores"), N_("SCORES") },
	{ "molperconc", 'm', 0, G_OPTION_ARG_INT, &mol_per_conc, N_("Number of molecules per concentration unit"), N_("MOL") },
	{ "fasttimemax", 't', 0, G_OPTION_ARG_DOUBLE, &fast_time_max, N_("Fast time limit"), N_("FASTTIME") },
	{ "faststepsfrac", 's', 0, G_OPTION_ARG_DOUBLE, &fast_steps_frac, N_("Fraction of fast steps"), N_("FRAC") },
	{ "verbose", 'v', 0, G_OPTION_ARG_NONE, &verbose, N_("verbose"), N_("VERBOSE") },
	{ "version", 'V', G_OPTION_FLAG_NO_ARG | G_OPTION_FLAG_HIDDEN, G_OPTION_ARG_CALLBACK, option_version_cb, NULL, NULL },
	{ NULL }
};

int main(int argc, char**argv)
{
	if (verbose) printf("multiscale_ssa start\n");
	GOptionContext *context;
	GError *gerror = NULL;
	context = g_option_context_new (_("- MSSA blastoderm"));
	g_option_context_add_main_entries(context, (const GOptionEntry *)entries, NULL);
	g_option_context_set_ignore_unknown_options(context, TRUE);
	if (!g_option_context_parse (context, &argc, &argv, &gerror)) {
		g_error (_("option parsing failed: %s\n"), gerror->message);
	}
	g_option_context_free (context);
	if (argc < 2) {
		g_error(_("%s called without data file"), g_get_prgname());
	}
	data_file = g_strdup(argv[1]);
	if (log_file == NULL) {
		g_warning(_("%s called without log file"), g_get_prgname());
	}
	if (out_file == NULL) {
		g_warning(_("%s called without out file"), g_get_prgname());
	}
	MSSA_Problem *problem = mssa_read_problem(data_file);
	mssa_mark_site_overlap (problem, problem->parameters->range);
	grand = g_rand_new ();
	if (verbose) printf("multiscale_ssa read problem\n");
	if (verbose) printf("multiscale_ssa nnucs %d\n", problem->n_nucs);
	if (verbose) printf("multiscale_ssa tfs %d\n", problem->n_tfs);
	if (verbose) printf("multiscale_ssa targets %d\n", problem->n_target_genes);
	setup_problem(problem);
#ifdef OPENCL
	if (parallel == 1 || parallel == 2) setup_device(problem);
#endif
	if (!g_strcmp0 (action, "batch")) {
		g_warning(_("%s called for operation batch"), g_get_prgname());
		char tag[10];
		FILE *fp = fopen(argv[argc - 1], "r");
		fscanf(fp, "%s", tag);
		for (int i = 0; i < problem->n_target_genes; i++) {
			for (int j = 0; j < problem->n_tfs; j++) {
				fscanf(fp, "%lf", &(problem->parameters->T[i * problem->n_tfs + j]));
			}
			fscanf(fp, "%lf", &(problem->parameters->translation[i]));
			fscanf(fp, "%lf", &(problem->parameters->protein_degradation[i]));
			fscanf(fp, "%lf", &(problem->parameters->mrna_degradation[i]));
			fscanf(fp, "%lf", &(problem->parameters->transport_mrna[i]));
			fscanf(fp, "%lf", &(problem->parameters->transport_protein[i]));
		}
		fclose (fp);
		g_list_foreach (problem->tc_list, (GFunc) integrate, (gpointer) problem);
		g_list_foreach (problem->tc_list, (GFunc) print_scores, (gpointer) problem);
	} else if (!g_strcmp0 (action, "repl")) {
		char tag[10];
/*		while(1) {
		scanf("%s", tag);
		fprintf(stderr, "%s\n", tag);
		printf("\r\nfflush\r\n");
		fprintf(stderr, "\r\n%s\n", tag);
			fflush (stderr);
			fflush (stdout);
		}
		scanf("%*s");
		scanf("%*s");
		scanf("%*s");
		scanf("%*s");
		printf("fflush\r\n");
		printf("fflush\r\n");*/
		g_warning(_("%s called for operation optimize"), g_get_prgname());
//		printf("\r\nfflush\r\n");
//		fflush (stdout);
		while (!feof(stdin)) {
			scanf("%s(", tag);
			for (int i = 0; i < problem->n_target_genes; i++) {
				for (int j = 0; j < problem->n_tfs; j++) {
					scanf("%lf,", &(problem->parameters->T[i * problem->n_tfs + j]));
				}
				scanf("%lf,", &(problem->parameters->translation[i]));
				scanf("%lf,", &(problem->parameters->protein_degradation[i]));
				scanf("%lf,", &(problem->parameters->mrna_degradation[i]));
				scanf("%lf,", &(problem->parameters->transport_mrna[i]));
				if (i < problem->n_target_genes - 1) {
					scanf("%lf,", &(problem->parameters->transport_protein[i]));
				} else {
					scanf("%lf)", &(problem->parameters->transport_protein[i]));
				}
				scanf("\r\n%s\r\n", tag);
			}
			g_list_foreach (problem->tc_list, (GFunc) integrate, (gpointer) problem);
			g_list_foreach (problem->tc_list, (GFunc) print_scores, (gpointer) problem);
			printf("\r\nfflush\r\n");
			fflush (stdout);
		}
	} else {
		g_list_foreach (problem->tc_list, (GFunc) integrate, (gpointer) problem);
		g_list_foreach (problem->tc_list, (GFunc) print_scores, (gpointer) problem);
	}
#ifdef OPENCL
	if (parallel == 1 || parallel == 2) clear_device(problem);
#endif
	return (0);
}
