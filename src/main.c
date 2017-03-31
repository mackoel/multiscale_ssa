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
	double *solution_protein;
	double *bound_protein;
	double *data_mrna; /* a table of concentrations/molecule numbers */
	double *solution_mrna;
	double *corr;
	double *chisq;
	int *mrna_production_kount;
	int *protein_production_kount;
	int *reaction_kount;
	double *mrna_production_time;
	double *protein_production_time;
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
	double *translation_delay;
	double *transcription_delay;
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
static int use_delays = 0;

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
		tc->solution_protein = g_new0(double, tc->n_nucs * problem->n_tfs * problem->repeats);
		tc->bound_protein = g_new0(double, tc->n_nucs * problem->n_tfs * problem->repeats);
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
		tc->solution_mrna = g_new0(double, tc->n_nucs * problem->n_tfs * problem->repeats);
		tc->mrna_production_kount = g_new0(int, tc->n_nucs * problem->n_target_genes * problem->repeats);
		tc->mrna_production_time = g_new0(double, tc->n_nucs * problem->n_target_genes * problem->repeats);
		tc->protein_production_kount = g_new0(int, tc->n_nucs * problem->n_target_genes * problem->repeats);
		tc->protein_production_time = g_new0(double, tc->n_nucs * problem->n_target_genes * problem->repeats);
		tc->reaction_kount = g_new0(int, tc->n_nucs * problem->repeats);
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
		if (use_delays == 1) {
			fscanf(fp, "%lf", &(problem->parameters->transcription_delay[i]));
			fscanf(fp, "%lf", &(problem->parameters->translation_delay[i]));
		}
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
		int m = problem->m_sites[i];
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

	return(problem);
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

void mssa_check_site_overlap (MSSA_site ***allele, int nuc, int gen, int cur, int type, int block)
{
	int kount_fw = (type == 0) ? allele[nuc][gen][cur].blocked_fw : allele[nuc][gen][cur].repressed_fw;
	int kount_bk = (type == 0) ? allele[nuc][gen][cur].blocked_bk : allele[nuc][gen][cur].repressed_bk;
/*	kount_fw = MIN(kount_fw, 1);
	kount_bk = MIN(kount_bk, 1);*/
	for (int l = 0; l < kount_fw; l++) {
		if (allele[nuc][gen][cur + 1 + l].status != 1)
			allele[nuc][gen][cur + 1 + l].status = block;
	}
	for (int l = 0; l < kount_bk; l++) {
		if (allele[nuc][gen][cur - 1 - l].status != 1)
			allele[nuc][gen][cur - 1 - l].status = block;
	}
}

void mssa_print_timeclass (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	for (int l = 0; l < problem->repeats; l++) {
		fprintf(stdout, "%d time %f %f protein:\n", problem->repeat, tc->t_start, tc->t_end);
		for (int i = 0; i < tc->n_nucs; i++) {
			for (int j = 0; j < problem->n_tfs; j++) {
				fprintf(stdout, "%f ", tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
			}
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "time %f %f bound:\n", tc->t_start, tc->t_end);
		for (int i = 0; i < tc->n_nucs; i++) {
			for (int j = 0; j < problem->n_tfs; j++) {
				fprintf(stdout, "%f ", tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
			}
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "time %f %f mrna:\n", tc->t_start, tc->t_end);
		for (int i = 0; i < tc->n_nucs; i++) {
			for (int j = 0; j < problem->n_tfs; j++) {
				fprintf(stdout, "%f ", tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
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
				fprintf(fp, " %.0f", tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k] + tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k]);
			}
			for (int j = 0; j < problem->n_target_genes; j++) {
				int k = problem->target_gene_index[j];
				fprintf(fp, " %.0f", tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k]);
			}
			for (int j = 0; j < problem->n_target_genes; j++) {
				int k = problem->target_gene_index[j];
				fprintf(fp, " %.0f", tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k]);
			}
			for (int j = 0; j < problem->n_target_genes; j++) {
				fprintf(fp, " %d", tc->mrna_production_kount[l * tc->n_nucs * problem->n_target_genes + i * problem->n_target_genes + j]);
			}
			for (int j = 0; j < problem->n_target_genes; j++) {
				fprintf(fp, " %.5f", tc->mrna_production_time[l * tc->n_nucs * problem->n_target_genes + i * problem->n_target_genes + j]);
			}
			for (int j = 0; j < problem->n_target_genes; j++) {
				fprintf(fp, " %d", tc->protein_production_kount[l * tc->n_nucs * problem->n_target_genes + i * problem->n_target_genes + j]);
			}
			for (int j = 0; j < problem->n_target_genes; j++) {
				fprintf(fp, " %.5f", tc->protein_production_time[l * tc->n_nucs * problem->n_target_genes + i * problem->n_target_genes + j]);
			}
			fprintf(fp, " %d", tc->reaction_kount[l * tc->n_nucs + i]);
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
				tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k] += ceil(tc->data_mrna[i * problem->n_tfs + k]);
				tc->solution_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k] += ceil(tc->data_protein[i * problem->n_tfs + k] * mol_per_conc);
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
	GRand* *seedrn = g_new0(GRand*, tc->n_nucs * problem->repeats);

#pragma omp parallel for collapse(2) schedule(static) default(none) shared(problem, tc, seedrn, seedrng)
	for (int rep = 0; rep < problem->repeats; rep++) {
		for (int ap = 0; ap < tc->n_nucs; ap++) {
			seedrn[rep * tc->n_nucs + ap] = g_rand_new_with_seed(seedrng + rep * tc->n_nucs + ap);
		}
	}
#ifndef NO_OPENMP_NESTED
	omp_set_nested(1);
/* OMP_NESTED=1*/
#pragma omp parallel for schedule(static) default(none) shared(problem, tc, propensity, propensity_fast, site_tab, t_start_slow, t_stop_slow, interphase, seedrn, fast_time_max)
#endif
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
							s = g_rand_int_range(seedrn[rep * tc->n_nucs + ap], 0, problem->n_sites[i] - 1);
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
						found = -1;
						aggregate = 0;
						random = prop_sum * g_rand_double(seedrn[rep * tc->n_nucs + ap]);
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
						tau_fast = -log(g_rand_double(seedrn[rep * tc->n_nucs + ap])) / prop_sum;
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
									if (tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] < -1e-8) {
	 									g_warning("fast bind reaction n %d t %d: ap %d tf %d target %d: %f < 0", allele, reaction_type, ap, tf, promoter_number, tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf]);
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
									if (tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] < -1e-8) {
	 									g_warning("fast unbind reaction n %d t %d: ap %d tf %d target %d: %f < 0", found, reaction_type, ap, tf, promoter_number, tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf]);
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
									if (tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] < -1e-8) {
	 									g_warning("fast bind reaction n %d t %d: ap %d tf %d target %d: %f < 0", allele, reaction_type, ap, tf, promoter_number, tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf]);
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
									if (tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf] < -1e-8) {
										g_warning("fast unbind reaction n %d t %d: ap %d tf %d target %d: %f < 0", found, reaction_type, ap, tf, promoter_number, tc->bound_protein[rep * tc->n_nucs * problem->n_tfs + ap * problem->n_tfs + tf]);
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
			found = 0;
			aggregate = 0;
			random = prop_sum_all * g_rand_double(seedrn[rep * tc->n_nucs]);
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
			tau_slow = -log(g_rand_double(seedrn[rep * tc->n_nucs + reaction_nuc])) / prop_sum_all;
			tc->reaction_kount[rep * tc->n_nucs + reaction_nuc]++;
			switch (reaction_index) {
				case 1: /* transcription 1 */
					k = problem->target_gene_index[reaction_target];
					tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k]++;
					tc->mrna_production_kount[rep * tc->n_nucs * problem->n_target_genes + reaction_nuc * problem->n_target_genes + reaction_target]++;
					tc->mrna_production_time[rep * tc->n_nucs * problem->n_target_genes + reaction_nuc * problem->n_target_genes + reaction_target] += tau_slow;
				break;
				case 2: /* transcription 2 */
					k = problem->target_gene_index[reaction_target];
					tc->solution_mrna[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k]++;
					tc->mrna_production_kount[rep * tc->n_nucs * problem->n_target_genes + reaction_nuc * problem->n_target_genes + reaction_target]++;
					tc->mrna_production_time[rep * tc->n_nucs * problem->n_target_genes + reaction_nuc * problem->n_target_genes + reaction_target] += tau_slow;
				break;
				case 3: /* translation */
					k = problem->target_gene_index[reaction_target];
					tc->solution_protein[rep * tc->n_nucs * problem->n_tfs + reaction_nuc * problem->n_tfs + k]++;
					tc->protein_production_kount[rep * tc->n_nucs * problem->n_target_genes + reaction_nuc * problem->n_target_genes + reaction_target]++;
					tc->protein_production_time[rep * tc->n_nucs * problem->n_target_genes + reaction_nuc * problem->n_target_genes + reaction_target] += tau_slow;
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
				double conc = ceil(tc->data_protein[i * problem->n_tfs + k] * mol_per_conc) - tc->bound_protein[l * tc->n_nucs * problem->n_tfs + i * problem->n_tfs + k];
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
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(problem, tc)
	for(int k = 0; k < problem->n_nucs; k++) {
		for (int i = 0; i < problem->n_target_genes; i++) {
			for (int j = 0; j < problem->n_sites[i]; j++) {
				problem->allele_0[k][i][j].status = 0;
				problem->allele_1[k][i][j].status = 0;
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
				tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + 2 * i * problem->n_tfs + j] = ceil(0.5 * tc_prev->solution_mrna[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
				tc->solution_protein[l * tc->n_nucs * problem->n_tfs + 2 * i * problem->n_tfs + j] = ceil(0.5 * tc_prev->solution_protein[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
			}
			for (int i = 1; i < tc_prev->n_nucs - 1; i++) {
				for (int j = 0; j < problem->n_tfs; j++) {
					tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + 2 * i * problem->n_tfs + j] = ceil(0.5 * tc_prev->solution_mrna[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
					tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + (2 * i + 1) * problem->n_tfs + j] = floor(0.5 * tc_prev->solution_mrna[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
					tc->solution_protein[l * tc->n_nucs * problem->n_tfs + 2 * i * problem->n_tfs + j] = ceil(0.5 * tc_prev->solution_protein[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
					tc->solution_protein[l * tc->n_nucs * problem->n_tfs + (2 * i + 1) * problem->n_tfs + j] = floor(0.5 * tc_prev->solution_protein[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
				}
			}
			i = tc_prev->n_nucs - 1;
			for (int j = 0; j < problem->n_tfs; j++) {
				tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + (2 * i + 1) * problem->n_tfs + j] = ceil(0.5 * tc_prev->solution_mrna[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
				tc->solution_protein[l * tc->n_nucs * problem->n_tfs + (2 * i + 1) * problem->n_tfs + j] = ceil(0.5 * tc_prev->solution_protein[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
			}
		}
	} else {
#pragma omp parallel for schedule(static) default(none) shared(problem, tc, tc_prev)
		for (int l = 0; l < problem->repeats; l++) {
			for (int i = 0; i < tc_prev->n_nucs; i++) {
				for (int j = 0; j < problem->n_tfs; j++) {
					tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + 2 * i * problem->n_tfs + j] = ceil(0.5 * tc_prev->solution_mrna[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
					tc->solution_mrna[l * tc->n_nucs * problem->n_tfs + (2 * i + 1) * problem->n_tfs + j] = floor(0.5 * tc_prev->solution_mrna[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
					tc->solution_protein[l * tc->n_nucs * problem->n_tfs + 2 * i * problem->n_tfs + j] = ceil(0.5 * tc_prev->solution_protein[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
					tc->solution_protein[l * tc->n_nucs * problem->n_tfs + (2 * i + 1) * problem->n_tfs + j] = floor(0.5 * tc_prev->solution_protein[l * tc_prev->n_nucs * problem->n_tfs + i * problem->n_tfs + j]);
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
/*
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(range, problem)
	for(int k = 0; k < problem->n_nucs; k++) {
*/
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
/*
	}
*/
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
			propagate_with_transport_4 (tc, problem);
		}
	}
/* Run the model */
	if (tc->type == 1) {
		if (tc->has_data == 1) {
			inject (tc, problem);
			score (tc, problem);
		}
		if (!dryrun) {
			propagate_with_transport_4 (tc, problem);
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
			propagate_with_transport_4 (tc, problem);
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
	{ "usedelays", 'u', 0, G_OPTION_ARG_INT, &use_delays, N_("Use delays"), N_("DELAYS") },
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
	return (0);
}
