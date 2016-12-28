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
	int sum_sites; /* total number of sites */
	MSSA_site ***allele_0; /* sites on the first allele */
	MSSA_site ***allele_1; /* sites on the second allele */
	MSSA_Parameters *parameters;
	int repeat;
} MSSA_Problem;

static char*data_file;
static char*out_file;
static char*log_file;
static char*action;
static gboolean verbose = FALSE;
static gboolean dryrun = FALSE;
static int repeat;
static int mol_per_conc = MOLECULES_PER_CONCENTRATION;
static double fast_time_max = FAST_TIME_MAX;
static double fast_steps_frac = FAST_STEPS_FRAC;

MSSA_Problem *mssa_read_problem(gchar*filename)
{
	MSSA_Problem *problem = g_new(MSSA_Problem, 1);
//	char y[10];
	FILE*fp = g_fopen(filename, "r");
	fscanf(fp, "%*s\n");
//		fscanf(fp, "%s\n", y);
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
		tc->solution_protein = g_new0(double, tc->n_nucs * problem->n_tfs);
		tc->bound_protein = g_new0(double, tc->n_nucs * problem->n_tfs);
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
		tc->solution_mrna = g_new0(double, tc->n_nucs * problem->n_tfs);
		problem->tc_list = g_list_append(problem->tc_list, tc);
	}
	problem->sum_sites = 0;
	problem->n_sites = g_new0(int, problem->n_target_genes);
	for (int i = 0; i < problem->n_target_genes; i++) {
		fscanf(fp, "%d", &(problem->n_sites[i]));
		problem->sum_sites += problem->n_sites[i];
	}
	fscanf(fp, "%*s");
	problem->allele_0 = g_new0(MSSA_site**, problem->n_nucs);
	problem->allele_1 = g_new0(MSSA_site**, problem->n_nucs);
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
	for(int k = 1; k < problem->n_nucs; k++) {
		problem->allele_0[k] = g_new0(MSSA_site*, problem->n_target_genes);
		problem->allele_1[k] = g_new0(MSSA_site*, problem->n_target_genes);
		for (int i = 0; i < problem->n_target_genes; i++) {
			problem->allele_0[k][i] = g_new0(MSSA_site, problem->n_sites[i]);
			problem->allele_1[k][i] = g_new0(MSSA_site, problem->n_sites[i]);
			for (int j = 0; j < problem->n_sites[i]; j++) {
				problem->allele_0[k][i][j] = problem->allele_0[0][i][j];
				problem->allele_1[k][i][j] = problem->allele_1[0][i][j];
			}
		}
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

int mssa_get_reaction_fast(double *solution,
                           double *bound,
                           MSSA_site ***allele,
                           double *tau,
                           int *reaction_type,
                           int *tf,
                           int *target,
                           int *site_number,
                           int *n_sites,
                           int n_tfs,
                           int n_target_genes,
                           int ap)
{
	int i, j, k, reaction_number, s, l;
	double *propensity;
	double *probability;
	double aggregate;
	double prop_sum = 0;
	double random = g_rand_double (grand);
//	fprintf(stdout, "fast %d!\n", ap);
	int number_of_reactions = 2 * n_tfs * n_target_genes; /* bind/unbind each tf in each promotor */
	propensity = g_new0(double, number_of_reactions);
	probability = g_new0(double, number_of_reactions);
/* Binding */
	for (i = 0; i < n_target_genes; i++) {
		for (j = 0; j < n_tfs; j++) {
			double prop = 0; /* Sum of energies of free bs */
			for (k = 0; k < n_sites[i]; k++) {
				prop += (allele[ap][i][k].status == 0 && allele[ap][i][k].tf_index == j) ? allele[ap][i][k].energy : 0;
			}
			propensity[i * n_tfs + j] = solution[ap * n_tfs + j] * prop;
			prop_sum += propensity[i * n_tfs + j];
		}

	}
/* Unbinding */
	for (i = 0; i < n_target_genes; i++) {
		for (j = 0; j < n_tfs; j++) {
			double prop = 0; /* Sum of energies of bound bs */
			for (k = 0; k < n_sites[i]; k++) {
				prop += (allele[ap][i][k].status == 1 && allele[ap][i][k].tf_index == j) ? allele[ap][i][k].energy : 0;
			}
//			propensity[n_tfs * n_target_genes + i * n_tfs + j] = solution[ap * n_tfs + j] * prop;
			propensity[n_tfs * n_target_genes + i * n_tfs + j] = bound[ap * n_tfs + j] * prop;
			prop_sum += propensity[n_tfs * n_target_genes + i * n_tfs + j];
		}

	}
	if (prop_sum <= 0.00000000001) {
		g_warning("fast %d: Sum of propensities is too small %g!", ap, prop_sum);
		(*tau) = TMAX;
		g_free(propensity);
		g_free(probability);
		return(-1);
	}
	reaction_number = -1;
	aggregate = 0;
/* Binding */
	for (i = 0; i < n_target_genes; i++) {
		for (j = 0; j < n_tfs; j++) {
			aggregate += propensity[i * n_tfs + j];
			probability[i * n_tfs + j] = aggregate / prop_sum;
			if (random < probability[i * n_tfs + j]) {
				reaction_number = i * n_tfs + j;
				(*target) = i;
				(*tf) = j;
				(*reaction_type) = 1;
				break;
			}
		}
		if (reaction_number > 0) break;
	}
/* Unbinding */
	if (reaction_number < 0) {
		for (i = 0; i < n_target_genes; i++) {
			for (j = 0; j < n_tfs; j++) {
				aggregate += propensity[n_tfs * n_target_genes + i * n_tfs + j];
				probability[n_tfs * n_target_genes + i * n_tfs + j] = aggregate / prop_sum;
				if (random < probability[n_tfs * n_target_genes + i * n_tfs + j]) {
					reaction_number = n_tfs * n_target_genes + i * n_tfs + j;
					(*target) = i;
					(*tf) = j;
					(*reaction_type) = 0;
					break;
				}
			}
			if (reaction_number > 0) break;
		}
	}
	(*tau) = -log(g_rand_double (grand)) / prop_sum;
	(*site_number) = -1;
	l = floor(g_rand_double (grand) * n_sites[(*target)]);
	if (reaction_number > 0) {
		s = l;
		for (k = 0; k < n_sites[(*target)]; k++) {
			/* Looking for a bound/unbound site */
			if (allele[ap][(*target)][s].status == (1 - (*reaction_type))) {
				(*site_number) = s;
				allele[ap][(*target)][s].status = (*reaction_type);
				solution[ap * n_tfs + (*tf)] += ((*reaction_type) == 1) ? -1 : 1;
				bound[ap * n_tfs + (*tf)] -= ((*reaction_type) == 1) ? -1 : 1;
				if (solution[ap * n_tfs + (*tf)] < 0) {
					g_warning("fast reaction n %d t %d: ap %d tf %d target %d < 0", reaction_number, (*reaction_type), ap, (*tf), (*target));
					solution[ap * n_tfs + (*tf)] = 0;
				}
				break;
			}
			s++;
			if (s > n_sites[(*target)] - 1) {
				s = 0;
			}
		}
	}
	g_free(propensity);
	g_free(probability);
	return (reaction_number);
}

int mssa_get_reaction_fast_with_buffers(double *solution,
                                        double *bound,
                                        MSSA_site ***allele,
                                        double *propensity,
                                        double *probability,
                                        double *tau,
                                        int *reaction_type,
                                        int *tf,
                                        int *target,
                                        int *site_number,
                                        int *n_sites,
                                        int n_tfs,
                                        int n_target_genes,
                                        int ap)
{
	int i, j, k, reaction_number, s, l;
	double aggregate;
	double prop_sum = 0;
	double random = g_rand_double (grand);
/* Binding */
	for (i = 0; i < n_target_genes; i++) {
		for (j = 0; j < n_tfs; j++) {
			double prop = 0; /* Sum of energies of free bs */
			for (k = 0; k < n_sites[i]; k++) {
				prop += (allele[ap][i][k].status == 0 && allele[ap][i][k].tf_index == j) ? allele[ap][i][k].energy : 0;
			}
			propensity[i * n_tfs + j] = solution[ap * n_tfs + j] * prop;
			prop_sum += propensity[i * n_tfs + j];
		}

	}
/* Unbinding */
	for (i = 0; i < n_target_genes; i++) {
		for (j = 0; j < n_tfs; j++) {
			double prop = 0; /* Sum of energies of bound bs */
			for (k = 0; k < n_sites[i]; k++) {
				prop += (allele[ap][i][k].status == 1 && allele[ap][i][k].tf_index == j) ? allele[ap][i][k].energy : 0;
			}
//			propensity[n_tfs * n_target_genes + i * n_tfs + j] = solution[ap * n_tfs + j] * prop;
			propensity[n_tfs * n_target_genes + i * n_tfs + j] = bound[ap * n_tfs + j] * prop;
			prop_sum += propensity[n_tfs * n_target_genes + i * n_tfs + j];
		}

	}
	if (prop_sum <= 0.00000000001) {
		g_warning("fast %d: Sum of propensities is too small %g!", ap, prop_sum);
		(*tau) = TMAX;
		(*target) = -1;
		(*tf) = -1;
		(*reaction_type) = -1;
		(*site_number) = -1;
		return(-1);
	}
	reaction_number = -1;
	aggregate = 0;
/* Binding */
	for (i = 0; i < n_target_genes; i++) {
		for (j = 0; j < n_tfs; j++) {
			aggregate += propensity[i * n_tfs + j];
			probability[i * n_tfs + j] = aggregate / prop_sum;
			if (random < probability[i * n_tfs + j]) {
				reaction_number = i * n_tfs + j;
				(*target) = i;
				(*tf) = j;
				(*reaction_type) = 1;
				break;
			}
		}
		if (reaction_number > 0) break;
	}
/* Unbinding */
	if (reaction_number < 0) {
		for (i = 0; i < n_target_genes; i++) {
			for (j = 0; j < n_tfs; j++) {
				aggregate += propensity[n_tfs * n_target_genes + i * n_tfs + j];
				probability[n_tfs * n_target_genes + i * n_tfs + j] = aggregate / prop_sum;
				if (random < probability[n_tfs * n_target_genes + i * n_tfs + j]) {
					reaction_number = n_tfs * n_target_genes + i * n_tfs + j;
					(*target) = i;
					(*tf) = j;
					(*reaction_type) = 0;
					break;
				}
			}
			if (reaction_number > 0) break;
		}
	}
	(*tau) = -log(g_rand_double (grand)) / prop_sum;
	(*site_number) = -1;
	l = floor(g_rand_double (grand) * n_sites[(*target)]);
	if (reaction_number > 0) {
		s = l;
		for (k = 0; k < n_sites[(*target)]; k++) {
			/* Looking for a bound/unbound site */
			if (allele[ap][(*target)][s].status == (1 - (*reaction_type))) {
				(*site_number) = s;
				allele[ap][(*target)][s].status = (*reaction_type);
				solution[ap * n_tfs + (*tf)] += ((*reaction_type) == 1) ? -1 : 1;
				bound[ap * n_tfs + (*tf)] -= ((*reaction_type) == 1) ? -1 : 1;
				if (solution[ap * n_tfs + (*tf)] < 0) {
					g_warning("fast reaction n %d t %d: ap %d tf %d target %d < 0", reaction_number, (*reaction_type), ap, (*tf), (*target));
					solution[ap * n_tfs + (*tf)] = 0;
					(*tau) = TMAX;
				}
				break;
			}
			s++;
			if (s > n_sites[(*target)] - 1) {
				s = 0;
			}
		}
	}
	return (reaction_number);
}

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

int mssa_get_reaction_fast_with_buffers_2(double *solution,
                                          double *bound,
                                          MSSA_site ***allele_0,
                                          MSSA_site ***allele_1,
                                          double *propensity,
                                          double *probability,
                                          int **site_tab_bind_0,
                                          int **site_tab_bind_1,
                                          int **site_tab_unbind_0,
                                          int **site_tab_unbind_1,
                                          double *tau,
                                          int *reaction_type,
                                          int *tf,
                                          int *target,
                                          int *site_number,
                                          int *n_sites,
                                          int n_tfs,
                                          int n_target_genes,
                                          double *parameters,
                                          double range,
                                          int ap)
{
	int i, j, k, reaction_number, s, l, found, allele;
	double aggregate;
	double prop_sum = 0;
	double random = g_rand_double (grand);
/* Binding */
	for (i = 0; i < n_target_genes; i++) {
		for (j = 0; j < n_tfs; j++) {
			double prop = 0; /* Sum of energies of free bs */
			reaction_number = i * n_tfs + j;
			l = g_rand_int_range (grand, 0, n_sites[i]);
			s = l;
			for (k = 0; k < n_sites[i]; k++) {
				if (allele_0[ap][i][s].status == 0 && allele_0[ap][i][s].tf_index == j) {
					prop += allele_0[ap][i][s].energy;
					site_tab_bind_0[i][j] = s;
				}
				s++;
				if (s > n_sites[i] - 1) {
					s = 0;
				}
			}
			propensity[reaction_number] = solution[ap * n_tfs + j] * prop;
			prop_sum += propensity[reaction_number];
		}
	}
	for (i = 0; i < n_target_genes; i++) {
		for (j = 0; j < n_tfs; j++) {
			double prop = 0; /* Sum of energies of free bs */
			reaction_number = n_tfs * n_target_genes + i * n_tfs + j;
			l = g_rand_int_range (grand, 0, n_sites[i]);
			s = l;
			for (k = 0; k < n_sites[i]; k++) {
				prop += (allele_1[ap][i][s].status == 0 && allele_1[ap][i][s].tf_index == j) ? allele_1[ap][i][k].energy : 0;
				if (allele_1[ap][i][s].status == 0 && allele_1[ap][i][s].tf_index == j) {
					prop += allele_1[ap][i][s].energy;
					site_tab_bind_1[i][j] = s;
				}
				s++;
				if (s > n_sites[i] - 1) {
					s = 0;
				}
			}
			propensity[reaction_number] = solution[ap * n_tfs + j] * prop;
			prop_sum += propensity[reaction_number];
		}
	}
/* Unbinding */
	for (i = 0; i < n_target_genes; i++) {
		for (j = 0; j < n_tfs; j++) {
			double prop = 0; /* Sum of energies of bound bs */
			reaction_number = 2 * n_tfs * n_target_genes + i * n_tfs + j;
			l = g_rand_int_range (grand, 0, n_sites[i]);
			s = l;
			for (k = 0; k < n_sites[i]; k++) {
				prop += (allele_0[ap][i][s].status == 1 && allele_0[ap][i][s].tf_index == j) ? allele_0[ap][i][s].energy : 0;
				if (allele_0[ap][i][s].status == 1 && allele_0[ap][i][s].tf_index == j) {
					prop += allele_0[ap][i][s].energy;
					site_tab_unbind_0[i][j] = s;
				}
				s++;
				if (s > n_sites[i] - 1) {
					s = 0;
				}
			}
			propensity[reaction_number] = bound[ap * n_tfs + j] * prop;
			prop_sum += propensity[reaction_number];
		}
	}
	for (i = 0; i < n_target_genes; i++) {
		for (j = 0; j < n_tfs; j++) {
			double prop = 0; /* Sum of energies of bound bs */
			reaction_number = 3 * n_tfs * n_target_genes + i * n_tfs + j;
			l = g_rand_int_range (grand, 0, n_sites[i]);
			s = l;
			for (k = 0; k < n_sites[i]; k++) {
				prop += (allele_1[ap][i][s].status == 1 && allele_1[ap][i][s].tf_index == j) ? allele_1[ap][i][s].energy : 0;
				if (allele_1[ap][i][s].status == 1 && allele_1[ap][i][s].tf_index == j) {
					prop += allele_1[ap][i][s].energy;
					site_tab_unbind_1[i][j] = s;
				}
				s++;
				if (s > n_sites[i] - 1) {
					s = 0;
				}
			}
			propensity[reaction_number] = bound[ap * n_tfs + j] * prop;
			prop_sum += propensity[reaction_number];
		}

	}
	if (prop_sum <= 0.00000000001) {
		if (verbose) {
			int y = 0;
			for (i = 0; i < n_target_genes; i++) {
				for (j = 0; j < n_tfs; j++) {
					if (bound[ap * n_tfs + j] > 0) y++;
				}
			}
			int z = 0;
			for (i = 0; i < n_target_genes; i++) {
				for (j = 0; j < n_tfs; j++) {
					if (solution[ap * n_tfs + j] > 0) z++;
				}
			}
			int p = 0;
			for (i = 0; i < 4 * n_tfs * n_target_genes; i++) {
				if (propensity[i] < 0) p++;
			}
			g_warning("fast %d: Sum of propensities is too small %g! bound %d sol %d p %d", ap, prop_sum, y, z, p);
		} else {
			g_warning("fast %d: Sum of propensities is too small %g!", ap, prop_sum);
		}
		(*tau) = TMAX;
		(*target) = -1;
		(*tf) = -1;
		(*reaction_type) = -1;
		(*site_number) = -1;
		return(-1);
	}
	found = 0;
	aggregate = 0;
/* Binding */
	for (i = 0; i < n_target_genes; i++) {
		for (j = 0; j < n_tfs; j++) {
			reaction_number = i * n_tfs + j;
			aggregate += propensity[reaction_number];
			probability[reaction_number] = aggregate / prop_sum;
			if (random < probability[reaction_number]) {
				found = 1;
				(*target) = i;
				(*tf) = j;
				(*reaction_type) = 1;
				allele = 0;
				break;
			}
		}
		if (found == 1) break;
	}
	if (found == 0) {
		for (i = 0; i < n_target_genes; i++) {
			for (j = 0; j < n_tfs; j++) {
				reaction_number = n_tfs * n_target_genes + i * n_tfs + j;
				aggregate += propensity[reaction_number];
				probability[reaction_number] = aggregate / prop_sum;
				if (random < probability[reaction_number]) {
					found = 1;
					(*target) = i;
					(*tf) = j;
					(*reaction_type) = 1;
					allele = 1;
					break;
				}
			}
			if (found == 1) break;
		}
	}
/* Unbinding */
	if (found == 0) {
		for (i = 0; i < n_target_genes; i++) {
			for (j = 0; j < n_tfs; j++) {
				reaction_number = 2 * n_tfs * n_target_genes + i * n_tfs + j;
				aggregate += propensity[reaction_number];
				probability[reaction_number] = aggregate / prop_sum;
				if (random < probability[reaction_number]) {
					found = 1;
					(*target) = i;
					(*tf) = j;
					(*reaction_type) = 0;
					allele = 0;
					break;
				}
			}
			if (found == 1) break;
		}
	}
	if (found == 0) {
		for (i = 0; i < n_target_genes; i++) {
			for (j = 0; j < n_tfs; j++) {
				reaction_number = 3 * n_tfs * n_target_genes + i * n_tfs + j;
				aggregate += propensity[reaction_number];
				probability[reaction_number] = aggregate / prop_sum;
				if (random < probability[reaction_number]) {
					found = 1;
					(*target) = i;
					(*tf) = j;
					(*reaction_type) = 0;
					allele = 1;
					break;
				}
			}
			if (found == 1) break;
		}
	}
	(*tau) = -log(g_rand_double (grand)) / prop_sum;
	(*site_number) = -1;
	if (found == 1) {
		if (allele == 0) {
			if ((*reaction_type) == 1) { /* Binding */
				s = site_tab_bind_0[(*target)][(*tf)];
				(*site_number) = s;
				allele_0[ap][(*target)][s].status = 1;
				solution[ap * n_tfs + (*tf)] += -1;
				bound[ap * n_tfs + (*tf)] -= -1;
				if (solution[ap * n_tfs + (*tf)] < 0) {
					g_warning("fast bind reaction n %d t %d: ap %d tf %d target %d < 0", reaction_number, (*reaction_type), ap, (*tf), (*target));
					solution[ap * n_tfs + (*tf)] = 0;
				}
			/* allelle, nuc, gen, cur, type, block */
				if (parameters[(*target) * n_tfs + (*tf)] < 0 && range > 0) {
					mssa_check_site_overlap (allele_0, ap, (*target), s, 0, -1);
				} else {
					mssa_check_site_overlap (allele_0, ap, (*target), s, 1, -1);
				}
			} else { /* Unbinding */
				s = site_tab_unbind_0[(*target)][(*tf)];
				(*site_number) = s;
				allele_0[ap][(*target)][s].status = 0;
				solution[ap * n_tfs + (*tf)] += 1;
				bound[ap * n_tfs + (*tf)] -= 1;
				if (bound[ap * n_tfs + (*tf)] < 0) {
					g_warning("fast unbind reaction n %d t %d: ap %d tf %d target %d < 0", reaction_number, (*reaction_type), ap, (*tf), (*target));
					bound[ap * n_tfs + (*tf)] = 0;
				}
			/* allelle, nuc, gen, cur, type, block */
				if (parameters[(*target) * n_tfs + (*tf)] < 0 && range > 0) {
					mssa_check_site_overlap (allele_0, ap, (*target), s, 0, 0);
				} else {
					mssa_check_site_overlap (allele_0, ap, (*target), s, 1, 0);
				}
			}
		}
		if (allele == 1) {
			if ((*reaction_type) == 1) { /* Binding */
				s = site_tab_bind_1[(*target)][(*tf)];
				(*site_number) = s;
				allele_1[ap][(*target)][s].status = 1;
				solution[ap * n_tfs + (*tf)] += -1;
				bound[ap * n_tfs + (*tf)] -= -1;
				if (solution[ap * n_tfs + (*tf)] < 0) {
					g_warning("fast bind reaction n %d t %d: ap %d tf %d target %d < 0", reaction_number, (*reaction_type), ap, (*tf), (*target));
					solution[ap * n_tfs + (*tf)] = 0;
				}
			/* allelle, nuc, gen, cur, type, block */
				if (parameters[(*target) * n_tfs + (*tf)] < 0 && range > 0) {
					mssa_check_site_overlap (allele_1, ap, (*target), s, 0, -1);
				} else {
					mssa_check_site_overlap (allele_1, ap, (*target), s, 1, -1);
				}
			} else { /* Unbinding */
				s = site_tab_unbind_1[(*target)][(*tf)];
				(*site_number) = s;
				allele_0[ap][(*target)][s].status = 0;
				solution[ap * n_tfs + (*tf)] += 1;
				bound[ap * n_tfs + (*tf)] -= 1;
				if (bound[ap * n_tfs + (*tf)] < 0) {
					g_warning("fast unbind reaction n %d t %d: ap %d tf %d target %d < 0", reaction_number, (*reaction_type), ap, (*tf), (*target));
					bound[ap * n_tfs + (*tf)] = 0;
				}
			/* allelle, nuc, gen, cur, type, block */
				if (parameters[(*target) * n_tfs + (*tf)] < 0 && range > 0) {
					mssa_check_site_overlap (allele_1, ap, (*target), s, 0, 0);
				} else {
					mssa_check_site_overlap (allele_1, ap, (*target), s, 1, 0);
				}
			}
		}
	}
	return (reaction_number);
}

int mssa_get_reaction_slow(double *solution_mrna,
                           double *solution_protein,
                           MSSA_site ***allele_0,
                           MSSA_site ***allele_1,
                           double *tau,
                           int *reaction_type,
                           int *target,
                           int *target_gene_index,
                           int *n_sites,
                           double *T,
                           double *protein_degradation,
                           double *mrna_degradation,
                           double *translation,
                           int n_tfs,
                           int n_target_genes,
                           int interphase,
                           int ap)
{
	int i, k, reaction_number;
	double *propensity;
	double *probability;
	double aggregate;
	double prop_sum = 0;
	double random = g_rand_double (grand);
	int number_of_reactions = 2 * n_target_genes + /* transcription */
		n_target_genes + /* translation */
		n_target_genes + /* degradation mrna */
		n_target_genes; /* degradation protein */
	propensity = g_new0(double, number_of_reactions);
	probability = g_new0(double, number_of_reactions);
	if (interphase == 1) {
/* transcription */
		for (i = 0; i < n_target_genes; i++) {
			double prop = 0; /* Product of T of bound bs */
			int zero = 1;
			for (k = 0; k < n_sites[i]; k++) {
				if (allele_0[ap][i][k].status == 1) {
					prop += T[i * n_tfs + allele_0[ap][i][k].tf_index];
					zero = 0;
				}

			}
			propensity[i] = (zero == 0) ? exp(prop) : 0;
			prop_sum += propensity[i];
		}
		for (i = 0; i < n_target_genes; i++) {
			double prop = 0; /* Product of T of bound bs */
			int zero = 1;
			for (k = 0; k < n_sites[i]; k++) {
				if (allele_1[ap][i][k].status == 1) {
					prop += T[i * n_tfs + allele_1[ap][i][k].tf_index];
					zero = 0;
				}

			}
			propensity[n_target_genes + i] = (zero == 0) ? exp(prop) : 0;
			prop_sum += propensity[n_target_genes + i];
		}
/* translation */
		for (i = 0; i < n_target_genes; i++) {
			k = target_gene_index[i];
			propensity[2 * n_target_genes + i] = solution_mrna[ap * n_tfs + k] * translation[i];
			prop_sum += propensity[2 * n_target_genes + i];
		}
	}
/* mrna degradation */
	for (i = 0; i < n_target_genes; i++) {
		k = target_gene_index[i];
		propensity[3 * n_target_genes + i] = solution_mrna[ap * n_tfs + k] * mrna_degradation[i];
		prop_sum += propensity[3 * n_target_genes + i];
	}
/* protein degradation */
	for (i = 0; i < n_target_genes; i++) {
		k = target_gene_index[i];
		propensity[4 * n_target_genes + i] = solution_protein[ap * n_tfs + k] * protein_degradation[i];
		prop_sum += propensity[4 * n_target_genes + i];
	}
	if (prop_sum <= 0.00000000001) {
		g_warning("slow Sum of propensities is too small %g!", prop_sum);
		(*tau) = TMAX;
		g_free(propensity);
		g_free(probability);
		return(-1);
	}
	reaction_number = -1;
	aggregate = 0;
	if (interphase == 1) {
/* transcription */
		for (i = 0; i < n_target_genes; i++) {
			aggregate += propensity[i];
			probability[i] = aggregate / prop_sum;
			if (random < probability[i]) {
				reaction_number = i;
				(*target) = i;
				(*reaction_type) = 1;
				break;
			}
		}
		if (reaction_number < 0) {
			for (i = 0; i < n_target_genes; i++) {
				aggregate += propensity[n_target_genes + i];
				probability[n_target_genes + i] = aggregate / prop_sum;
				if (random < probability[n_target_genes + i]) {
					reaction_number = n_target_genes + i;
					(*target) = i;
					(*reaction_type) = 2;
					break;
				}
			}
		}
/* translation */
		if (reaction_number < 0) {
			for (i = 0; i < n_target_genes; i++) {
				aggregate += propensity[2 * n_target_genes + i];
				probability[2 * n_target_genes + i] = aggregate / prop_sum;
				if (random < probability[2 * n_target_genes + i]) {
					reaction_number = 2 * n_target_genes + i;
					(*target) = i;
					(*reaction_type) = 3;
					break;
				}
			}
		}
	}
	if (reaction_number < 0) {
/* mrna degradation */
		for (i = 0; i < n_target_genes; i++) {
			aggregate += propensity[3 * n_target_genes + i];
			probability[3 * n_target_genes + i] = aggregate / prop_sum;
			if (random < probability[3 * n_target_genes + i]) {
				reaction_number = 3 * n_target_genes + i;
				(*target) = i;
				(*reaction_type) = 4;
				break;
			}
		}
	}
	if (reaction_number < 0) {
/* protein degradation */
		for (i = 0; i < n_target_genes; i++) {
			aggregate += propensity[4 * n_target_genes + i];
			probability[4 * n_target_genes + i] = aggregate / prop_sum;
			if (random < probability[4 * n_target_genes + i]) {
				reaction_number = 4 * n_target_genes + i;
				(*target) = i;
				(*reaction_type) = 5;
				break;
			}
		}
	}
	(*tau) = -log(g_rand_double (grand)) / prop_sum;
	switch ((*reaction_type)) {
		case 1:
			k = target_gene_index[(*target)];
			solution_mrna[ap * n_tfs + k]++;
		break;
		case 2:
			k = target_gene_index[(*target)];
			solution_mrna[ap * n_tfs + k]++;
		break;
		case 3:
			k = target_gene_index[(*target)];
			solution_protein[ap * n_tfs + k]++;
		break;
		case 4:
			k = target_gene_index[(*target)];
			solution_mrna[ap * n_tfs + k]--;
			if (solution_mrna[ap * n_tfs + k] < 0) {
				g_warning("slow mrna: ap %d tf %d < 0", ap, k);
				solution_mrna[ap * n_tfs + k] = 0;
			}
		break;
		case 5:
			k = target_gene_index[(*target)];
			solution_protein[ap * n_tfs + k]--;
			if (solution_protein[ap * n_tfs + k] < 0) {
				g_warning("slow prot: ap %d tf %d < 0", ap, k);
				solution_protein[ap * n_tfs + k] = 0;
			}
		break;
	}
	g_free(propensity);
	g_free(probability);
	return (reaction_number);
}

int mssa_get_reaction_slow_with_transport(double *solution_mrna,
                                          double *solution_protein,
                                          MSSA_site ***allele_0,
                                          MSSA_site ***allele_1,
                                          double *tau,
                                          int *reaction_type,
                                          int *target,
                                          int *nuc,
                                          int *target_gene_index,
                                          int *n_sites,
                                          double *propensity,
                                          double *probability,
                                          int number_of_reactions_per_nuc,
                                          int number_of_reactions,
                                          double *T,
                                          double *protein_degradation,
                                          double *mrna_degradation,
                                          double *translation,
                                          double *transport_mrna,
                                          double *transport_protein,
                                          int n_tfs,
                                          int n_target_genes,
                                          int n_nucs,
                                          int interphase)
{
	int i, k, reaction_number, ap;
	int reaction_index; /* local reaction type */
	int reaction_target; /* local reaction target */
	int reaction_nuc; /* local reaction ap */
	int found = 0;
	double aggregate;
	double prop_sum = 0;
	double random = g_rand_double (grand);
	if (interphase == 1) {
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(n_nucs, n_target_genes, number_of_reactions_per_nuc, T, translation, solution_mrna, n_tfs, allele_0, allele_1, n_sites, propensity, target_gene_index) reduction(+:prop_sum)
		for (ap = 0; ap < n_nucs; ap++) {
			for (i = 0; i < n_target_genes; i++) {
/* transcription */
				int reaction_number = ap * number_of_reactions_per_nuc + i;
				double prop = 0; /* Product of T of bound bs */
				int zero = 1; /* flag */
				int nact = 0;
				int nrep = 0;
		/* Allele_0 */
				for (int k = 0; k < n_sites[i]; k++) {
					if (allele_0[ap][i][k].status == 1) {
						prop += T[i * n_tfs + allele_0[ap][i][k].tf_index];
						zero = 0;
						nact += (T[i * n_tfs + allele_0[ap][i][k].tf_index] > 0) ? 1:0;
						nrep += (T[i * n_tfs + allele_0[ap][i][k].tf_index] < 0) ? 1:0;
					}

				}
				propensity[reaction_number] = (nact > 0) ? exp(prop) : 0;
				prop_sum += propensity[reaction_number];
		/* Allele_1 */
				reaction_number = ap * number_of_reactions_per_nuc + n_target_genes + i;
				prop = 0; /* Product of T of bound bs */
				zero = 1; /* flag */
				nact = 0;
				nrep = 0;
				for (int k = 0; k < n_sites[i]; k++) {
					if (allele_1[ap][i][k].status == 1) {
						prop += T[i * n_tfs + allele_1[ap][i][k].tf_index];
						zero = 0;
						nact += (T[i * n_tfs + allele_1[ap][i][k].tf_index] > 0) ? 1:0;
						nrep += (T[i * n_tfs + allele_1[ap][i][k].tf_index] < 0) ? 1:0;
					}

				}
				propensity[reaction_number] = (nact > 0) ? exp(prop) : 0;
				prop_sum += propensity[reaction_number];
/* translation */
				reaction_number = ap * number_of_reactions_per_nuc + 2 * n_target_genes + i;
				int k = target_gene_index[i];
				propensity[reaction_number] = solution_mrna[ap * n_tfs + k] * translation[i];
				prop_sum += propensity[reaction_number];
			}
		}
	}
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(solution_mrna, solution_protein, n_nucs, n_target_genes, number_of_reactions_per_nuc, n_tfs, target_gene_index, n_sites, propensity, mrna_degradation, protein_degradation) reduction(+:prop_sum)
	for (ap = 0; ap < n_nucs; ap++) {
		for (i = 0; i < n_target_genes; i++) {
			int reaction_number = ap * number_of_reactions_per_nuc + 3 * n_target_genes + i;
/* mrna degradation */
			int k = target_gene_index[i];
			propensity[reaction_number] = solution_mrna[ap * n_tfs + k] * mrna_degradation[i];
			prop_sum += propensity[reaction_number];
			reaction_number = ap * number_of_reactions_per_nuc + 4 * n_target_genes + i;
/* protein degradation */
			k = target_gene_index[i];
			propensity[reaction_number] = solution_protein[ap * n_tfs + k] * protein_degradation[i];
			prop_sum += propensity[reaction_number];
		}
	}
/* transport */
	ap = 0;
	for (i = 0; i < n_target_genes; i++) {
/* right mrna*/
		reaction_number = ap * number_of_reactions_per_nuc + 6 * n_target_genes + i;
		int k = target_gene_index[i];
		propensity[reaction_number] = (solution_mrna[(ap + 1) * n_tfs + k] - solution_mrna[ap * n_tfs + k]) * transport_mrna[i];
		prop_sum += propensity[reaction_number];
/* right protein */
		reaction_number = ap * number_of_reactions_per_nuc + 8 * n_target_genes + i;
		k = target_gene_index[i];
		propensity[reaction_number] = (solution_protein[(ap + 1) * n_tfs + k] - solution_protein[ap * n_tfs + k]) * transport_protein[i];
		prop_sum += propensity[reaction_number];
	}
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(transport_mrna, transport_protein, solution_mrna, solution_protein, n_nucs, n_target_genes, number_of_reactions_per_nuc, n_tfs, target_gene_index, n_sites, propensity, mrna_degradation, protein_degradation) reduction(+:prop_sum)
	for (ap = 1; ap < n_nucs - 1; ap++) {
		for (i = 0; i < n_target_genes; i++) {
/* left mrna*/
			int reaction_number = ap * number_of_reactions_per_nuc + 5 * n_target_genes + i;
			int k = target_gene_index[i];
			propensity[reaction_number] = (solution_mrna[(ap - 1) * n_tfs + k] - solution_mrna[ap * n_tfs + k]) * transport_mrna[i];
			prop_sum += propensity[reaction_number];
/* right mrna*/
			reaction_number = ap * number_of_reactions_per_nuc + 6 * n_target_genes + i;
			k = target_gene_index[i];
			propensity[reaction_number] = (solution_mrna[(ap + 1) * n_tfs + k] - solution_mrna[ap * n_tfs + k]) * transport_mrna[i];
			prop_sum += propensity[reaction_number];
/* left protein */
			reaction_number = ap * number_of_reactions_per_nuc + 7 * n_target_genes + i;
			k = target_gene_index[i];
			propensity[reaction_number] = (solution_protein[(ap - 1) * n_tfs + k] - solution_protein[ap * n_tfs + k]) * transport_protein[i];
			prop_sum += propensity[reaction_number];
/* right protein */
			reaction_number = ap * number_of_reactions_per_nuc + 8 * n_target_genes + i;
			k = target_gene_index[i];
			propensity[reaction_number] = (solution_protein[(ap + 1) * n_tfs + k] - solution_protein[ap * n_tfs + k]) * transport_protein[i];
			prop_sum += propensity[reaction_number];
		}
	}
	ap = n_nucs - 1;
	for (i = 0; i < n_target_genes; i++) {
/* right mrna*/
		reaction_number = ap * number_of_reactions_per_nuc + 5 * n_target_genes + i;
		int k = target_gene_index[i];
		propensity[reaction_number] = (solution_mrna[(ap - 1) * n_tfs + k] - solution_mrna[ap * n_tfs + k]) * transport_mrna[i];
		prop_sum += propensity[reaction_number];
/* right protein */
		reaction_number = ap * number_of_reactions_per_nuc + 7 * n_target_genes + i;
		k = target_gene_index[i];
		propensity[reaction_number] = (solution_protein[(ap - 1) * n_tfs + k] - solution_protein[ap * n_tfs + k]) * transport_protein[i];
		prop_sum += propensity[reaction_number];
	}
	if (prop_sum <= 0.00000000001) {
		g_warning("slow Sum of propensities is too small %g!", prop_sum);
		(*tau) = TMAX;
		(*target) = -1;
		(*reaction_type) = -1;
		(*nuc) = -1;
		return(-1);
	}
	found = 0;
	aggregate = 0;
	if (interphase == 1) {
		for (ap = 0; ap < n_nucs; ap++) {
			for (i = 0; i < n_target_genes; i++) {
/* transcription 1 */
				reaction_number = ap * number_of_reactions_per_nuc + i;
				aggregate += propensity[reaction_number];
				probability[reaction_number] = aggregate / prop_sum;
				if (random < probability[reaction_number]) {
					reaction_nuc = ap;
					reaction_target = i;
					reaction_index = 1;
					found = 1;
					break;
				}
/* transcription 2 */
				reaction_number = ap * number_of_reactions_per_nuc + n_target_genes + i;
				aggregate += propensity[reaction_number];
				probability[reaction_number] = aggregate / prop_sum;
				if (random < probability[reaction_number]) {
					reaction_nuc = ap;
					reaction_target = i;
					reaction_index = 2;
					found = 1;
					break;
				}
/* translation */
				reaction_number = ap * number_of_reactions_per_nuc + 2 * n_target_genes + i;
				aggregate += propensity[reaction_number];
				probability[reaction_number] = aggregate / prop_sum;
				if (random < probability[reaction_number]) {
					reaction_nuc = ap;
					reaction_target = i;
					reaction_index = 3;
					found = 1;
					break;
				}
			}
			if (found == 1) {
				break;
			}
		}
	}
	if (found == 0) {
		for (ap = 0; ap < n_nucs; ap++) {
			for (i = 0; i < n_target_genes; i++) {
/* mrna degradation */
				reaction_number = ap * number_of_reactions_per_nuc + 3 * n_target_genes + i;
				aggregate += propensity[reaction_number];
				probability[reaction_number] = aggregate / prop_sum;
				if (random < probability[reaction_number]) {
					reaction_nuc = ap;
					reaction_target = i;
					reaction_index = 4;
					found = 1;
					break;
				}
/* protein degradation */
				reaction_number = ap * number_of_reactions_per_nuc + 4 * n_target_genes + i;
				aggregate += propensity[reaction_number];
				probability[reaction_number] = aggregate / prop_sum;
				if (random < probability[reaction_number]) {
					reaction_nuc = ap;
					reaction_target = i;
					reaction_index = 5;
					found = 1;
					break;
				}
			}
			if (found == 1) {
				break;
			}
		}
	}
	if (found == 0) {
		ap = 0;
		for (i = 0; i < n_target_genes; i++) {
/* right mrna*/
			reaction_number = ap * number_of_reactions_per_nuc + 6 * n_target_genes + i;
			aggregate += propensity[reaction_number];
			probability[reaction_number] = aggregate / prop_sum;
			if (random < probability[reaction_number]) {
				reaction_nuc = ap;
				reaction_target = i;
				reaction_index = 7;
				found = 1;
				break;
			}
/* right protein */
			reaction_number = ap * number_of_reactions_per_nuc + 8 * n_target_genes + i;
			aggregate += propensity[reaction_number];
			probability[reaction_number] = aggregate / prop_sum;
			if (random < probability[reaction_number]) {
				reaction_nuc = ap;
				reaction_target = i;
				reaction_index = 9;
				found = 1;
				break;
			}
		}
		if (found == 0) {
			for (ap = 1; ap < n_nucs - 1; ap++) {
				for (i = 0; i < n_target_genes; i++) {
					/* left mrna*/
					reaction_number = ap * number_of_reactions_per_nuc + 5 * n_target_genes + i;
					aggregate += propensity[reaction_number];
					probability[reaction_number] = aggregate / prop_sum;
					if (random < probability[reaction_number]) {
						reaction_nuc = ap;
						reaction_target = i;
						reaction_index = 6;
						found = 1;
						break;
					}
					/* right mrna*/
					reaction_number = ap * number_of_reactions_per_nuc + 6 * n_target_genes + i;
					aggregate += propensity[reaction_number];
					probability[reaction_number] = aggregate / prop_sum;
					if (random < probability[reaction_number]) {
						reaction_nuc = ap;
						reaction_target = i;
						reaction_index = 7;
						found = 1;
						break;
					}
					/* left protein */
					reaction_number = ap * number_of_reactions_per_nuc + 7 * n_target_genes + i;
					aggregate += propensity[reaction_number];
					probability[reaction_number] = aggregate / prop_sum;
					if (random < probability[reaction_number]) {
						reaction_nuc = ap;
						reaction_target = i;
						reaction_index = 8;
						found = 1;
						break;
					}
					/* right protein */
					reaction_number = ap * number_of_reactions_per_nuc + 8 * n_target_genes + i;
					aggregate += propensity[reaction_number];
					probability[reaction_number] = aggregate / prop_sum;
					if (random < probability[reaction_number]) {
						reaction_nuc = ap;
						reaction_target = i;
						reaction_index = 9;
						found = 1;
						break;
					}
				}
				if (found == 1) {
					break;
				}
			}
		}
		if (found == 0) {
			ap = n_nucs - 1;
			for (i = 0; i < n_target_genes; i++) {
				/* right mrna*/
				reaction_number = ap * number_of_reactions_per_nuc + 5 * n_target_genes + i;
				aggregate += propensity[reaction_number];
				probability[reaction_number] = aggregate / prop_sum;
				if (random < probability[reaction_number]) {
					reaction_nuc = ap;
					reaction_target = i;
					reaction_index = 6;
					found = 1;
					break;
				}
				/* right protein */
				reaction_number = ap * number_of_reactions_per_nuc + 7 * n_target_genes + i;
				aggregate += propensity[reaction_number];
				probability[reaction_number] = aggregate / prop_sum;
				if (random < probability[reaction_number]) {
					reaction_nuc = ap;
					reaction_target = i;
					reaction_index = 8;
					found = 1;
					break;
				}
			}
		}
	}
	(*tau) = -log(g_rand_double (grand)) / prop_sum;
	switch (reaction_index) {
		case 1: /* transcription 1 */
			k = target_gene_index[reaction_target];
			solution_mrna[reaction_nuc * n_tfs + k]++;
		break;
		case 2: /* transcription 2 */
			k = target_gene_index[reaction_target];
			solution_mrna[reaction_nuc * n_tfs + k]++;
		break;
		case 3: /* translation */
			k = target_gene_index[reaction_target];
			solution_protein[reaction_nuc * n_tfs + k]++;
		break;
		case 4: /* mrna degradation */
			k = target_gene_index[reaction_target];
			solution_mrna[reaction_nuc * n_tfs + k]--;
			if (solution_mrna[reaction_nuc * n_tfs + k] < 0) {
				g_warning("slow mrna: ap %d tf %d < 0", ap, k);
				solution_mrna[reaction_nuc * n_tfs + k] = 0;
			}
		break;
		case 5: /* protein degradation */
			k = target_gene_index[reaction_target];
			solution_protein[reaction_nuc * n_tfs + k]--;
			if (solution_protein[reaction_nuc * n_tfs + k] < 0) {
				g_warning("slow prot: ap %d tf %d < 0", ap, k);
				solution_protein[reaction_nuc * n_tfs + k] = 0;
			}
		break;
		case 6: /* left transport mrna */
			k = target_gene_index[reaction_target];
			if ((reaction_nuc - 1) > -1) {
				solution_mrna[(reaction_nuc - 1) * n_tfs + k]--;
				solution_mrna[reaction_nuc * n_tfs + k]++;
				if (solution_mrna[(reaction_nuc - 1) * n_tfs + k] < 0) {
					g_warning("slow prot: ap %d tf %d < 0", ap, k);
					solution_mrna[(reaction_nuc - 1) * n_tfs + k] = 0;
				}
			} else {
				g_warning("slow prot: %d ap %d bad %d", reaction_index, reaction_nuc, n_nucs);
			}
		break;
		case 7: /* right transport mrna */
			k = target_gene_index[reaction_target];
			if ((reaction_nuc + 1) < n_nucs) {
				solution_mrna[(reaction_nuc + 1) * n_tfs + k]--;
				solution_mrna[reaction_nuc * n_tfs + k]++;
				if (solution_mrna[(reaction_nuc + 1) * n_tfs + k] < 0) {
					g_warning("slow prot: ap %d tf %d < 0", ap, k);
					solution_mrna[(reaction_nuc + 1) * n_tfs + k] = 0;
				}
			} else {
				g_warning("slow prot: %d ap %d bad %d", reaction_index, reaction_nuc, n_nucs);
			}
		break;
		case 8: /* left transport protein */
			k = target_gene_index[reaction_target];
			if ((reaction_nuc - 1) > -1) {
				solution_protein[(reaction_nuc - 1) * n_tfs + k]--;
				solution_protein[reaction_nuc * n_tfs + k]++;
				if (solution_protein[(reaction_nuc - 1) * n_tfs + k] < 0) {
					g_warning("slow prot: ap %d tf %d < 0", ap, k);
					solution_protein[(reaction_nuc - 1) * n_tfs + k] = 0;
				}
			} else {
				g_warning("slow prot: %d ap %d bad %d", reaction_index, reaction_nuc, n_nucs);
			}
		break;
		case 9: /* right transport protein */
			k = target_gene_index[reaction_target];
			if ((reaction_nuc + 1) < n_nucs) {
				solution_protein[(reaction_nuc + 1) * n_tfs + k]--;
				solution_protein[reaction_nuc * n_tfs + k]++;
				if (solution_protein[(reaction_nuc + 1) * n_tfs + k] < 0) {
					g_warning("slow prot: ap %d tf %d < 0", ap, k);
					solution_protein[(reaction_nuc + 1) * n_tfs + k] = 0;
				}
			} else {
				g_warning("slow prot: %d ap %d bad %d", reaction_index, reaction_nuc, n_nucs);
			}
		break;
	}
	(*target) = reaction_target;
	(*reaction_type) = reaction_index;
	(*nuc) = reaction_nuc;
	return (reaction_number);
}

void mssa_print_timeclass (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	fprintf(stdout, "%d time %f %f protein:\n", problem->repeat, tc->t_start, tc->t_end);
	for (int i = 0; i < tc->n_nucs; i++) {
		for (int j = 0; j < problem->n_tfs; j++) {
			fprintf(stdout, "%f ", tc->solution_protein[i * problem->n_tfs + j]);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "time %f %f bound:\n", tc->t_start, tc->t_end);
	for (int i = 0; i < tc->n_nucs; i++) {
		for (int j = 0; j < problem->n_tfs; j++) {
			fprintf(stdout, "%f ", tc->bound_protein[i * problem->n_tfs + j]);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "time %f %f mrna:\n", tc->t_start, tc->t_end);
	for (int i = 0; i < tc->n_nucs; i++) {
		for (int j = 0; j < problem->n_tfs; j++) {
			fprintf(stdout, "%f ", tc->solution_mrna[i * problem->n_tfs + j]);
		}
		fprintf(stdout, "\n");
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
	for (int i = 0; i < tc->n_nucs; i++) {
		fprintf(fp, "%d %d %d %6.3f ", problem->repeat, tc->kounter, i, tc->t_end);
		for (int j = 0; j < problem->n_target_genes; j++) {
			int k = problem->target_gene_index[j];
			fprintf(fp, "%.0f ", tc->solution_protein[i * problem->n_tfs + k] + tc->bound_protein[i * problem->n_tfs + k]);
		}
		for (int j = 0; j < problem->n_target_genes; j++) {
			int k = problem->target_gene_index[j];
			fprintf(fp, "%.0f ", tc->solution_mrna[i * problem->n_tfs + k]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void add_bias (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa %d add bias %d\n", problem->repeat, tc->kounter);
	for (int i = 0; i < tc->n_nucs; i++) {
		for (int j = 0; j < problem->n_target_genes; j++) {
			int k = problem->target_gene_index[j];
			tc->solution_mrna[i * problem->n_tfs + k] += tc->data_mrna[i * problem->n_tfs + k];
			tc->solution_protein[i * problem->n_tfs + k] += tc->data_protein[i * problem->n_tfs + k] * mol_per_conc;
		}
	}
	if (verbose) mssa_print_timeclass (tc, problem);
}

void score (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	int n = tc->n_nucs;
/*	printf("multiscale_ssa %d %d chisq ", problem->repeat, tc->kounter);
	for (int j = 0; j < problem->n_target_genes; j++) {
		printf("g%d =", j);
		double chisq = 0;
		int k = problem->target_gene_index[j];
#pragma omp parallel for schedule(static) default(none) shared(tc, problem, k, n) reduction(+:chisq)
		for (int i = 0; i < n; i++) {
			double difference = tc->solution_protein[i * problem->n_tfs + k] + tc->bound_protein[i * problem->n_tfs + k] - tc->data_protein[i * problem->n_tfs + k] * MOLECULES_PER_CONCENTRATION;
			chisq += difference * difference;
		}
		printf("%15.6f ", chisq);
	}
	printf("\n");*/
	printf("multiscale_ssa %d %d corr ", problem->repeat, tc->kounter);
	for (int j = 0; j < problem->n_target_genes; j++) {
		printf("g%d =", j);
		int k = problem->target_gene_index[j];
		double sx, sy, mx, my, xx, yy, nx, ny, xy, sxy, cxy;
		sx = sy = mx = my = xx = yy = nx = ny = xy = sxy = cxy = 0;
#pragma omp parallel for schedule(static) default(none) shared(tc, problem, k, n) reduction(+:nx) reduction(+:ny) reduction(+:xx) reduction(+:xy) reduction(+:yy)
		for (int i = 0; i < n; i++) {
			double x = tc->solution_protein[i * problem->n_tfs + k] + tc->bound_protein[i * problem->n_tfs + k];
			double y = tc->data_protein[i * problem->n_tfs + k];
			nx += x;
			ny += y;
			xx += x * x;
			yy += y * y;
			xy += x * y;
		}
/*		mx = nx / (double)n;
		my = ny / (double)n;
		sx = (xx - 2 * mx * nx + mx * mx * n) / (double)(n - 1);
		sy = (yy - 2 * my * ny + my * my * n) / (double)(n - 1);
		sxy = (xy - my * nx - mx * ny + mx * my * n) / (double)(n - 1);
		cxy = (sx > 0 && sy > 0) ? sxy / (sx * sy) : 0;*/
		sx = (double)n * xx - nx * nx;
		sy = (double)n * yy - ny * ny;
		sxy = (double)n * xy - nx * ny;
		cxy = (sx > 0 && sy > 0) ? sxy / sqrt(sx * sy) : 0;
		printf("%15.6f ", 1 - cxy);
/*		printf("%15.6f %15.6f %15.6f ", sx, sy, 1 - cxy);*/
	}
	printf("\n");
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


void propagate (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa propagate %d\n", tc->kounter);
	double t_start_slow;
	double t_stop_slow;
	double t_slow;
/* Set simulation time */
	t_start_slow = tc->t_start;
	t_stop_slow = tc->t_end;
/* Simulate */
	int iter_kounter = 0;
	t_slow = t_start_slow;
	while (t_slow < t_stop_slow) {
		double t_synch_start = t_slow;
		double t_synch_stop = t_slow + (t_stop_slow - t_start_slow) / SYNCH_STEPS_PER_SLOW;
#pragma omp parallel for schedule(static) default(none) shared(problem, tc, t_synch_start, t_synch_stop, fast_time_max)
		for (int ap = 0; ap < tc->n_nucs; ap++) {
			double t_synch = t_synch_start;
			int synch_iter_kounter = 0;
			while (t_synch < t_synch_stop && synch_iter_kounter < MAX_SYNCH_ITER) {
				double tau_synch;
				int reaction_number_slow, inner_iter_kounter;
				int reaction_type_slow, promoter_number_slow;
				double t_start_fast;
				double t_stop_fast;
				double t_fast;
				t_start_fast = 0;
				t_stop_fast = fast_time_max;
				/* Allele_0 */
				t_fast = t_start_fast;
				inner_iter_kounter = 0;
				while (t_fast < t_stop_fast) {
					double tau_fast;
					int promoter_number, tf;
					int site_number;
					int reaction_type;
					int reaction_number;
					reaction_number = mssa_get_reaction_fast(tc->solution_protein,
					                                         tc->bound_protein,
					                                         problem->allele_0,
					                                         &tau_fast,
					                                         &reaction_type,
					                                         &tf,
					                                         &promoter_number,
					                                         &site_number,
					                                         problem->n_sites,
					                                         problem->n_tfs,
					                                         problem->n_target_genes,
					                                         ap);
//					printf("multiscale_ssa %d a0 %f %f %f %d\n", ap, t_synch, t_fast, tau_fast, reaction_number);
					t_fast += tau_fast;
					inner_iter_kounter++;
				}
				/* Allele_1 */
				t_fast = t_start_fast;
				inner_iter_kounter = 0;
				while (t_fast < t_stop_fast) {
					double tau_fast;
					int promoter_number, tf;
					int site_number;
					int reaction_type;
					int reaction_number;
					reaction_number = mssa_get_reaction_fast(tc->solution_protein,
					                                         tc->bound_protein,
					                                         problem->allele_1,
					                                         &tau_fast,
					                                         &reaction_type,
					                                         &tf,
					                                         &promoter_number,
					                                         &site_number,
					                                         problem->n_sites,
					                                         problem->n_tfs,
					                                         problem->n_target_genes,
					                                         ap);
//					printf("multiscale_ssa %d a1 %f %f %f %d\n", ap, t_synch, t_fast, tau_fast, reaction_number);
					t_fast += tau_fast;
					inner_iter_kounter++;
				}
				reaction_number_slow = mssa_get_reaction_slow(tc->solution_mrna,
				                                              tc->solution_protein,
				                                              problem->allele_0,
				                                              problem->allele_1,
				                                              &tau_synch,
				                                              &reaction_type_slow,
				                                              &promoter_number_slow,
				                                              problem->target_gene_index,
				                                              problem->n_sites,
				                                              problem->parameters->T,
				                                              problem->parameters->protein_degradation,
				                                              problem->parameters->mrna_degradation,
				                                              problem->parameters->translation,
				                                              problem->n_tfs,
				                                              problem->n_target_genes,
				                                              1,
				                                              ap);
				printf("multiscale_ssa synch %d %d t %f %f %d %d fast %d\n", ap, synch_iter_kounter, t_synch, tau_synch, reaction_number_slow, reaction_type_slow, inner_iter_kounter);
				t_synch += tau_synch;
				synch_iter_kounter++;
			} /* end of synch loop */
		} /* end of nuc loop */
		t_slow = t_synch_stop;
	} /* end of slow loop */
	mssa_print_timeclass (tc, problem);
}

void propagate_slow_only (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa propagate %d\n", tc->kounter);
	double t_start_slow;
	double t_stop_slow;
	double t_slow;
/* Set simulation time */
	t_start_slow = tc->t_start;
	t_stop_slow = tc->t_end;
/* Simulate */
	int iter_kounter = 0;
	t_slow = t_start_slow;
	while (t_slow < t_stop_slow) {
		double t_synch_start = t_slow;
		double t_synch_stop = t_slow + (t_stop_slow - t_start_slow) / SYNCH_STEPS_PER_SLOW;
#pragma omp parallel for schedule(static) default(none) shared(problem, tc, t_synch_start, t_synch_stop)
		for (int ap = 0; ap < tc->n_nucs; ap++) {
			double t_synch = t_synch_start;
			int synch_iter_kounter = 0;
			while (t_synch < t_synch_stop && synch_iter_kounter < MAX_SYNCH_ITER) {
				double tau_synch;
				int reaction_number_slow;
				int reaction_type_slow, promoter_number_slow;
				reaction_number_slow = mssa_get_reaction_slow(tc->solution_mrna,
				                                              tc->solution_protein,
				                                              problem->allele_0,
				                                              problem->allele_1,
				                                              &tau_synch,
				                                              &reaction_type_slow,
				                                              &promoter_number_slow,
				                                              problem->target_gene_index,
				                                              problem->n_sites,
				                                              problem->parameters->T,
				                                              problem->parameters->protein_degradation,
				                                              problem->parameters->mrna_degradation,
				                                              problem->parameters->translation,
				                                              problem->n_tfs,
				                                              problem->n_target_genes,
				                                              0,
				                                              ap);
				printf("multiscale_ssa msynch %d %d t %f %f %d %d\n", ap, synch_iter_kounter, t_synch, tau_synch, reaction_number_slow, reaction_type_slow);
				t_synch += tau_synch;
				synch_iter_kounter++;
			} /* end of synch loop */
		} /* end of nuc loop */
		t_slow = t_synch_stop;
	} /* end of slow loop */
	if (verbose) mssa_print_timeclass (tc, problem);
}

void propagate_with_transport (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa %d propagate %d\n", problem->repeat, tc->kounter);
	double t_start_slow;
	double t_stop_slow;
	double t_slow, tau_slow;
	double *propensity;
	double *probability;
	int number_of_reactions_per_nuc;
	int number_of_reactions;
	number_of_reactions_per_nuc = 2 * problem->n_target_genes + /* transcription */
		problem->n_target_genes + /* translation */
		problem->n_target_genes + /* degradation mrna */
		problem->n_target_genes + /* degradation protein */
		problem->n_target_genes + /* left transport mrna */
		problem->n_target_genes + /* right transport mrna */
		problem->n_target_genes + /* left transport protein */
		problem->n_target_genes; /* right transport protein */
	number_of_reactions = number_of_reactions_per_nuc * tc->n_nucs; // - 2 * problem->n_target_genes;
	int slow_only = (tc->type == 3) ? 1 : 0;
	propensity = g_new0(double, number_of_reactions);
	probability = g_new0(double, number_of_reactions);
	double **propensity_fast;
	double **probability_fast;
	int number_of_reactions_fast_per_nuc = 4 * problem->n_tfs * problem->n_target_genes; /* bind/unbind each tf in each promotor */
	propensity_fast = g_new0(double*, tc->n_nucs);
	probability_fast = g_new0(double*, tc->n_nucs);
	int ***site_tab_bind_0 = g_new0(int**, tc->n_nucs);
	int ***site_tab_bind_1 = g_new0(int**, tc->n_nucs);
	int ***site_tab_unbind_0 = g_new0(int**, tc->n_nucs);
	int ***site_tab_unbind_1 = g_new0(int**, tc->n_nucs);
	for (int ap = 0; ap < tc->n_nucs; ap++) {
		propensity_fast[ap] = g_new0(double, number_of_reactions_fast_per_nuc);
		probability_fast[ap] = g_new0(double, number_of_reactions_fast_per_nuc);
		site_tab_bind_0[ap] = g_new0(int*, problem->n_target_genes);
		site_tab_bind_1[ap] = g_new0(int*, problem->n_target_genes);
		site_tab_unbind_0[ap] = g_new0(int*, problem->n_target_genes);
		site_tab_unbind_1[ap] = g_new0(int*, problem->n_target_genes);
		for (int i = 0; i < problem->n_target_genes; i++) {
			site_tab_bind_0[ap][i] = g_new0(int, problem->n_tfs);
			site_tab_bind_1[ap][i] = g_new0(int, problem->n_tfs);
			site_tab_unbind_0[ap][i] = g_new0(int, problem->n_tfs);
			site_tab_unbind_1[ap][i] = g_new0(int, problem->n_tfs);
		}
	}
/* Set simulation time */
	t_start_slow = tc->t_start;
	t_stop_slow = tc->t_end;
/* Simulate */
	int iter_kounter = 0;
	t_slow = t_start_slow;
	while (t_slow < t_stop_slow) {
		int reaction_number_slow, inner_iter_kounter = 0;
		int reaction_type_slow, promoter_number_slow;
		int nuc_number_slow;
		if (slow_only == 0) { /* interphase */
#pragma omp parallel for schedule(static) default(none) shared(problem, tc, propensity_fast, probability_fast, site_tab_bind_0, site_tab_bind_1, site_tab_unbind_0, site_tab_unbind_1,fast_time_max, fast_steps_frac) reduction(+:inner_iter_kounter)
			for (int ap = 0; ap < tc->n_nucs; ap++) {
				double t_start_fast;
				double t_stop_fast;
				double t_fast;
				t_start_fast = 0;
				t_stop_fast = fast_time_max;
				t_fast = t_start_fast;
				for (int l = 0; l < (int)(problem->sum_sites * fast_steps_frac); l++) {
/*				while (t_fast < t_stop_fast) {*/
					double tau_fast;
					int promoter_number, tf;
					int site_number;
					int reaction_type;
					int reaction_number;
					reaction_number = mssa_get_reaction_fast_with_buffers_2 (tc->solution_protein,
					                                                         tc->bound_protein,
					                                                         problem->allele_0,
					                                                         problem->allele_1,
					                                                         propensity_fast[ap],
					                                                         probability_fast[ap],
									                                         site_tab_bind_0[ap],
									                                         site_tab_bind_1[ap],
									                                         site_tab_unbind_0[ap],
									                                         site_tab_unbind_1[ap],
					                                                         &tau_fast,
					                                                         &reaction_type,
					                                                         &tf,
					                                                         &promoter_number,
					                                                         &site_number,
					                                                         problem->n_sites,
					                                                         problem->n_tfs,
					                                                         problem->n_target_genes,
					                                                         problem->parameters->T,
					                                                         problem->parameters->range,
					                                                         ap);
					//					printf("multiscale_ssa %d a0 %f %f %f %d\n", ap, t_synch, t_fast, tau_fast, reaction_number);
					t_fast += tau_fast;
					inner_iter_kounter++;
					if (t_fast > t_stop_fast) break;
				}
			}
		}
		reaction_number_slow = mssa_get_reaction_slow_with_transport (tc->solution_mrna,
		                                                              tc->solution_protein,
		                                                              problem->allele_0,
		                                                              problem->allele_1,
		                                                              &tau_slow,
		                                                              &reaction_type_slow,
		                                                              &promoter_number_slow,
		                                                              &nuc_number_slow,
		                                                              problem->target_gene_index,
		                                                              problem->n_sites,
		                                                              propensity,
		                                                              probability,
		                                                              number_of_reactions_per_nuc,
		                                                              number_of_reactions,
		                                                              problem->parameters->T,
		                                                              problem->parameters->protein_degradation,
		                                                              problem->parameters->mrna_degradation,
		                                                              problem->parameters->translation,
		                                                              problem->parameters->transport_mrna,
		                                                              problem->parameters->transport_protein,
		                                                              problem->n_tfs,
		                                                              problem->n_target_genes,
		                                                              tc->n_nucs,
		                                                              1 - slow_only);
		iter_kounter++;
		t_slow += tau_slow;
/* Print the configuration of reg region to the log file */
		if (log_file != NULL && g_strcmp0 (log_file, "nolog")) {
			FILE*fp = fopen(log_file, "a");
			fprintf(fp, "%d %d slow %d %f %f %d %d %d %d", problem->repeat, tc->kounter, iter_kounter, t_slow, tau_slow, reaction_number_slow, nuc_number_slow, promoter_number_slow, reaction_type_slow);
			if (slow_only == 0) {
				fprintf(fp, " fast %d allele_0", inner_iter_kounter);
				for (int ap = 0; ap < tc->n_nucs; ap++) {
					for(int i = 0; i < problem->n_target_genes; i++) {
						for(int j = 0; j < problem->n_sites[i]; j++) {
							fprintf(fp, " %d", problem->allele_0[ap][i][j].status);
						}
					}
				}
				fprintf(fp, " allele_1");
				for (int ap = 0; ap < tc->n_nucs; ap++) {
					for(int i = 0; i < problem->n_target_genes; i++) {
						for(int j = 0; j < problem->n_sites[i]; j++) {
							fprintf(fp, " %d", problem->allele_1[ap][i][j].status);
						}
					}
				}
			}
			fprintf(fp, "\n");
			fclose (fp);
		} else if (log_file == NULL) {
/* if log file is not given print log to the screen */
			fprintf(stderr, "%d %d slow %d %f %f %d %d %d %d", problem->repeat, tc->kounter, iter_kounter, t_slow, tau_slow, reaction_number_slow, nuc_number_slow, promoter_number_slow, reaction_type_slow);
			if (slow_only == 0) {
				fprintf(stderr, " fast %d", inner_iter_kounter);
			}
			fprintf(stderr, "\n");
		}
	} /* end of slow loop */
	for (int ap = 0; ap < tc->n_nucs; ap++) {
		g_free (propensity_fast[ap]);
		g_free (probability_fast[ap]);
		for (int i = 0; i < problem->n_target_genes; i++) {
			g_free(site_tab_bind_0[ap][i]);
			g_free(site_tab_bind_1[ap][i]);
			g_free(site_tab_unbind_0[ap][i]);
			g_free(site_tab_unbind_1[ap][i]);
		}
		g_free(site_tab_bind_0[ap]);
		g_free(site_tab_bind_1[ap]);
		g_free(site_tab_unbind_0[ap]);
		g_free(site_tab_unbind_1[ap]);
	}
	g_free(propensity_fast);
	g_free(probability_fast);
	g_free(propensity);
	g_free(probability);
	g_free(site_tab_bind_0);
	g_free(site_tab_bind_1);
	g_free(site_tab_unbind_0);
	g_free(site_tab_unbind_1);
	if (verbose) mssa_print_timeclass (tc, problem);
}

/* External inputs
 * are added
 */

void inject (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa inject %d\n", tc->kounter);
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(problem, tc, mol_per_conc)
	for (int i = 0; i < tc->n_nucs; i++) {
		for (int j = 0; j < problem->n_external_genes; j++) {
			int k = problem->external_gene_index[j];
			double conc = tc->data_protein[i * problem->n_tfs + k] * mol_per_conc - tc->bound_protein[i * problem->n_tfs + k];
			tc->solution_protein[i * problem->n_tfs + k] = MAX(conc, 0);
		}
	}
	if (verbose) mssa_print_timeclass (tc, problem);
}

void zero_structure (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa zero_structure %d\n", tc->kounter);
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(problem, tc)
	for (int i = 0; i < tc->n_nucs; i++) {
		for (int j = 0; j < problem->n_tfs; j++) {
			tc->solution_protein[i * problem->n_tfs + j] = 0;
			tc->solution_mrna[i * problem->n_tfs + j] = 0;
			tc->bound_protein[i * problem->n_tfs + j] = 0;
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
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(problem, tc)
	for (int i = 0; i < tc->n_nucs; i++) {
		for (int j = 0; j < problem->n_tfs; j++) {
			tc->solution_protein[i * problem->n_tfs + j] += tc->bound_protein[i * problem->n_tfs + j];
			tc->bound_protein[i * problem->n_tfs + j] = 0;
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

/* Nucleaar division
 * All protein is assumed to be already unbound!
 */

void divide (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa divide %d\n", tc->kounter);
	MSSA_Timeclass *tc_prev = (MSSA_Timeclass *) g_list_nth_data (problem->tc_list, (tc->kounter - 1));
	if (!tc_prev) return;
	if (tc->n_nucs > tc_prev->n_nucs * 2) {
		int i = 0;
		for (int j = 0; j < problem->n_tfs; j++) {
			tc->solution_mrna[2 * i * problem->n_tfs + j] = tc_prev->solution_mrna[i * problem->n_tfs + j] / 2;
			tc->solution_protein[2 * i * problem->n_tfs + j] = tc_prev->solution_protein[i * problem->n_tfs + j] / 2;
		}
		for (int i = 1; i < tc_prev->n_nucs - 1; i++) {
			for (int j = 0; j < problem->n_tfs; j++) {
				tc->solution_mrna[2 * i * problem->n_tfs + j] = tc_prev->solution_mrna[i * problem->n_tfs + j] / 2;
				tc->solution_mrna[(2 * i + 1) * problem->n_tfs + j] = tc_prev->solution_mrna[i * problem->n_tfs + j] / 2;
				tc->solution_protein[2 * i * problem->n_tfs + j] = tc_prev->solution_protein[i * problem->n_tfs + j] / 2;
				tc->solution_protein[(2 * i + 1) * problem->n_tfs + j] = tc_prev->solution_protein[i * problem->n_tfs + j] / 2;
			}
		}
		i = tc_prev->n_nucs - 1;
		for (int j = 0; j < problem->n_tfs; j++) {
			tc->solution_mrna[(2 * i + 1) * problem->n_tfs + j] = tc_prev->solution_mrna[i * problem->n_tfs + j] / 2;
			tc->solution_protein[(2 * i + 1) * problem->n_tfs + j] = tc_prev->solution_protein[i * problem->n_tfs + j] / 2;
		}
	} else {
		for (int i = 0; i < tc_prev->n_nucs; i++) {
			for (int j = 0; j < problem->n_tfs; j++) {
				tc->solution_mrna[2 * i * problem->n_tfs + j] = tc_prev->solution_mrna[i * problem->n_tfs + j] / 2;
				tc->solution_mrna[(2 * i + 1) * problem->n_tfs + j] = tc_prev->solution_mrna[i * problem->n_tfs + j] / 2;
				tc->solution_protein[2 * i * problem->n_tfs + j] = tc_prev->solution_protein[i * problem->n_tfs + j] / 2;
				tc->solution_protein[(2 * i + 1) * problem->n_tfs + j] = tc_prev->solution_protein[i * problem->n_tfs + j] / 2;
			}
		}
	}
	if (verbose) mssa_print_timeclass (tc, problem);
}

void connect (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	if (verbose) printf("multiscale_ssa connect %d\n", tc->kounter);
	MSSA_Timeclass *tc_prev = (MSSA_Timeclass *) g_list_nth_data (problem->tc_list, (tc->kounter - 1));
	if (!tc_prev) return;
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(problem, tc, tc_prev)
	for (int i = 0; i < tc_prev->n_nucs; i++) {
		for (int j = 0; j < problem->n_tfs; j++) {
			tc->solution_mrna[i * problem->n_tfs + j] = tc_prev->solution_mrna[i * problem->n_tfs + j];
			tc->solution_protein[i * problem->n_tfs + j] = tc_prev->solution_protein[i * problem->n_tfs + j];
			tc->bound_protein[i * problem->n_tfs + j] = tc_prev->bound_protein[i * problem->n_tfs + j];
		}
	}
	if (verbose) mssa_print_timeclass (tc, problem);
}

void mssa_mark_site_overlap (MSSA_Problem *problem, int range)
{
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(range, problem)
	for(int k = 0; k < problem->n_nucs; k++) {
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
}

void integrate (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
/* Debug tag */
	if (verbose) printf("multiscale_ssa %d integrate %d\n", problem->repeat, tc->kounter);
/* Print initial condition */
	if (tc->kounter == 0 && verbose) mssa_print_timeclass (tc, problem);
/* Copy result from previous TC */
	if (tc->type > 0) connect (tc, problem);
/* Add bias or set the initial cond */
	if (tc->type == 2) {
		inject (tc, problem);
		add_bias (tc, problem);
		if (!dryrun) propagate_with_transport (tc, problem);
	}
/* Run the model */
	if (tc->type == 1) {
		if (tc->has_data == 1) {
			inject (tc, problem);
			score (tc, problem);
		}
		if (!dryrun) propagate_with_transport (tc, problem);
	}
/* Nuclear division, All protein is to be unbound already */
	if (tc->type == 0) divide (tc, problem);
/* Run the model in the mitosis mode - without fast reactions,
 * translation and transcription.
 * Firstly, unbound all proteins!
 */
	if (tc->type == 3) {
		unbound (tc, problem);
		if (!dryrun) propagate_with_transport (tc, problem);
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
	{ "repeat", 'r', 0, G_OPTION_ARG_INT, &repeat, N_("Number of repeats"), N_("REPEAT") },
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
		for (int r = 0; r < repeat - 1; r++) {
			problem->repeat = r;
			g_list_foreach (problem->tc_list, (GFunc) integrate, (gpointer) problem);
			zero_structure (g_list_nth_data (problem->tc_list, 0), problem);
		}
		problem->repeat = repeat - 1;
		g_list_foreach (problem->tc_list, (GFunc) integrate, (gpointer) problem);
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
		while (!feof(stdout)) {
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
			for (int r = 0; r < repeat - 1; r++) {
				problem->repeat = r;
				g_list_foreach (problem->tc_list, (GFunc) integrate, (gpointer) problem);
				zero_structure (g_list_nth_data (problem->tc_list, 0), problem);
			}
			problem->repeat = repeat - 1;
			g_list_foreach (problem->tc_list, (GFunc) integrate, (gpointer) problem);
			printf("\r\nfflush\r\n");
			fflush (stdout);
		}
	} else {
		for (int r = 0; r < repeat - 1; r++) {
			problem->repeat = r;
			g_list_foreach (problem->tc_list, (GFunc) integrate, (gpointer) problem);
			zero_structure (g_list_nth_data (problem->tc_list, 0), problem);
		}
		problem->repeat = repeat - 1;
		g_list_foreach (problem->tc_list, (GFunc) integrate, (gpointer) problem);
	}
	return (0);
}
