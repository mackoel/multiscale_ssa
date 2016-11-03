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
#define MOLECULES_PER_CONCENTRATION 1e2
#define SYNCH_STEPS_PER_SLOW 20
#define FAST_TIME_MAX 0.50

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
	int status; /* bound = 1, free = 0, blocked = -1 */
} MSSA_site;

typedef struct {
	double *T; /* activation coefficient */
	double *protein_degradation; /* n_target_genes length */
	double *mrna_degradation; /* n_target_genes length */
	double *translation;
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
	MSSA_site ***allele_0; /* sites on the first allele */
	MSSA_site ***allele_1; /* sites on the second allele */
	MSSA_Parameters *parameters;
} MSSA_Problem;

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
	problem->n_sites = g_new0(int, problem->n_target_genes);
	for (int i = 0; i < problem->n_target_genes; i++) {
		fscanf(fp, "%d", &(problem->n_sites[i]));
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
	for (int i = 0; i < problem->n_target_genes; i++) {
		for (int j = 0; j < problem->n_tfs; j++) {
			fscanf(fp, "%lf", &(problem->parameters->T[i * problem->n_tfs + j]));
		}
		fscanf(fp, "%lf", &(problem->parameters->translation[i]));
		fscanf(fp, "%lf", &(problem->parameters->protein_degradation[i]));
		fscanf(fp, "%lf", &(problem->parameters->mrna_degradation[i]));
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

void mssa_print_timeclass (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	fprintf(stdout, "time %f %f protein:\n", tc->t_start, tc->t_end); 
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

void add_bias (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	printf("multiscale_ssa add bias %d\n", tc->kounter);
	for (int i = 0; i < tc->n_nucs; i++) {
		for (int j = 0; j < problem->n_target_genes; j++) {
			int k = problem->target_gene_index[j];
			tc->solution_mrna[i * problem->n_tfs + k] += tc->data_mrna[i * problem->n_tfs + k];
			tc->solution_protein[i * problem->n_tfs + k] += tc->data_protein[i * problem->n_tfs + k] * MOLECULES_PER_CONCENTRATION;
		}
	}
	mssa_print_timeclass (tc, problem);
}

void score (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	printf("multiscale_ssa score %d", tc->kounter);
	double score = 0;
	for (int i = 0; i < tc->n_nucs; i++) {
		for (int j = 0; j < problem->n_target_genes; j++) {
			int k = problem->target_gene_index[j];
			double difference = tc->solution_protein[i * problem->n_tfs + k] - tc->data_protein[i * problem->n_tfs + k] * MOLECULES_PER_CONCENTRATION; 
			score += difference * difference;
		}
	}
	printf("multiscale_ssa score %d=%15.6f", tc->kounter, score);
}

void propagate (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	printf("multiscale_ssa propagate %d\n", tc->kounter);
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
				int reaction_number_slow, inner_iter_kounter;
				int reaction_type_slow, promoter_number_slow;
				double t_start_fast;
				double t_stop_fast;
				double t_fast;
				t_start_fast = 0;
				t_stop_fast = FAST_TIME_MAX;
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
	printf("multiscale_ssa propagate %d\n", tc->kounter);
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
	mssa_print_timeclass (tc, problem);
}

/* External inputs
 * are added
 */

void inject (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	printf("multiscale_ssa inject %d\n", tc->kounter);
	for (int i = 0; i < tc->n_nucs; i++) {
		for (int j = 0; j < problem->n_external_genes; j++) {
			int k = problem->external_gene_index[j];
			double conc = tc->data_protein[i * problem->n_tfs + k] * MOLECULES_PER_CONCENTRATION - tc->bound_protein[i * problem->n_tfs + k];
			tc->solution_protein[i * problem->n_tfs + k] = MAX(conc, 0);
		}
	}
	mssa_print_timeclass (tc, problem);
}

void unbound (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	printf("multiscale_ssa unbound %d\n", tc->kounter);
	for (int i = 0; i < tc->n_nucs; i++) {
		for (int j = 0; j < problem->n_tfs; j++) {
			tc->solution_protein[i * problem->n_tfs + j] += tc->bound_protein[i * problem->n_tfs + j];
			tc->bound_protein[i * problem->n_tfs + j] = 0;
		}
	}
	mssa_print_timeclass (tc, problem);
}

/* Nucleaar division
 * All protein is assumed to be already unbound!
 */

void divide (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	printf("multiscale_ssa divide %d\n", tc->kounter);
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
	for(int k = 1; k < problem->n_nucs; k++) {
		for (int i = 0; i < problem->n_target_genes; i++) {
			for (int j = 0; j < problem->n_sites[i]; j++) {
				problem->allele_0[k][i][j].status = 0;
				problem->allele_1[k][i][j].status = 0;
			}
		}
	}
	mssa_print_timeclass (tc, problem);
}

void connect (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
	printf("multiscale_ssa connect %d\n", tc->kounter);
	MSSA_Timeclass *tc_prev = (MSSA_Timeclass *) g_list_nth_data (problem->tc_list, (tc->kounter - 1));
	if (!tc_prev) return;
	for (int i = 0; i < tc_prev->n_nucs; i++) {
		for (int j = 0; j < problem->n_tfs; j++) {
			tc->solution_mrna[i * problem->n_tfs + j] = tc_prev->solution_mrna[i * problem->n_tfs + j];
			tc->solution_protein[i * problem->n_tfs + j] = tc_prev->solution_protein[i * problem->n_tfs + j];
			tc->bound_protein[i * problem->n_tfs + j] = tc_prev->bound_protein[i * problem->n_tfs + j];
		}
	}
	mssa_print_timeclass (tc, problem);
}

void integrate (MSSA_Timeclass *tc, MSSA_Problem *problem)
{
/* Debug tag */
	printf("multiscale_ssa integrate %d\n", tc->kounter);
/* Print initial condition */
	if (tc->kounter == 0) mssa_print_timeclass (tc, problem);
/* Copy result from previous TC */
	if (tc->type > 0) connect (tc, problem);
/* Add bias or set the initial cond */ 
	if (tc->type == 2) {
		inject (tc, problem);
		add_bias (tc, problem);
		propagate (tc, problem);
	}
/* Run the model */
	if (tc->type == 1) {
		inject (tc, problem);
		if (tc->has_data == 1) score (tc, problem);
		propagate (tc, problem);
	}
/* Nuclear division, All protein is to be unbound already */
	if (tc->type == 0) divide (tc, problem);
/* Run the model in the mitosis mode - without fast reactions, 
 * translation of transcription.
 * Firstly, unbound all proteins!
 */
	if (tc->type == 3) {
		unbound (tc, problem);
		propagate_slow_only (tc, problem);
	}
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

static char*data_file;
static char*log_file;
static char*operation;

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
	{ "datafile", 0, 0, G_OPTION_ARG_STRING, &data_file, N_("File name for everything and default group names"), N_("FILENAME") },
	{ "logfile", 0, 0, G_OPTION_ARG_STRING, &log_file, N_("File name where model is described"), N_("FILENAME") },
	{ "operation", 0, 0, G_OPTION_ARG_STRING, &operation, N_("What to do"), N_("OPERATION") },
	{ "version", 0, G_OPTION_FLAG_NO_ARG | G_OPTION_FLAG_HIDDEN, G_OPTION_ARG_CALLBACK, option_version_cb, NULL, NULL },
	{ NULL }
};

int main(int argc, char**argv)
{
	printf("multiscale_ssa start\n");
	GOptionContext *context;
	GError *gerror = NULL;
	context = g_option_context_new (_("- DEEP optimizer"));
	g_option_context_add_main_entries(context, (const GOptionEntry *)entries, NULL);
	g_option_context_set_ignore_unknown_options(context, TRUE);
	if (!g_option_context_parse (context, &argc, &argv, &gerror)) {
		g_error (_("option parsing failed: %s\n"), gerror->message);
	}
	g_option_context_free (context);
	if (data_file == NULL) {
		g_error(_("%s called with wrong options for model"), g_get_prgname());
	}
	MSSA_Problem *problem = mssa_read_problem(data_file);
	grand = g_rand_new ();
	printf("multiscale_ssa read problem\n");
	printf("multiscale_ssa nnucs %d\n", problem->n_nucs);
	printf("multiscale_ssa tfs %d\n", problem->n_tfs);
	printf("multiscale_ssa targets %d\n", problem->n_target_genes);
	g_list_foreach (problem->tc_list, (GFunc) integrate, (gpointer) problem);
	return (0);
}
