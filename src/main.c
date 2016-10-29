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
#include <glib.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <limits.h>

#define TMAX 1000000
#define LOWBITS 0x330E
#define BYTESIZE 8

#define SECONDS_PER_DAY 86400

long seed = 1;
unsigned short ysubj[3];

int ssa_prepare(int *species, int *stoichiometric_matrix, double *tau, double *parameter, int number_of_species, int number_of_reactions)
{
	int i, j, reaction_number;
	double *propensity;
	double *probability;
	double prop_sum = 0;
	double random = erand48(ysubj);
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
	(*tau) = -log(erand48(ysubj)) / prop_sum;
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
	if (reaction_number < 0) return;
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
				break;
				case 1:
/*					species[i] = (sum > 0) ? 1:0;*/
/*					species[i] = (inhibitors < activators) ? 1:0;*/
					species[i] = (inhibitors == 0 && 0 < activators) ? 1:0;
				break;
			}
		}
	}
}

int main(int argc, char**argv)
{
	printf("multiscale_ssa\n");
	int number_of_reactions_slow;
	int number_of_reactions_fast;
	int number_of_species_slow;
	int number_of_species_fast;
	int number_of_binding_sites;
	int number_of_promoters;
	int *M;
	double*parameter_slow;
	double*parameter_fast;
	int*stoichiometric_matrix_slow;
	int*stoichiometric_matrix_fast;
	int *reaction_type_fast;
	double t_start_slow;
	double t_stop_slow;
	double t_start_fast;
	double t_stop_fast;
	int iter_kounter, inner_iter_kounter;
	double t_slow, t_fast;
	int *common_spices_indices;
	int *slow_species_promoters_indices;
	int *slow_matrix_reaction;
	int *fast_species_promoters_indices;
	int *fast_reactions_promoters_indices;
	int *slow_species_promoters_indices_type;
	int *fast_species_promoters_indices_type;
	int *species_slow;
	int *species_fast;
	int i;

	ysubj[0] = LOWBITS;
	ysubj[1] = (unsigned short)seed;
	ysubj[2] = (unsigned short) ( seed >> (BYTESIZE * sizeof(unsigned short)) );

	number_of_reactions_slow = 11;
	number_of_reactions_fast = 6;
	number_of_species_slow = 7;
	number_of_species_fast = 7;
	number_of_binding_sites = 20;
	number_of_promoters = 2;
	parameter_slow = g_new0(double, number_of_reactions_slow);
	parameter_fast = g_new0(double, number_of_reactions_fast);
	stoichiometric_matrix_slow = g_new0(int, number_of_reactions_slow * number_of_species_slow);
	slow_matrix_reaction = g_new0(int, number_of_reactions_slow);
	stoichiometric_matrix_fast = g_new0(int, number_of_reactions_fast * number_of_species_fast);
	reaction_type_fast = g_new0(int, number_of_reactions_fast);
	fast_reactions_promoters_indices = g_new0(int, number_of_reactions_fast);
	t_start_slow = 0;
	t_stop_slow = 20 * SECONDS_PER_DAY;
	t_start_fast = 0;
	t_stop_fast = 1000 * SECONDS_PER_DAY;
	common_spices_indices = g_new0(int, number_of_species_slow);
	slow_species_promoters_indices = g_new0(int, number_of_species_slow);
	slow_species_promoters_indices_type = g_new0(int, number_of_species_slow);
	fast_species_promoters_indices = g_new0(int, number_of_species_fast);
	fast_species_promoters_indices_type = g_new0(int, number_of_species_fast);
	species_slow = g_new0(int, number_of_species_slow);
	species_fast = g_new0(int, number_of_species_fast);
/*
 * Slow species
 * CpAST, mRNA16, QpAST, mRNA1, EBNA1, dimEBNA1, Oct2
 */
	for (i = 0; i < number_of_species_slow; i++) {
		common_spices_indices[i] = -1;
		slow_species_promoters_indices[i] = -1;
		slow_species_promoters_indices_type[i] = 1;
		species_slow[i] = 0;
	}
	species_slow[6] = 1.5 * 1e4;
	species_slow[4] = 0.125 * 1e4;
/*
 * Fast species
 * Cp, EBNA1, CpAST, Oct2, CpI, Qp, QpI
 */
	for (i = 0; i < number_of_species_fast; i++) {
		species_fast[i] = 0;
		fast_species_promoters_indices[i] = -1;
		fast_species_promoters_indices_type[i] = 1;
	}
	common_spices_indices[4] = 1;
	common_spices_indices[6] = 3;
	slow_species_promoters_indices[0] = 0;
	slow_species_promoters_indices[2] = 1;
	slow_species_promoters_indices_type[0] = 1;
	slow_species_promoters_indices_type[2] = 1;
	fast_species_promoters_indices[0] = 0;
	fast_species_promoters_indices[2] = 0;
	fast_species_promoters_indices[4] = 0;
	fast_species_promoters_indices[5] = 1;
	fast_species_promoters_indices[6] = 1;
	fast_species_promoters_indices_type[0] = 0;
	fast_species_promoters_indices_type[2] = 1;
	fast_species_promoters_indices_type[4] = -1;
	fast_species_promoters_indices_type[5] = 0;
	fast_species_promoters_indices_type[6] = -1;
/*
 * -2 - unbind an inhibitor
 * -1 - bind an inhibitor
 * 0 - not binding or unbinding
 * +1 - bind activator
 * +2 - unbind activator
 */
	reaction_type_fast[0] = 1;
	reaction_type_fast[1] = 2;
	reaction_type_fast[2] = -1;
	reaction_type_fast[3] = -2;
	reaction_type_fast[4] = -1;
	reaction_type_fast[5] = -2;
	fast_reactions_promoters_indices[0] = 0;
	fast_reactions_promoters_indices[1] = 0;
	fast_reactions_promoters_indices[2] = 0;
	fast_reactions_promoters_indices[3] = 0;
	fast_reactions_promoters_indices[4] = 1;
	fast_reactions_promoters_indices[5] = 1;
/*
 * Slow species
 * Reaction, CpAST, mRNA16, QpAST, mRNA1, EBNA1, dimEBNA1, Oct2
 * 0       , -1   , 1,      0,      0,    0,     0,         0
 * 1       , 0   , 0,      -1,      1,    0,     0,         0
 * 2       , 0   , -1,      0,      0,    1,     0,         0
 * 3       , 0   , 0,      0,      -1,    1,     0,         0
 * 4       , 0   , 0,      0,      0,    -2,     1,         0
 * 5       , 0   , 0,      0,      0,    2,     -1,         0
 * 6       , 0   , -1,      0,      0,    0,     0,         0
 * 7       , 0   , 0,      0,      -1,    0,     0,         0
 * 8       , 0   , 0,      0,      0,    -1,     0,         0
 * 9       , 0   , 0,      0,      0,    0,     -1,         0
 * 10      , 0   , 0,      0,      0,    0,     0,         -1
 */

	stoichiometric_matrix_slow[0] = -1;stoichiometric_matrix_slow[1] = 1;stoichiometric_matrix_slow[2] = 0;stoichiometric_matrix_slow[3] = 0;stoichiometric_matrix_slow[4] = 0;stoichiometric_matrix_slow[5] = 0;stoichiometric_matrix_slow[6] = 0;
	stoichiometric_matrix_slow[7] = 0;stoichiometric_matrix_slow[8] = 0;stoichiometric_matrix_slow[9] = -1;stoichiometric_matrix_slow[10] = 1;stoichiometric_matrix_slow[11] = 0;stoichiometric_matrix_slow[12] = 0;stoichiometric_matrix_slow[13] = 0;
	stoichiometric_matrix_slow[14] = 0;stoichiometric_matrix_slow[15] = -1;stoichiometric_matrix_slow[16] = 0;stoichiometric_matrix_slow[17] = 0;stoichiometric_matrix_slow[18] = 1;stoichiometric_matrix_slow[19] = 0;stoichiometric_matrix_slow[20] = 0;
	stoichiometric_matrix_slow[21] = 0;stoichiometric_matrix_slow[22] = 0;stoichiometric_matrix_slow[23] = 0;stoichiometric_matrix_slow[24] = -1;stoichiometric_matrix_slow[25] = 1;stoichiometric_matrix_slow[26] = 0;stoichiometric_matrix_slow[27] = 0;
	stoichiometric_matrix_slow[28] = 0;stoichiometric_matrix_slow[29] = 0;stoichiometric_matrix_slow[30] = 0;stoichiometric_matrix_slow[31] = 0;stoichiometric_matrix_slow[32] = -2;stoichiometric_matrix_slow[33] = 1;stoichiometric_matrix_slow[34] = 0;
	stoichiometric_matrix_slow[35] = 0;stoichiometric_matrix_slow[36] = 0;stoichiometric_matrix_slow[36] = 0;stoichiometric_matrix_slow[37] = 0;stoichiometric_matrix_slow[38] = 2;stoichiometric_matrix_slow[39] = -1;stoichiometric_matrix_slow[41] = 0;
	stoichiometric_matrix_slow[42] = 0;stoichiometric_matrix_slow[43] = -1;stoichiometric_matrix_slow[44] = 0;stoichiometric_matrix_slow[45] = 0;stoichiometric_matrix_slow[46] = 0;stoichiometric_matrix_slow[47] = 0;stoichiometric_matrix_slow[48] = 0;
	stoichiometric_matrix_slow[49] = 0;stoichiometric_matrix_slow[50] = 0;stoichiometric_matrix_slow[51] = 0;stoichiometric_matrix_slow[52] = -1;stoichiometric_matrix_slow[53] = 0;stoichiometric_matrix_slow[54] = 0;stoichiometric_matrix_slow[55] = 0;
	stoichiometric_matrix_slow[56] = 0;stoichiometric_matrix_slow[57] = 0;stoichiometric_matrix_slow[58] = 0;stoichiometric_matrix_slow[59] = 0;stoichiometric_matrix_slow[60] = -1;stoichiometric_matrix_slow[61] = 0;stoichiometric_matrix_slow[62] = 0;
	stoichiometric_matrix_slow[63] = 0;stoichiometric_matrix_slow[64] = 0;stoichiometric_matrix_slow[65] = 0;stoichiometric_matrix_slow[66] = 0;stoichiometric_matrix_slow[67] = 0;stoichiometric_matrix_slow[68] = -1;stoichiometric_matrix_slow[69] = 0;
	stoichiometric_matrix_slow[70] = 0;stoichiometric_matrix_slow[71] = 0;stoichiometric_matrix_slow[72] = 0;stoichiometric_matrix_slow[73] = 0;stoichiometric_matrix_slow[74] = 0;stoichiometric_matrix_slow[75] = 0;stoichiometric_matrix_slow[76] = -1;

	slow_matrix_reaction[0] = 1;
	slow_matrix_reaction[1] = 1;
	slow_matrix_reaction[2] = 1;
	slow_matrix_reaction[3] = 1;
	slow_matrix_reaction[4] = 0;
	slow_matrix_reaction[5] = 0;
	slow_matrix_reaction[6] = 0;
	slow_matrix_reaction[7] = 0;
	slow_matrix_reaction[8] = 0;
	slow_matrix_reaction[9] = 0;
	slow_matrix_reaction[10] = 0;
	
/*
 * Fast species
 * Reaction, Cp, EBNA1, CpAST, Oct2, CpI, Qp, QpI
 * 0       , -1,   -1,      1,    0,   0,  0,  0
 * 1       ,  1,    1,     -1,    0,   0,  0,  0
 * 2       , -1,    0,      0 ,  -1,   1,  0,  0
 * 3       ,  1,    0,      0,    1,  -1,  0,  0
 * 4       ,  0,   -1,      0,    0,   0, -1,  1
 * 5       ,  0,    1,      0,    0,   0,  1, -1
 */
stoichiometric_matrix_fast[0] = -1;stoichiometric_matrix_fast[1] = -1;stoichiometric_matrix_fast[2] = 1;stoichiometric_matrix_fast[3] = 0;stoichiometric_matrix_fast[4] = 0;stoichiometric_matrix_fast[5] = 0;stoichiometric_matrix_fast[6] = 0;
stoichiometric_matrix_fast[7] = 1;stoichiometric_matrix_fast[8] = 1;stoichiometric_matrix_fast[9] = -1;stoichiometric_matrix_fast[10] = 0;stoichiometric_matrix_fast[11] = 0;stoichiometric_matrix_fast[12] = 0;stoichiometric_matrix_fast[13] = 0;
stoichiometric_matrix_fast[14] = -1;stoichiometric_matrix_fast[15] = 0;stoichiometric_matrix_fast[16] = 0;stoichiometric_matrix_fast[17] = -1;stoichiometric_matrix_fast[18] = 1;stoichiometric_matrix_fast[19] = 0;stoichiometric_matrix_fast[20] = 0;
stoichiometric_matrix_fast[21] = 1;stoichiometric_matrix_fast[22] = 0;stoichiometric_matrix_fast[23] = 0;stoichiometric_matrix_fast[24] = 1;stoichiometric_matrix_fast[25] = -1;stoichiometric_matrix_fast[26] = 0;stoichiometric_matrix_fast[27] = 0;
stoichiometric_matrix_fast[28] = 0;stoichiometric_matrix_fast[29] = -1;stoichiometric_matrix_fast[30] = 0;stoichiometric_matrix_fast[31] = 0;stoichiometric_matrix_fast[32] = 0;stoichiometric_matrix_fast[33] = -1;stoichiometric_matrix_fast[34] = 1;
stoichiometric_matrix_fast[35] = 0;stoichiometric_matrix_fast[36] = 1;stoichiometric_matrix_fast[37] = 0;stoichiometric_matrix_fast[38] = 0;stoichiometric_matrix_fast[39] = 0;stoichiometric_matrix_fast[40] = 1;stoichiometric_matrix_fast[41] = -1;

	parameter_slow[0] = 0.002069;
	parameter_slow[1] = 0.006240;
	parameter_slow[2] = 0.005172;
	parameter_slow[3] = 0.0156;
	parameter_slow[4] = 1.8492 * 1e-10;
	parameter_slow[5] = 1e-8;
	parameter_slow[6] = 1.9254 * 1e-5;
	parameter_slow[7] = 1.9254 * 1e-5;
	parameter_slow[8] = 9.627 * 1e-6;
	parameter_slow[9] = 9.627 * 1e-6;
	parameter_slow[10] = 9.627 * 1e-6;
	parameter_fast[0] = 9.2462 * 1e-12;
	parameter_fast[1] = 1.5 * 1e-11;
	parameter_fast[2] = 9.2462 * 1e-12;
	parameter_fast[3] = 2.5 * 1e-9;
	parameter_fast[4] = 9.2462 * 1e-12;
	parameter_fast[5] = 2.1 * 1e-10;
	M = g_new0(int, number_of_promoters * number_of_binding_sites);
	for (i = 0; i < number_of_promoters * number_of_binding_sites; i++) {
		M[i] = 0;
	}
	iter_kounter = 0;
	t_slow = t_start_slow;
	int flag = 0;
	while (t_slow < t_stop_slow) {
		double tau_slow;
		int reaction_number_slow;
// Quantify the number of activated, repressed,
// and free binding sites for each type of promoters
// x_2[p] = EvaluateBindingSites(M_p(j))
//		evaluate_binding_sites(species_fast, M);
		if (t_slow > 0.5 * t_stop_slow && flag == 0) {
			species_slow[6] += 5 * 1e5;
			flag = 1;
		}
		evaluate_promotor_state(species_fast, M, fast_species_promoters_indices, number_of_species_fast, number_of_binding_sites, fast_species_promoters_indices_type);
		for (i = 0; i < number_of_species_slow; i++) {
			if (common_spices_indices[i] > 0) {
				species_fast[common_spices_indices[i]] = species_slow[i];
			}
		}
		t_fast = t_start_fast;
		inner_iter_kounter = 0;
		while (t_fast < t_stop_fast) {
			double tau_fast;
			int promoter_number;
			int site_number;
			int reaction_type;
			int reaction_number;
			int k, l, s;
			reaction_number = ssa_prepare(species_fast, stoichiometric_matrix_fast, &tau_fast, parameter_fast, number_of_species_fast, number_of_reactions_fast);
			reaction_type = reaction_type_fast[reaction_number];
			promoter_number = fast_reactions_promoters_indices[reaction_number];
			site_number = -1;
			switch (reaction_type) {
				case -2: /* Unbind an inhibitor */
					l = floor(erand48(ysubj) * number_of_binding_sites);
					s = l;
					for (k = 0; k < number_of_binding_sites; k++) {
						/* Looking for a bound inhbitor */
						if (M[promoter_number * number_of_binding_sites + s] == -1) {
							site_number = s;
							break;
						}
						s++;
						if (s > number_of_binding_sites - 1) {
							s = 0;
						}
					}
					if (site_number > -1) {
						M[promoter_number * number_of_binding_sites + site_number] = 0;
					}
				break;
				case -1: /* Bind an inhibitor */
					l = floor(erand48(ysubj) * number_of_binding_sites);
					s = l;
					for (k = 0; k < number_of_binding_sites; k++) {
						/* Looking for a free site */
						if (M[promoter_number * number_of_binding_sites + s] == 0) {
							site_number = s;
							break;
						}
						s++;
						if (s > number_of_binding_sites - 1) {
							s = 0;
						}
					}
					if (site_number > -1) {
						M[promoter_number * number_of_binding_sites + site_number] = -1;
					}
				break;
				case 1: /* Bind activator */
					l = floor(erand48(ysubj) * number_of_binding_sites);
					s = l;
					for (k = 0; k < number_of_binding_sites; k++) {
						/* Looking for a free site */
						if (M[promoter_number * number_of_binding_sites + s] == 0) {
							site_number = s;
							break;
						}
						s++;
						if (s > number_of_binding_sites - 1) {
							s = 0;
						}
					}
					if (site_number > -1) {
						M[promoter_number * number_of_binding_sites + site_number] = 1;
					}
				break;
				case 2: /* Unbind an activator */
					l = floor(erand48(ysubj) * number_of_binding_sites);
					s = l;
					for (k = 0; k < number_of_binding_sites; k++) {
						/* Looking for a bound activator */
						if (M[promoter_number * number_of_binding_sites + s] == 1) {
							site_number = s;
							break;
						}
						s++;
						if (s > number_of_binding_sites - 1) {
							s = 0;
						}
					}
					if (site_number > -1) {
						M[promoter_number * number_of_binding_sites + site_number] = 0;
					}
				break;
			}
			if (site_number > -1) {
/*				for (i = 0; i < number_of_species_fast; i++) {
					if (i != 1 && i != 3) continue;	
					species_fast[i] += stoichiometric_matrix_fast[reaction_number * number_of_species_fast + i];
					if (species_fast[i] < 0) {
						g_warning("Specie fast %d < 0", i);
						species_fast[i] = 0;
					}
				}*/
				for (i = 0; i < number_of_species_slow; i++) {
					k = common_spices_indices[i];
					if (k > 0) {
						species_fast[k] += stoichiometric_matrix_fast[reaction_number * number_of_species_fast + k];
						if (species_fast[k] < 0) {
							g_warning("Specie fast %d < 0", k);
							species_fast[k] = 0;
						}
					}
				}
			}
/*
			fprintf(stdout, "%f %d %f %f fast", t_slow, inner_iter_kounter, t_fast, tau_fast);
			for (i = 0; i < number_of_species_fast; i++) {
				fprintf(stdout, " %d", species_fast[i]);
			}
			fprintf(stdout, " reaction %d %d %d %d\n", reaction_number, reaction_type, promoter_number, site_number);
*/  
			t_fast += tau_fast;
			inner_iter_kounter++;
		}
		evaluate_promotor_state(species_slow, M, slow_species_promoters_indices, number_of_species_slow, number_of_binding_sites, slow_species_promoters_indices_type);
		for (i = 0; i < number_of_species_slow; i++) {
			if (common_spices_indices[i] > 0) {
				species_slow[i] = species_fast[common_spices_indices[i]];
			}
		}
		reaction_number_slow = ssa_step(species_slow, stoichiometric_matrix_slow, &tau_slow, parameter_slow, number_of_species_slow, number_of_reactions_slow, slow_matrix_reaction);
		fprintf(stdout, "%d %d %f %f slow", iter_kounter, reaction_number_slow, t_slow/SECONDS_PER_DAY, t_slow);
		for (i = 0; i < number_of_species_slow; i++) {
			fprintf(stdout, " %d", species_slow[i]);
		}
		fprintf(stdout, " fast");
		for (i = 0; i < number_of_species_fast; i++) {
			fprintf(stdout, " %d", species_fast[i]);
		}
		fprintf(stdout, "\n");
		t_slow += tau_slow;
		iter_kounter++;
	}
	return (0);
}
