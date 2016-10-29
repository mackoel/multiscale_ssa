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

#define TMAX 10000
#define LOWBITS 0x330E
#define BYTESIZE 8

long seed = 1;
unsigned short ysubj[3];

int ssa_prepare(int *spieces, int *stoichiometric_matrix, double *tau, double *parameter, int number_of_spieces, int number_of_reactions)
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
		for (j = 0; j < number_of_spieces; j++) {
			if (stoichiometric_matrix[i * number_of_spieces + j] > 0) {
				zero = 0;
				prop *= spieces[j];
			}
		}
		propensity[i] = (zero == 1) ? 0 : (parameter[i] * prop);
		prop_sum += propensity[i];
	}
	if (prop_sum <= 0.00000000001) {
		(*tau) = TMAX;
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
	return (reaction_number);
}

void ssa_step(int *spieces, int *stoichiometric_matrix, double *tau, double *parameter, int number_of_spieces, int number_of_reactions)
{
	int j, reaction_number;
	reaction_number = ssa_prepare(spieces, stoichiometric_matrix, tau, parameter, number_of_spieces, number_of_reactions);
	for (j = 0; j < number_of_spieces; j++) {
		spieces[j] += stoichiometric_matrix[reaction_number * number_of_spieces + j];
		if (spieces[j] == 0) {
			g_warning("Specie %d < 0", j);
		}
		spieces[j] = 0;
	}
}

void evaluate_promotor_state(int *spieces, int *M, int *spieces_promoters_indices, int number_of_spieces, int number_of_binding_sites, int *opposite_spieces)
{
	int i, j;
	for (i = 0; i < number_of_spieces; i++) {
		int k = spieces_promoters_indices[i];
		if (k > 0) {
			double sum = 0;
			for (j = 0; j < number_of_binding_sites; j++) {
				sum += M[k * number_of_binding_sites + j];
			}
			spieces[i] = (opposite_spieces[i] * sum > 0) ? 1:0;
		}
	}
}

int main(int argc, char**argv)
{
	printf("multiscale_ssa\n");
	int number_of_reactions_slow;
	int number_of_reactions_fast;
	int number_of_spieces_slow;
	int number_of_spieces_fast;
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
	int iter_kounter;
	double t_slow, t_fast;
	int *common_spices_indices;
	int *slow_spieces_promoters_indices;
	int *fast_spieces_promoters_indices;
	int *fast_reactions_promoters_indices;
	int *slow_spieces_promoters_indices_opposite;
	int *fast_spieces_promoters_indices_opposite;
	int *spieces_slow;
	int *spieces_fast;
	int i;

	ysubj[0] = LOWBITS;
	ysubj[1] = (unsigned short)seed;
	ysubj[2] = (unsigned short) ( seed >> (BYTESIZE * sizeof(unsigned short)) );

	number_of_reactions_slow = 11;
	number_of_reactions_fast = 6;
	number_of_spieces_slow = 7;
	number_of_spieces_fast = 7;
	number_of_binding_sites = 20;
	number_of_promoters = 2;
	parameter_slow = g_new0(double, number_of_reactions_slow);
	parameter_fast = g_new0(double, number_of_reactions_fast);
	stoichiometric_matrix_slow = g_new0(int, number_of_reactions_slow * number_of_spieces_slow);
	stoichiometric_matrix_fast = g_new0(int, number_of_reactions_fast * number_of_spieces_fast);
	reaction_type_fast = g_new0(int, number_of_reactions_fast);
	fast_reactions_promoters_indices = g_new0(int, number_of_reactions_fast);
	t_start_slow = 0;
	t_stop_slow = 2000;
	t_start_fast = 0;
	t_stop_fast = 100;
	common_spices_indices = g_new0(int, number_of_spieces_slow);
	slow_spieces_promoters_indices = g_new0(int, number_of_spieces_slow);
	slow_spieces_promoters_indices_opposite = g_new0(int, number_of_spieces_slow);
	fast_spieces_promoters_indices = g_new0(int, number_of_spieces_fast);
	fast_spieces_promoters_indices_opposite = g_new0(int, number_of_spieces_fast);
	spieces_slow = g_new0(int, number_of_spieces_slow);
	spieces_fast = g_new0(int, number_of_spieces_fast);
/*
 * Slow spieces
 * CpAST, mRNA16, QpAST, mRNA1, EBNA1, dimEBNA1, Oct2
 */
	for (i = 0; i < number_of_spieces_slow; i++) {
		common_spices_indices[i] = -1;
		slow_spieces_promoters_indices[i] = -1;
		slow_spieces_promoters_indices_opposite[i] = 1;
		spieces_slow[i] = 0;
	}
/*
 * Fast spieces
 * Cp, EBNA1, CpAST, Oct2, CpI, Qp, QpI
 */
	for (i = 0; i < number_of_spieces_fast; i++) {
		spieces_fast[i] = 0;
		fast_spieces_promoters_indices[i] = -1;
		fast_spieces_promoters_indices_opposite[i] = 1;
	}
	common_spices_indices[4] = 1;
	common_spices_indices[6] = 3;
	slow_spieces_promoters_indices[0] = 0;
	slow_spieces_promoters_indices[2] = 1;
	fast_spieces_promoters_indices[0] = 0;
	fast_spieces_promoters_indices[5] = 1;
	fast_spieces_promoters_indices_opposite[4] = -1;
	fast_spieces_promoters_indices_opposite[6] = -1;
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
 * Slow spieces
 * Reaction, CpAST, mRNA16, QpAST, mRNA1, EBNA1, dimEBNA1, Oct2
 * 0       , 0   , 1,      0,      0,    0,     0,         0
 * 1       , 0   , 0,      0,      1,    0,     0,         0
 * 2       , 0   , 0,      0,      0,    1,     0,         0
 * 3       , 0   , 0,      0,      0,    1,     0,         0
 * 4       , 0   , 0,      0,      0,    -2,     1,         0
 * 5       , 0   , 0,      0,      0,    2,     -1,         0
 * 6       , 0   , -1,      0,      0,    0,     0,         0
 * 7       , 0   , 0,      0,      -1,    0,     0,         0
 * 8       , 0   , 0,      0,      0,    -1,     0,         0
 * 9       , 0   , 0,      0,      0,    0,     -1,         0
 * 10      , 0   , 0,      0,      0,    0,     0,         -1
 */

 stoichiometric_matrix_slow[0] = 0;stoichiometric_matrix_slow[1] = 1;stoichiometric_matrix_slow[2] = 0;stoichiometric_matrix_slow[3] = 0;stoichiometric_matrix_slow[4] = 0;stoichiometric_matrix_slow[5] = 0;stoichiometric_matrix_slow[6] = 0;
 stoichiometric_matrix_slow[7] = 0;stoichiometric_matrix_slow[8] = 0;stoichiometric_matrix_slow[9] = 0;stoichiometric_matrix_slow[10] = 1;stoichiometric_matrix_slow[11] = 0;stoichiometric_matrix_slow[12] = 0;stoichiometric_matrix_slow[13] = 0;
 stoichiometric_matrix_slow[14] = 0;stoichiometric_matrix_slow[15] = 0;stoichiometric_matrix_slow[16] = 0;stoichiometric_matrix_slow[17] = 0;stoichiometric_matrix_slow[18] = 1;stoichiometric_matrix_slow[19] = 0;stoichiometric_matrix_slow[20] = 0;
 stoichiometric_matrix_slow[21] = 0;stoichiometric_matrix_slow[22] = 0;stoichiometric_matrix_slow[23] = 0;stoichiometric_matrix_slow[24] = 0;stoichiometric_matrix_slow[25] = 1;stoichiometric_matrix_slow[26] = 0;stoichiometric_matrix_slow[27] = 0;
 stoichiometric_matrix_slow[28] = 0;stoichiometric_matrix_slow[29] = 0;stoichiometric_matrix_slow[30] = 0;stoichiometric_matrix_slow[31] = 0;stoichiometric_matrix_slow[32] = -2;stoichiometric_matrix_slow[33] = 1;stoichiometric_matrix_slow[34] = 0;
 stoichiometric_matrix_slow[35] = 0;stoichiometric_matrix_slow[36] = 0;stoichiometric_matrix_slow[36] = 0;stoichiometric_matrix_slow[37] = 0;stoichiometric_matrix_slow[38] = 2;stoichiometric_matrix_slow[39] = -1;stoichiometric_matrix_slow[41] = 0;
 stoichiometric_matrix_slow[42] = 0;stoichiometric_matrix_slow[43] = -1;stoichiometric_matrix_slow[44] = 0;stoichiometric_matrix_slow[45] = 0;stoichiometric_matrix_slow[46] = 0;stoichiometric_matrix_slow[47] = 0;stoichiometric_matrix_slow[48] = 0;
 stoichiometric_matrix_slow[49] = 0;stoichiometric_matrix_slow[50] = 0;stoichiometric_matrix_slow[51] = 0;stoichiometric_matrix_slow[52] = -1;stoichiometric_matrix_slow[53] = 0;stoichiometric_matrix_slow[54] = 0;stoichiometric_matrix_slow[55] = 0;
 stoichiometric_matrix_slow[56] = 0;stoichiometric_matrix_slow[57] = 0;stoichiometric_matrix_slow[58] = 0;stoichiometric_matrix_slow[59] = 0;stoichiometric_matrix_slow[60] = -1;stoichiometric_matrix_slow[61] = 0;stoichiometric_matrix_slow[62] = 0;
 stoichiometric_matrix_slow[63] = 0;stoichiometric_matrix_slow[64] = 0;stoichiometric_matrix_slow[65] = 0;stoichiometric_matrix_slow[66] = 0;stoichiometric_matrix_slow[67] = 0;stoichiometric_matrix_slow[68] = -1;stoichiometric_matrix_slow[69] = 0;
 stoichiometric_matrix_slow[70] = 0;stoichiometric_matrix_slow[71] = 0;stoichiometric_matrix_slow[72] = 0;stoichiometric_matrix_slow[73] = 0;stoichiometric_matrix_slow[74] = 0;stoichiometric_matrix_slow[75] = 0;stoichiometric_matrix_slow[76] = -1;
/*
 * Fast spieces
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
	while (t_slow < t_stop_slow) {
		double tau_slow;
// Quantify the number of activated, repressed,
// and free binding sites for each type of promoters
// x_2[p] = EvaluateBindingSites(M_p(j))
//		evaluate_binding_sites(spieces_fast, M);
		evaluate_promotor_state(spieces_fast, M, fast_spieces_promoters_indices, number_of_spieces_fast, number_of_binding_sites, fast_spieces_promoters_indices_opposite);
		t_fast = t_start_fast;
		while (t_fast < t_stop_fast) {
			double tau_fast;
			int promoter_number;
			int site_number;
			int reaction_type;
			int reaction_number;
			int k, l, s;
			reaction_number = ssa_prepare(spieces_fast, stoichiometric_matrix_fast, &tau_fast, parameter_fast, number_of_spieces_fast, number_of_reactions_fast);
			reaction_type = reaction_type_fast[reaction_number];
			promoter_number = fast_reactions_promoters_indices[reaction_number];
			switch (reaction_type) {
				case -2: /* Unbind an inhibitor */
					l = floor(erand48(ysubj) * number_of_binding_sites);
					s = l;
					site_number = -1;
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
					site_number = -1;
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
					site_number = -1;
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
					site_number = -1;
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
			M[promoter_number * number_of_binding_sites + site_number] = reaction_type;
			t_fast += tau_fast;
		}
		evaluate_promotor_state(spieces_slow, M, slow_spieces_promoters_indices, number_of_spieces_slow, number_of_binding_sites, slow_spieces_promoters_indices_opposite);
		for (i = 0; i < number_of_spieces_slow; i++) {
			if (common_spices_indices[i] > 0) {
				spieces_slow[i] = spieces_fast[common_spices_indices[i]];
			}
		}
		ssa_step(spieces_slow, stoichiometric_matrix_slow, &tau_slow, parameter_slow, number_of_spieces_slow,number_of_reactions_slow);
		t_slow += tau_slow;
		fprintf(stdout, "%d %f", iter_kounter, t_slow);
		for (i = 0; i < number_of_spieces_slow; i++) {
			fprintf(stdout, " %d", spieces_slow[i]);
		}
		for (i = 0; i < number_of_spieces_fast; i++) {
			fprintf(stdout, " %d", spieces_fast[i]);
		}
		fprintf(stdout, "\n");
	}
	return (0);
}

