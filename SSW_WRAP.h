#pragma once

#include "ssw.h"

static const int8_t nt_table[128] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

class SSW_LIB
{
public:
	// this is for ssw
	s_profile* profile;
	int8_t* num;	// the read sequence represented in numbers
	int8_t* ref_num;	// the read sequence represented in numbers
	s_align* result;
	int8_t* mat;
	/* This table is used to transform nucleotide letters into numbers. */

	SSW_LIB(int32_t match, int32_t mismatch, int32_t len)
	{
		int32_t l, m, k;

		num = (int8_t*)malloc(len);	// the read sequence represented in numbers
		ref_num = (int8_t*)malloc(len);	// the read sequence represented in numbers

											// initialize scoring matrix for genome sequences
											//  A  C  G  T	N (or other ambiguous code)
											//  2 -2 -2 -2 	0	A
											// -2  2 -2 -2 	0	C
											// -2 -2  2 -2 	0	G
											// -2 -2 -2  2 	0	T
											//	0  0  0  0  0	N (or other ambiguous code)
		mat = (int8_t*)calloc(25, sizeof(int8_t));
		for (l = k = 0; l < 4; ++l) {
			for (m = 0; m < 4; ++m) mat[k++] = l == m ? match : -mismatch;	/* weight_match : -weight_mismatch */
			mat[k++] = 0; // ambiguous base: no penalty
		}
		for (m = 0; m < 5; ++m) mat[k++] = 0;
	}

	int32 SSW(char *read_seq, char *ref_seq, int32_t len, int32_t gap_open, int32_t gap_extension)
	{
		int32_t m;
		// reference sequence
		//static const char ref_seq[40] = { 'C', 'A', 'G', 'C', 'C', 'T', 'T', 'T', 'C', 'T', 'G', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'A', 'T',
		//	'C', 'A', 'A', 'A', 'A', 'T', 'A', 'G', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'A', 'A', 'A', '\0' };
		//static const char read_seq[16] = { 'C', 'T', 'G', 'A', 'G', 'C', 'C', 'G', 'G', 'T', 'A', 'A', 'A', 'T', 'C', '\0' };	// read sequence

		for (m = 0; m < len; ++m) num[m] = nt_table[(int)read_seq[m]];
		profile = ssw_init(num, len, mat, 5, 2);
		for (m = 0; m < len; ++m) ref_num[m] = nt_table[(int)ref_seq[m]];

		result = ssw_align(profile, ref_num, len, gap_open, gap_extension, 1, 0, 0, len);
		return (int32)result->score1;
		//ssw_write(result, ref_seq, read_seq, nt_table);
	}

	~SSW_LIB()
	{
		align_destroy(result);
		init_destroy(profile);
		free(mat);
		free(ref_num);
		free(num);
	}

};
