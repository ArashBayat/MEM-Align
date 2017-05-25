#pragma once

#include "define.h"
#include "seqence.h"
#include "mem.h"
#include "SSW_WRAP.h"

enum ALGO
{
	MEM_ALIGN,
	SMITH_WATERMAN
};

class CIGAR
{
public:
	char str[500];
	int str_pos;
	char last_added;
	int last_cnt;
	CIGAR()
	{
		str_pos = 0;
		last_added = ' ';
		last_cnt = 0;
	}
	void Add(char ch)
	{
		if (ch == last_added)
			last_cnt++;
		else
		{
			if (last_cnt != 0)
			{
				str_pos += sprintf(&str[str_pos], "%d%c", last_cnt, last_added);
			}
			last_added = ch;
			last_cnt = 1;
		}
	}
	void Close()
	{
		if (last_cnt != 0)
		{
			str_pos += sprintf(&str[str_pos], "%d%c", last_cnt, last_added);
		}
		str[str_pos] = 0;
		str_pos = 0;
		last_added = ' ';
		last_cnt = 0;
	}
};

class ALIGNER
{
public:
	SSW_LIB *SSW;
	ALGO processed_by;

	// T and Q should be of the same lenght in this implementation of MEM-Align
	SEQUENCE T;		// Target Sequence
	SEQUENCE Q;		// Query Sequence
	int_t len;		// length of sequences 

	#define Rm 1	// match score
	int_t Px;		// mismatch penalty
	int_t Po;		// gap open penalty
	int_t Pe;		// gap extend penalty
	
	int_t sl;		// short MEM threashold
	int_t gl;		// banded alignment limits
	int_t TD;		// MEM distance threshold
	int_t TM;		// Maximum number of MEMs threashold
	int_t TS;		// minimum optimal alignment score threashold

	//int_t mmlioa;	// max mismatch len in optimal alignment. (to be reliable)
	//bool not_reliable; 	// not reliable if there is a mismatch len over mmlioa in optimal alignment.
	
	MEM *mems;	// Extracted MEMs are stored in this array
	MEM *smems;	// Sorted    MEMs are stored in this array
	int_t num_mems; // number of MEMs in the array
	
	// optimal alignment (oa)
	int_t	oa_score;
	int_t	oa_last_mem_idx;
	int_t	oa_first_mem_idx;
	int_t	oa_num_mems; // the number of MEMs that forms optimal alignment (not total number of MEMs)

	// after extending the first and the last mem of the optimal alignment
	// the begining and ending of alignment in T and Q are stored in these variables
	int_t	oa_target_start_pos;
	int_t	oa_query_start_pos;
	int_t	oa_target_end_pos;
	int_t	oa_query_end_pos;


	int_t num_ofs; // 2*gl+1; when gap limited optimisation is applied this is number of row in Figure 7 of the paper
	
	// it is the length of padding aroung bitvector of T
	// it is consider to avoid losing T information when we shift it to left ot right
	int_t pad_word; // ceil(gl/BASE_PER_WORD); 

	// there is one flag for each row of Figure 7 of the paper.
	// Basicaly when ofs_flag(ofs) it is set to (i) for (Mi) this means all member 
	// of the set OMEGA in H^i_ofs are processed. Thus the computation of
	// S^j_i for MEMs with OFS=ofs is not necessary (should be avoided)
	int_t *ofs_flag;

	SEQUENCE mask_seq;	// Used in MEM extraction. See Section 5 of the paper (look for the word mask)
	SEQUENCE cmp_seq;	// E variable in Algorithm 1 of the paper
	
	// Temporary data for counting sort
	int_t *csort_cummulative_cnt;
	int_t *csort_temp;

	//statistics
	uint64 bypass_SSW_by_TM;			// number of sequence pairs bypassed to Smith-Waterman by TM
	uint64 bypass_SSW_by_TS;			// number of sequence pairs bypassed to Smith-Waterman by optimal alignemnt score
	uint64 bypass_SSW_by_noMEM;			// number of sequence pairs bypassed to Smith-Waterman by num_MEM=0
	uint64 bypass_SSW_total;			// number of sequence pairs bypassed to Smith-Waterman total
	//int_t bypass_SW_by_mmlioa;		// number of sequence pairs bypassed to Smith-Waterman by mmlioa
	uint64 num_SJI_normal;				// number of SJI to be computed
	uint64 num_SJI_Omega;				// number of SJI after Omega optimisation
	uint64 num_SJI_Omega_TD;				// number of SJI after Distant MEM optimisation
	uint64 numNeededSeqCmp;				// number of time needed to look back REG and compare sequences
	uint64 numSeqCmpAfterPreCompute;	// number of sequence compare after the optimisation that compute minimum PIJ

	CIGAR cstr;

	~ALIGNER()
	{
		//delete[] mems;
		//delete[] smems;
		//delete[] ofs_flag;
		//delete[] csort_cummulative_cnt;
		//delete[] csort_temp;
	}

	// the lenght of seuences and the size of array of MEM is given for initialisation
 	ALIGNER(int_t p_seq_len, int_t max_mem)
	{
		len = p_seq_len;
		mems = new MEM[max_mem];
		smems = new MEM[max_mem];
		ofs_flag = new int_t[len *2+1+10];
		csort_cummulative_cnt = new int_t[len +10];
		csort_temp = new int_t[len +10];

		// BWA Default values
		//match_score = 1;
		Px = 4;
		Po = 6;
		Pe = 1;
		
		sl = 5;
		gl = 5;
		TD = 20;
		TM = 100;

		num_ofs = (2 * gl) + 1;
		pad_word = gl / BASE_PER_WORD;
		if (gl % BASE_PER_WORD)
			pad_word++;

		num_mems = 0;

		oa_score = 0;
		oa_num_mems = 0;
		oa_last_mem_idx = -1;
		oa_first_mem_idx = -1;

		bypass_SSW_by_TM = 0;
		bypass_SSW_by_TS = 0;
		bypass_SSW_by_noMEM = 0;
		bypass_SSW_total = 0;
		num_SJI_normal = 0;
		num_SJI_Omega = 0;
		num_SJI_Omega_TD = 0;
		numNeededSeqCmp = 0;
		numSeqCmpAfterPreCompute = 0;

		T.Init(len, pad_word);
		mask_seq.Init(len, pad_word);
		Q.Init(len, 0);
		temp_seq[0].Init(len, 0);
		temp_seq[1].Init(len, 0);
		cmp_seq.Init(len, 0);
	}

	inline void Set_Penalties(int_t mismatch, int_t gap_open, int_t gap_extend)
	{
		Px = mismatch;
		Po = gap_open;
		Pe = gap_extend;
		SSW = new SSW_LIB(Rm, Px, len);
	}

	inline void Set_Parameter(int_t p_gl, int_t p_sl, int_t p_TD, int_t p_TM)
	{
		Set_Parameter(p_gl, p_sl, p_TD, p_TM, 0);
	}

	inline void Set_Parameter(int_t p_gl, int_t p_sl, int_t p_TD, int_t p_TM, int_t p_TS)
	{
		gl = p_gl;
		sl = p_sl;
		TD = p_TD;
		TM = p_TM;
		TS = p_TS;

		num_ofs = (2 * gl) + 1;
		pad_word = gl / BASE_PER_WORD;
		if (gl % BASE_PER_WORD)
			pad_word++;

		T.Init(len, pad_word);
		mask_seq.Init(len, pad_word);
	}

	// extract MEMs for all ofset
	inline void Extract_MEM()
	{
		mask_seq.Fill_Match();

		//target.Print_Seq();
		T.Long_Shift_Right_Logical(gl * BIT_PER_BASE);
		//target.Print_Seq();
		//mask_seq.Print_Hex_All();
		mask_seq.Long_Shift_Right_Insert_1(gl * BIT_PER_BASE);
		//mask_seq.Print_Hex_All();

		num_mems = 0;
		int_t current_offset = -gl;
		int_t current_offset_idx = 0;// 0 means -gl

		// for all offsets
		for (int_t i = 0; i < num_ofs; i++)
		{
			// Algorithm 1 and Algorithm 2 of the paper and the masking process
			cmp_seq = Q;
			cmp_seq.Compare(T);
			cmp_seq |= mask_seq;
			cmp_seq.Mask_Short_MEM(sl);
			cmp_seq.Mark_Edge();
			cmp_seq.Extract_MEM(mems, num_mems, current_offset, current_offset_idx, TM);
			if (num_mems > TM) // skip this to Smith-Waterman
			{
				bypass_SSW_by_TM++;
				return;
			}
			current_offset++;
			current_offset_idx++;

			T.Shift_Left_Logical(BIT_PER_BASE);
			mask_seq.Shift_Left_Insert_1(BIT_PER_BASE);
		}
	}
	
	// sorting MEM in a counting sort (see suplementary data of the paper) look for counting sort
	inline void Sort_MEM()
	{
		// clear sort bucket
		for (int_t i = 0; i < Q.seq_len; i++)
			csort_cummulative_cnt[i] = csort_temp[i] = 0;

		// fill sort bucket
		for (int_t i = 0; i < num_mems; i++)
			csort_cummulative_cnt[mems[i].EQ]++;

		// Cummulative item count
		int_t cnt = 0;
		for (int_t i = 0; i < Q.seq_len; i++)
			if (csort_cummulative_cnt[i])
			{
				int_t temp = csort_cummulative_cnt[i];
				csort_cummulative_cnt[i] = cnt;
				cnt += temp;
			}

		// sort loop
		for (int_t i = 0; i < num_mems; i++)
		{
			int_t end = mems[i].EQ;
			int_t idx = csort_cummulative_cnt[end] + csort_temp[end];
			csort_temp[end]++;
			smems[idx] = mems[i];
		}
	}
	
	inline void Quick_Sort_MEM()
	{
		Quick_Sort_MEM(0, num_mems - 1);
		for (int_t i = 0; i < num_mems; i++)
			smems[i] = mems[i];
	}

	inline void Quick_Sort_MEM(int_t left, int_t right)
	{
		int mid = (left + right) / 2;

		int i = left;
		int j = right;
		int pivot = mems[mid].EQ;

		while (left<j || i<right)
		{
			while (mems[i].EQ < pivot)
				i++;
			while (mems[j].EQ > pivot)
				j--;

			if (i <= j)
			{
				if (mems[i].EQ != mems[j].EQ)
				{
					MEM temp = mems[i];
					mems[i] = mems[j];
					mems[j] = temp;
				}
				i++;
				j--;
			}
			else
			{
				if (left<j)
					Quick_Sort_MEM(left, j);
				if (i<right)
					Quick_Sort_MEM(i, right);
				return;
			}
		}
	}

	// This is the main body of MEM-Align with all optimisation
	// Given a sorted list of all MEMs it find the last MEM in the optimal alignment
	inline void Align_MEM()
	{
		oa_score = 0;
		oa_last_mem_idx = -1;

		//init ofs_flag
		for (int_t i = 0; i < num_ofs; i++)
			ofs_flag[i] = 0;
		// To compute Si for each MEM
		for (int_t i = 0; i < num_mems; i++)
		{
			// First assume the alignment is begin with this MEM
			smems[i].Si = smems[i].L * Rm;
			smems[i].Wi = i;
			// Computation of S^j_i
			// the direction of this loop efffect the performance it shold go backward; other wise ofs_flag does not work
			for (int_t j = i - 1; j >= 0; j--)
			{
				num_SJI_normal++;
				// check ofs_flag (optimisation in Section 6.1 of the paper)
				if (ofs_flag[smems[j].OFS] == i)
					continue; // if the flag is set avoid checking

							  // find mismatches and gaps. Return false if fully overlapped (optimisation in Section 6.1 of the paper)
				if (smems[i].Extend_To(smems[j]))
				{
					num_SJI_Omega++;
					if (smems[j].EQ < smems[i].BQ && smems[j].ET < smems[i].BT)
						ofs_flag[smems[j].OFS] = i; // set the flag for this offset to avoid checking other mem with same offset

					if (smems[i].mismatch_len < TD)
					{
						num_SJI_Omega_TD++;
						// compute score if extend Mi to Mj (in paper it is j to i)
						int_t Sij = smems[j].Si;
						if (smems[i].gap_len)
						{
							Sij -= smems[i].gap_len * Pe;
							Sij -= Po;
						}
						Sij += (smems[i].EQ - smems[i].XBQ + 1) * Rm;

						// before applying mismatch penalty first we check if ammend could change the result
						// Then we Amend and finally we apply mismatch penalty
						{
							// Amend mismatch penalty if length of mismatch is larger than one
							if (smems[i].mismatch_len > 1)
							{
								numNeededSeqCmp++;
								// check if the maximum recovered mismatch could leads to a new Maximum Score for the MEM
								// there should be at least one mismatch in every sl bases to be removed
								// for exampel if sl is 4 and mismatch len is 11 the XMMMXMMMXMM shwos the maximum matches M
								// and minumum mismatches X that results in elimination of whole region.
								// note that if there is no gap there should be a mismatch on the either side so XMMMXMMMXMX would be the acceptable pattern
								// however if there is a gap we could consider a mismatch only on one side
								int_t min_mismatch = smems[i].mismatch_len / sl;
								int_t mod = smems[i].mismatch_len % sl;
								if (mod)
									min_mismatch++;
								if (smems[i].gap_len == 0 && mod != 1)
									min_mismatch++;
								int_t max_match = smems[i].mismatch_len - min_mismatch;
								int_t temp_score = Sij + (max_match * Rm) - (min_mismatch * Px);

								int_t real_match = 0;
								if (temp_score > smems[i].Si)
								{
									numSeqCmpAfterPreCompute++;
									real_match = Amend_Mismatch_Len(smems[i], smems[j]);
								}
								// apply amended mismatch penalty
								Sij += real_match * Rm;
								Sij -= (smems[i].mismatch_len - real_match) * Px;
							}
							else
								Sij -= smems[i].mismatch_len * Px;
						}
						// check if the new score is maximum
						if (Sij > smems[i].Si)
						{
							smems[i].Si = Sij;
							smems[i].Wi = j;
						}
					}
					else // if (smems[i].mismatch_len >= D)
					{
						// do nothing
					}
				}
			}
			// check if S*i is larger than oa_score
			if (smems[i].Si > oa_score)
			{
				oa_score = smems[i].Si;
				oa_last_mem_idx = i;
			}
		}
	}

	inline int_t Amend_Mismatch_Len(MEM &Mi, MEM &Mj)
	{
		// number of match extending Mi and Mj
		int_t match_i = 0;
		int_t match_j = 0;

		// first consider extending Mi
		int_t iq = Mi.BQ - 1;
		int_t it = Mi.BT - 1;
		for (int_t i = 0; i < Mi.mismatch_len; i++)
		{
			if (Q.str[iq] == T.str[it])
				match_i++;
			iq--;
			it--;
		}

		if (Mi.gap_len == 0)
			return match_i;

		// then consider extending Mj
		iq = Mj.EQ + 1;
		it = Mj.ET + 1;
		for (int_t i = 0; i < Mi.mismatch_len; i++)
		{
			if (Q.str[iq] == T.str[it])
				match_j++;
			iq++;
			it++;
		}

		if (match_i > match_j)
		{
			Mi.amended_match_stick_to_pervious_MEM = false;
			return match_i;
		}
		else
		{
			Mi.amended_match_stick_to_pervious_MEM = true;
			return match_j;
		}
	}

	inline void Back_Track()
	{
		if (oa_last_mem_idx == -1)
		{
			oa_first_mem_idx = -1;
			return;
		}
		//not_reliable = false;

		int_t current = oa_last_mem_idx;

		oa_num_mems = 1;

		smems[current].WRi = current; // close the end of alignment

		// Backtracking
		while (smems[current].Wi != current)
		{
			//added for max mismatch len
			smems[current].Extend_To(smems[smems[current].Wi]);
			//if (smems[current].mismatch_len > mmlioa)
			//{
			//	oa_first_mem_idx = -1;
			//	not_reliable = true;
			//	return;
			//}


			oa_num_mems++;
			smems[smems[current].Wi].WRi = current;
			current = smems[current].Wi;
		}

		oa_first_mem_idx = current;

		oa_query_start_pos = smems[oa_first_mem_idx].BQ;
		oa_target_start_pos = smems[oa_first_mem_idx].BT;
		oa_query_end_pos = smems[oa_last_mem_idx].EQ;
		oa_target_end_pos = smems[oa_last_mem_idx].ET;

		//if(oa_query_start_pos>mmlioa || oa_query_end_pos<(Q.seq_len-mmlioa))
		//{
		//	oa_first_mem_idx = -1;
		//	not_reliable = true;
		//	return;
		//}
	}

	inline void Head_Tail_Extend()
	{
		int_t max = 0;
		int_t max_it = 0;
		int_t max_iq = 0;

		// Head Extend
		if ((smems[oa_first_mem_idx].BQ > 1) && (smems[oa_first_mem_idx].BT > 1))
		{
			int_t iq = smems[oa_first_mem_idx].BQ - 1;
			int_t it = smems[oa_first_mem_idx].BT - 1;
			int_t temp_score = 0;
			while(it>=0 && iq>=0)
			{
				if (Q.str[iq] == T.str[it])
					temp_score += Rm;
				else
					temp_score -= Px;
				if (temp_score > max)
				{
					max = temp_score;
					max_it = it;
					max_iq = iq;
				}
				iq--;
				it--;
			}
		}
		if (max > 0)
		{
			oa_score += max;
			oa_target_start_pos = max_it;
			oa_query_start_pos = max_iq;
		}

		// Tail Extend
		max = 0;
		if ((smems[oa_last_mem_idx].ET < (T.seq_len - 2)) && (smems[oa_last_mem_idx].EQ < (Q.seq_len - 2)))
		{
			int_t iq = smems[oa_last_mem_idx].EQ + 1;
			int_t it = smems[oa_last_mem_idx].ET + 1;
			int_t temp_score = 0;
			while (it < T.seq_len && iq < Q.seq_len)
			{
				if (Q.str[iq] == T.str[it])
					temp_score += Rm;
				else
					temp_score -= Px;
				if (temp_score > max)
				{
					max = temp_score;
					max_it = it;
					max_iq = iq;
				}
				iq++;
				it++;
			}
		}
		if (max > 0)
		{
			oa_score += max;
			oa_target_end_pos = max_it;
			oa_query_end_pos = max_iq;
		}
	}

	inline void Align(char *t, char *q)
	{
		T.Set(t);
		//printf("\n\n>%s\n>", t);
		//target.Print_Seq();
		Q.Set(q);
		//printf("\n\n>%s\n>", q);
		//query.Print_Seq();
		add_ticks(ET_CONVERT);
		Extract_MEM();
		add_ticks(ET_EXTRACT);
		oa_score = -1;
		if ((num_mems <= TM) && (num_mems != 0)) // (num_mems>1) to be checked : might result in better accuracy
		{
			Sort_MEM();
			add_ticks(ET_SORT);
			Align_MEM();
			add_ticks(ET_ALN);
			Back_Track();
			Head_Tail_Extend();
			add_ticks(ET_BACKTRACK);
			if (oa_score <= TS)
			{
				oa_score = SSW->SSW(q, t, len, Po + Pe, Pe);
				add_ticks(ET_SW);
				processed_by = SMITH_WATERMAN;
				bypass_SSW_total++;
				bypass_SSW_by_TS++;
			}
			else
			{
				processed_by = MEM_ALIGN;
			}
		}
		else
		{
			oa_score = SSW->SSW(q, t, len, Po+Pe, Pe);
			add_ticks(ET_SW);
			processed_by = SMITH_WATERMAN;
			bypass_SSW_total++;
		}
	}

	inline void Align_CS_A(char *t, char *q)
	{
		T.Set(t);
		//printf("\n\n>%s\n>", t);
		//target.Print_Seq();
		Q.Set(q);
		//printf("\n\n>%s\n>", q);
		//query.Print_Seq();
		add_ticks(ET_CONVERT);
		Extract_MEM();
		add_ticks(ET_EXTRACT);
		oa_score = -1;
		if (/*(num_mems <= TM) &&*/ (num_mems != 0)) // (num_mems>1) to be checked : might result in better accuracy
		{
			Sort_MEM();
			add_ticks(ET_SORT);
			//Align_MEM();
			//add_ticks(ET_ALN);
			//Back_Track();
			//Head_Tail_Extend();
			//add_ticks(ET_BACKTRACK);
			//if (oa_score <= TS)
			//{
			//	oa_score = SSW->SSW(q, t, len, Po + Pe, Pe);
			//	add_ticks(ET_SW);
			//	processed_by = SMITH_WATERMAN;
			//	bypass_SSW_total++;
			//	bypass_SSW_by_TS++;
			//}
			//else
			//{
			//	processed_by = MEM_ALIGN;
			//}
		}
		else
		{
			//oa_score = SSW->SSW(q, t, len, Po + Pe, Pe);
			//add_ticks(ET_SW);
			//processed_by = SMITH_WATERMAN;
			//bypass_SSW_total++;
		}
	}

	inline void Align_CS_B(char *t, char *q)
	{
		T.Set_Loop(t);
		//printf("\n\n>%s\n>", t);
		//target.Print_Seq();
		Q.Set_Loop(q);
		//printf("\n\n>%s\n>", q);
		//query.Print_Seq();
		add_ticks(ET_CONVERT);
		Extract_MEM();
		add_ticks(ET_EXTRACT);
		oa_score = -1;
		if (/*(num_mems <= TM) &&*/ (num_mems != 0)) // (num_mems>1) to be checked : might result in better accuracy
		{
			Quick_Sort_MEM();
			add_ticks(ET_SORT);
			//Align_MEM();
			//add_ticks(ET_ALN);
			//Back_Track();
			//Head_Tail_Extend();
			//add_ticks(ET_BACKTRACK);
			//if (oa_score <= TS)
			//{
			//	oa_score = SSW->SSW(q, t, len, Po + Pe, Pe);
			//	add_ticks(ET_SW);
			//	processed_by = SMITH_WATERMAN;
			//	bypass_SSW_total++;
			//	bypass_SSW_by_TS++;
			//}
			//else
			//{
			//	processed_by = MEM_ALIGN;
			//}
		}
		else
		{
			//oa_score = SSW->SSW(q, t, len, Po + Pe, Pe);
			//add_ticks(ET_SW);
			//processed_by = SMITH_WATERMAN;
			//bypass_SSW_total++;
		}
	}

	inline void Align_ExtractMEM(char *t, char *q)
	{
		T.Set(t);
		Q.Set(q);
		Extract_MEM();
	}

	inline void Align_Convert(char *t, char *q)
	{
		T.Set(t);
		Q.Set(q);
	}

	inline void Align_Convert_Loop(char *t, char *q)
	{
		T.str = t;
		T.Convert_Loop();
		Q.str = q;
		Q.Convert_Loop();
	}

	inline void Num_Possible_EMs_From_MEMs()
	{
		num_mems = 0;
		for (int_t i = 0; i < num_mems; i++)
		{
			for (int_t j = 1; j <= mems[i].L; j++)
				num_mems += mems[i].L - (j - 1);
		}
	}

	inline void Generate_Alignment()
	{
		if (processed_by == SMITH_WATERMAN)
		{
			Printf("\nThis alignment is generated by SSW\n");
			return;
		}
		int_t cur, next;
		smems[oa_first_mem_idx].XBQ = smems[oa_first_mem_idx].BQ;
		smems[oa_first_mem_idx].XBT = smems[oa_first_mem_idx].BT;
		// fix gaps and mismatches in optimal alignment path
		if (oa_num_mems > 1)
		{
			cur = oa_first_mem_idx;
			next = smems[cur].WRi;
			while (next != cur)
			{
				smems[next].Extend_To(smems[cur]);
				//smems[cur].Print();
				if (smems[next].mismatch_len > 1 && smems[next].gap_len)
					Amend_Mismatch_Len(smems[next], smems[cur]);
				cur = next;
				next = smems[cur].WRi;
			}
			//smems[cur].Print();
			//Printf("\n");
		}
		else
		{
			//smems[oa_first_mem_idx].Print();
			//Printf("\n");
		}

		Print_Target();
		printf("\n");
		Print_Query();
		printf("\n");
		Print_Cigar();
		printf("\nOptimal Alignment Score: %d\nCigar String: %s\n\n", oa_score, cstr.str);
	}

	inline void Print_Target()
	{
		int_t cur, next;
		int_t color = 0;
		// print initial gap
		if (oa_target_start_pos < oa_query_start_pos)
			for (int_t i = oa_target_start_pos; i < oa_query_start_pos; i++)
				printf(" ");

		// print initial clipped region
		if (oa_target_start_pos > 0)
			for (int_t i = 0; i < oa_target_start_pos; i++)
				Print_C(T.str[i]);

		// print initial recovered region
		if (oa_target_start_pos < smems[oa_first_mem_idx].XBT)
		{
			int_t j = oa_query_start_pos;
			for (int_t i = oa_target_start_pos; i < smems[oa_first_mem_idx].XBT; i++)
			{
				if (T.str[i] == Q.str[j])
					Print_RM(T.str[i])
				else
					Print_X(T.str[i]);
				j++;
			}
		}

		// print all mem and their distance except last MEM
		cur = oa_first_mem_idx;
		next = smems[cur].WRi;
		while (cur != next)
		{
			// print current MEM
			for (int_t i = smems[cur].XBT; i <= smems[cur].ET; i++)
				Print_M(T.str[i], color);
			color++;

			// print distance
			if (smems[next].mismatch_len > 0)
			{
				if (smems[next].gap_len > 0)
				{
					// first print gaps then print mismatch
					if (smems[next].amended_match_stick_to_pervious_MEM == false)
					{
						if (smems[next].insertion_len > 0) // if insertion print dash (-)
						{
							for (int_t i = 0; i < smems[next].gap_len; i++)
								Print_G("-");
						}
						else // if deletion print bases
						{
							int_t j = smems[cur].ET + 1;
							for (int_t i = 0; i < smems[next].gap_len; i++)
							{
								Print_G(T.str[j]);
								j++;
							}
						}


						int_t j = smems[next].XBQ - smems[next].mismatch_len;
						for (int_t i = (smems[next].XBT - smems[next].mismatch_len); i < smems[next].XBT; i++)
						{
							if (T.str[i] == Q.str[j])
								Print_RM(T.str[i])
							else
								Print_X(T.str[i]);
							j++;
						}

					}
					// first print mismatch, then print gaps
					else
					{
						int_t j = smems[cur].EQ + 1;
						for (int_t i = smems[cur].ET + 1; i <= smems[cur].ET + smems[next].mismatch_len; i++)
						{
							if (T.str[i] == Q.str[j])
								Print_RM(T.str[i])
							else
								Print_X(T.str[i]); 
							j++;
						}

						if (smems[next].insertion_len > 0) // if insertion print dash (-)
						{
							for (int_t i = 0; i < smems[next].gap_len; i++)
								Print_G("-");
						}
						else // if deletion print bases
						{
							int_t j = smems[next].XBT - smems[next].mismatch_len;
							for (int_t i = 0; i < smems[next].gap_len; i++)
							{
								Print_G(T.str[j]);
								j++;
							}
						}
					}

				}
				else // gap_len == 0
				{
					int_t j = smems[cur].EQ + 1;
					for (int_t i = smems[cur].ET + 1; i < smems[next].XBT; i++)
					{
						if (T.str[i] == Q.str[j])
							Print_RM(T.str[i])
						else
							Print_X(T.str[i]); 
						j++;
					}
				}
			}
			else if (smems[next].gap_len > 0) // mismatch_len==0
			{
				if (smems[next].insertion_len > 0) // if insertion print dash (-)
				{
					for (int_t i = 0; i < smems[next].gap_len; i++)
						Print_G("-");
				}
				else // if deletion print bases
				{
					int_t j = smems[cur].ET + 1;
					for (int_t i = 0; i < smems[next].gap_len; i++)
					{
						Print_G(T.str[j]);
						j++;
					}
				}
			}

			cur = next;
			next = smems[cur].WRi;
		}

		// print last mem
		for (int_t i = smems[oa_last_mem_idx].XBT; i <= smems[oa_last_mem_idx].ET; i++)
			Print_M(T.str[i],color);
		color++;

		// print final recovered region
		if (oa_target_end_pos > smems[oa_last_mem_idx].ET)
		{
			int_t j = smems[oa_last_mem_idx].EQ + 1;
			for (int_t i = smems[oa_last_mem_idx].ET + 1; i <= oa_target_end_pos; i++)
			{
				if (T.str[i] == Q.str[j])
					Print_RM(T.str[i])
				else
					Print_X(T.str[i]);
				j++;
			}
		}

		// print final clipped region
		if (oa_target_end_pos < T.seq_len - 1)
			for (int_t i = oa_target_end_pos + 1; i < T.seq_len; i++)
				Print_C(T.str[i]);

		return;
	}

	inline void Print_Query()
	{
		int_t cur, next;
		int_t color = 0;
		// print initial gap
		if (oa_target_start_pos > oa_query_start_pos)
			for (int_t i = oa_query_start_pos; i < oa_target_start_pos; i++)
				printf(" ");

		// print initial clipped region
		if (oa_query_start_pos > 0)
			for (int_t i = 0; i < oa_query_start_pos; i++)
				Print_C(Q.str[i]);

		// print initial recovered region
		if (oa_query_start_pos < smems[oa_first_mem_idx].XBQ)
		{
			int_t j = oa_target_start_pos;
			for (int_t i = oa_query_start_pos; i < smems[oa_first_mem_idx].XBQ; i++)
			{
				if (T.str[j] == Q.str[i])
					Print_RM(Q.str[i])
				else
					Print_X(Q.str[i]);
				j++;
			}
		}

		// print all mem and their distance except last MEM
		cur = oa_first_mem_idx;
		next = smems[cur].WRi;
		while (cur != next)
		{
			// print current MEM
			for (int_t i = smems[cur].XBQ; i <= smems[cur].EQ; i++)
				Print_M(Q.str[i], color);
			color++;

			// print distance
			if (smems[next].mismatch_len > 0)
			{
				if (smems[next].gap_len > 0)
				{
					// first print gaps then mismatch
					if (smems[next].amended_match_stick_to_pervious_MEM == false)
					{
						if (smems[next].deletion_len > 0) // if deletion print dash (-)
						{
							for (int_t i = 0; i < smems[next].gap_len; i++)
								Print_G("-");
						}
						else // if insertion print bases
						{
							int_t j = smems[cur].EQ + 1;
							for (int_t i = 0; i < smems[next].gap_len; i++)
							{
								Print_G(Q.str[j]);
								j++;
							}
						}

						int_t j = smems[next].XBT - smems[next].mismatch_len;
						for (int_t i = (smems[next].XBQ - smems[next].mismatch_len); i < smems[next].XBQ; i++)
						{
							if (T.str[j] == Q.str[i])
								Print_RM(Q.str[i])
							else
								Print_X(Q.str[i]);
							j++;
						}

					}
					// first print mismatch then print gaps
					else
					{
						int_t j = smems[cur].ET + 1;
						for (int_t i = smems[cur].EQ + 1; i < smems[cur].EQ + 1 + smems[next].mismatch_len; i++)
						{
							if (T.str[j] == Q.str[i])
								Print_RM(Q.str[i])
							else
								Print_X(Q.str[i]);
							j++;
						}

						if (smems[next].deletion_len > 0) // if deletion print dash (-)
						{
							for (int_t i = 0; i < smems[next].gap_len; i++)
								Print_G("-");
						}
						else // if insertion print bases
						{
							int_t j = smems[next].XBQ - smems[next].mismatch_len;
							for (int_t i = 0; i < smems[next].gap_len; i++)
							{
								Print_G(Q.str[j]);
								j++;
							}
						}
					}

				}
				else // insertion_len == 0
				{
					int_t j = smems[cur].ET + 1;
					for (int_t i = smems[cur].EQ + 1; i < smems[next].XBQ; i++)
					{
						if (T.str[j] == Q.str[i])
							Print_RM(Q.str[i])
						else
							Print_X(Q.str[i]);
						j++;
					}
				}
			}
			else if (smems[next].gap_len > 0) // mismatch_len==0
			{
				if (smems[next].deletion_len > 0) // if deletion print dash (-)
				{
					for (int_t i = 0; i < smems[next].gap_len; i++)
						Print_G("-");
				}
				else // if insertion print bases
				{
					int_t j = smems[cur].EQ + 1;
					for (int_t i = 0; i < smems[next].gap_len; i++)
					{
						Print_G(Q.str[j]);
						j++;
					}
				}
			}

			cur = next;
			next = smems[cur].WRi;
		}

		// print last mem
		for (int_t i = smems[oa_last_mem_idx].XBQ; i <= smems[oa_last_mem_idx].EQ; i++)
			Print_M(Q.str[i], color);
		color++;
		// print final recovered region
		if (oa_query_end_pos > smems[oa_last_mem_idx].EQ)
		{
			int_t j = smems[oa_last_mem_idx].ET + 1;
			for (int_t i = smems[oa_last_mem_idx].EQ + 1; i <= oa_query_end_pos; i++)
			{
				if (T.str[j] == Q.str[i])
					Print_RM(Q.str[i])
				else
					Print_X(Q.str[i]);
				j++;
			}
		}

		// print final clipped region
		if (oa_query_end_pos < Q.seq_len - 1)
			for (int_t i = oa_query_end_pos + 1; i < Q.seq_len; i++)
				Print_C(Q.str[i]);

		return;
	}

	inline void Print_Cigar()
	{
		int_t cur, next;
		int_t color = 0;
		// print initial gap
		if (oa_target_start_pos > oa_query_start_pos)
			for (int_t i = oa_query_start_pos; i < oa_target_start_pos; i++)
				printf(" ");

		// print initial clipped region
		if (oa_query_start_pos > 0)
			for (int_t i = 0; i < oa_query_start_pos; i++)
			{
				Print_C("S"); cstr.Add('S');
			}

		// print initial recovered region
		if (oa_query_start_pos < smems[oa_first_mem_idx].XBQ)
		{
			int_t j = oa_target_start_pos;
			for (int_t i = oa_query_start_pos; i < smems[oa_first_mem_idx].XBQ; i++)
			{
				if (T.str[j] == Q.str[i])
				{
					Print_RM("M"); cstr.Add('M');
				}
				else
				{
					Print_X("X"); cstr.Add('X');
				}
				j++;
			}
		}

		// print all mem and their distance except last MEM
		cur = oa_first_mem_idx;
		next = smems[cur].WRi;
		while (cur != next)
		{
			// print current MEM
			for (int_t i = smems[cur].XBQ; i <= smems[cur].EQ; i++)
			{
				Print_M("M", color); cstr.Add('M');
			}
			color++;

			// print distance
			if (smems[next].mismatch_len > 0)
			{
				if (smems[next].gap_len > 0)
				{
					// first print gaps then mismatch
					if (smems[next].amended_match_stick_to_pervious_MEM == false)
					{
						if (smems[next].deletion_len > 0) // if deletion print dash (-)
						{
							for (int_t i = 0; i < smems[next].gap_len; i++)
							{
								Print_G("D"); cstr.Add('D');
							}
						}
						else // if insertion print bases
						{
							int_t j = smems[cur].EQ + 1;
							for (int_t i = 0; i < smems[next].gap_len; i++)
							{
								Print_G("I"); cstr.Add('I');
								j++;
							}
						}

						int_t j = smems[next].XBT - smems[next].mismatch_len;
						for (int_t i = (smems[next].XBQ - smems[next].mismatch_len); i < smems[next].XBQ; i++)
						{
							if (T.str[j] == Q.str[i])
							{
								Print_RM("M"); cstr.Add('M');
							}
							else
							{
								Print_X("X"); cstr.Add('X');
							}
							j++;
						}

					}
					// first print mismatch then print gaps
					else
					{
						int_t j = smems[cur].ET + 1;
						for (int_t i = smems[cur].EQ + 1; i < smems[cur].EQ + 1 + smems[next].mismatch_len; i++)
						{
							if (T.str[j] == Q.str[i])
							{
								Print_RM("M"); cstr.Add('M');
							}
							else
							{
								Print_X("X"); cstr.Add('X');
							}
							j++;
						}

						if (smems[next].deletion_len > 0) // if deletion print dash (-)
						{
							for (int_t i = 0; i < smems[next].gap_len; i++)
							{
								Print_G("D"); cstr.Add('D');
							}
						}
						else // if insertion print bases
						{
							int_t j = smems[next].XBQ - smems[next].mismatch_len;
							for (int_t i = 0; i < smems[next].gap_len; i++)
							{
								Print_G("I"); cstr.Add('I');
								j++;
							}
						}
					}

				}
				else // insertion_len == 0
				{
					int_t j = smems[cur].ET + 1;
					for (int_t i = smems[cur].EQ + 1; i < smems[next].XBQ; i++)
					{
						if (T.str[j] == Q.str[i])
						{
							Print_RM("M"); cstr.Add('M');
						}
						else
						{
							Print_X("X"); cstr.Add('X');
						}
						j++;
					}
				}
			}
			else if (smems[next].gap_len > 0) // mismatch_len==0
			{
				if (smems[next].deletion_len > 0) // if deletion print dash (-)
				{
					for (int_t i = 0; i < smems[next].gap_len; i++)
					{
						Print_G("D"); cstr.Add('D');
					}
				}
				else // if insertion print bases
				{
					int_t j = smems[cur].EQ + 1;
					for (int_t i = 0; i < smems[next].gap_len; i++)
					{
						Print_G("I"); cstr.Add('I');
						j++;
					}
				}
			}

			cur = next;
			next = smems[cur].WRi;
		}

		// print last mem
		for (int_t i = smems[oa_last_mem_idx].XBQ; i <= smems[oa_last_mem_idx].EQ; i++)
		{
			Print_M("M", color);  cstr.Add('M');
		}
		color++;
		// print final recovered region
		if (oa_query_end_pos > smems[oa_last_mem_idx].EQ)
		{
			int_t j = smems[oa_last_mem_idx].ET + 1;
			for (int_t i = smems[oa_last_mem_idx].EQ + 1; i <= oa_query_end_pos; i++)
			{
				if (T.str[j] == Q.str[i])
				{
					Print_RM("M"); cstr.Add('M');
				}
				else
				{
					Print_X("X"); cstr.Add('X');
				}
				j++;
			}
		}

		// print final clipped region
		if (oa_query_end_pos < Q.seq_len - 1)
			for (int_t i = oa_query_end_pos + 1; i < Q.seq_len; i++)
			{
				Print_C("S"); cstr.Add('S');
			}

		cstr.Close();
		return;
	}

};

