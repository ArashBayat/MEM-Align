#pragma once

#include "define.h"
#include "mem.h"
class SEQUENCE;

extern SEQUENCE temp_seq[2];

class SEQUENCE
{
public:
	char *str;					// sequence data in ASCII
	word_t *bit_vector;			// sequence data in bitvector
								// This also accomodate pading that allow the sequence to be shifted to right and left with out losing information
	int_t seq_len;				// number of bases in the sequence

	int_t bit_last_word;		// number of bit in last word
	int_t base_last_word;		// number of base in last word
	int_t num_word;				// number of word including the last one
	int_t num_filed_word;		// number of word excluding the last one
	int_t total_word;			// total word inlcuding pads

	word_t mask_last_word_0;	// mask remaing bit with 0
	word_t mask_last_word_1;	// mask remaing bit with 1

	int_t pad_word;				// number of padded word on the other side
	int_t sw_idx;				// start word idx (excluding pad words)
	int_t ew_idx;				// end word idx (excluding pad words)

	SEQUENCE()
	{
	};
	
	~SEQUENCE()
	{
	};

	inline void Init(int_t len, int_t p_pad_word)
	{
		seq_len = len;
		pad_word = p_pad_word;
		num_filed_word = seq_len / BASE_PER_WORD;
		num_word = num_filed_word;
		base_last_word = seq_len % BASE_PER_WORD;
		bit_last_word = base_last_word * BIT_PER_BASE;
		mask_last_word_1 = 0;
		if (base_last_word != 0)
		{
			num_word++;
			mask_last_word_1 = POWER2(((BASE_PER_WORD - base_last_word)*BIT_PER_BASE)) - 1;
		}
		mask_last_word_0 = ~mask_last_word_1;
		
		// len include number of word needed for convert and pad_words
		int_t t_len = seq_len / BYTE_PER_WORD;
		if (seq_len % BYTE_PER_WORD)
			t_len++;
		t_len += pad_word * 2;

		bit_vector = 0;
		bit_vector = new word_t[t_len];

		sw_idx = pad_word;
		ew_idx = sw_idx + num_word - 1;
		total_word = num_word + (pad_word * 2);
	}

	inline word_t Get_Word(int_t i)
	{
		return bit_vector[i + sw_idx];
	}

	inline void Mask_1() { bit_vector[ew_idx] |= mask_last_word_1; }

	inline void Mask_0() { bit_vector[ew_idx] &= mask_last_word_0; }

	inline void Fill_Match() // fill the sequence all with matches (00) and remaining bits of last word and padding all with mismatch (11)
	{
		for (int_t i = sw_idx; i <= ew_idx; i++)
			bit_vector[i] = 0;
		Fill_Pad();
	}

	inline void Fill_Pad()  // fill the remaining bits of last word and padding all with mismatch (11)
	{
		Mask_1();
		for (int_t i = 0; i < sw_idx; i++)
			bit_vector[i] = ALL_F;

		for (int_t i = ew_idx + 1; i < total_word; i++)
			bit_vector[i] = ALL_F;
	}

	inline void Convert() // Bit level parallelized convert method
	{
		memcpy(&bit_vector[sw_idx], str, seq_len);

		// convert 8 base in each word first
		for (int_t i = sw_idx; i < (sw_idx + (num_word * 4)); i++)
		{
			bit_vector[i] &= ALL_06;
			bit_vector[i] >>= 1;
			bit_vector[i] |= (bit_vector[i] << 10) | (bit_vector[i] << 20) | (bit_vector[i] << 30);
			bit_vector[i] &= 0xFF000000FF000000;
			bit_vector[i] = (bit_vector[i] >> 8) | (bit_vector[i] << 32);
			bit_vector[i] &= 0xFFFF000000000000;
		}

		// merge every 4 words together
		for (int_t i = 0; i < num_word; i++)
		{
			int_t idx = (i * 4) + pad_word;
			bit_vector[i + pad_word] = bit_vector[idx] | bit_vector[idx + 1] >> 16 | bit_vector[idx + 2] >> 32 | bit_vector[idx + 3] >> 48;
		}
		Fill_Pad();
	}

	inline void Convert_Loop() // traditional loop based convert method (not to be used)
	{
		int_t k = 0;
		for (int_t i = sw_idx; i < ew_idx; i++)
		{
			for (int_t j = 0; j < BASE_PER_WORD; j++)
			{
				bit_vector[i] <<= BIT_PER_BASE;
				switch (str[k])
				{
				case 'A':
				case 'a':
					break;
				case 'C':
				case'c':
					bit_vector[i] |= 1;
					break;
				case 'T':
				case't':
					bit_vector[i] |= 2;
					break;
				case 'G':
				case 'g':
					bit_vector[i] |= 3;
					break;
				default:
					bit_vector[i] |= 0;
				}
				k++;
			}
		}
		// converting last word
		if (base_last_word)
		{
			for (int_t j = 0; j < base_last_word; j++)
			{
				bit_vector[ew_idx] <<= BIT_PER_BASE;
				switch (str[k])
				{
				case 'A':
				case 'a':
					break;
				case 'C':
				case'c':
					bit_vector[ew_idx] |= 1;
					break;
				case 'T':
				case't':
					bit_vector[ew_idx] |= 2;
					break;
				case 'G':
				case 'g':
					bit_vector[ew_idx] |= 3;
					break;
				default:
					bit_vector[ew_idx] |= 0;
				}
				k++;
			}
			bit_vector[ew_idx] <<= BIT_PER_WORD - bit_last_word;
		}
		Fill_Pad();
	}

	inline void Set(char *seq_str) // the Init function must be called before. does not copy the string. only takes pointer and does the convert
	{
		str = seq_str;
		Convert();
	}

	inline void Set_Loop(char *seq_str) // the Init function must be called before. does not copy the string. only takes pointer and does the convert
	{
		str = seq_str;
		Convert_Loop();
	}

	inline void Shift_Left_Logical(int_t shift) // less than the size of a word
	{
		if (shift == 0) return;

		int_t shift_inverse = BIT_PER_WORD - shift;
		for (int_t i = 0; i < total_word - 1; i++)
		{
			bit_vector[i] <<= shift;
			bit_vector[i] |= bit_vector[i + 1] >> shift_inverse;
		}
		bit_vector[total_word - 1] <<= shift;
	}

	inline void Shift_Right_Logical(int_t shift) // less than the size of a word
	{
		if (shift == 0) return;

		int_t shift_inverse = BIT_PER_WORD - shift;
		for (int_t i = total_word - 1; i > 0; i--)
		{
			bit_vector[i] >>= shift;
			bit_vector[i] |= bit_vector[i - 1] << shift_inverse;
		}
		bit_vector[0] >>= shift;
	}

	inline void Long_Shift_Left_Logical(int_t shift) // more than the size of a word
	{
		if (shift == 0) return;

		int_t word_shift = shift >> BIT_PER_WORD_LOG;

		if (word_shift)
		{
			for (int_t i = 0; i < total_word - word_shift; i++)
				bit_vector[i] = bit_vector[i + word_shift];
			for (int_t i = total_word - word_shift; i < num_word; i++)
				bit_vector[i] = 0;
		}
		else
		{
			Shift_Left_Logical(shift);
			return;
		}

		int_t word_ofs = shift & BIT_PER_WORD_MASK;
		if (word_ofs)
			Shift_Left_Logical(word_ofs);
	}

	inline void Long_Shift_Right_Logical(int_t shift) // more than the size of a word
	{
		if (shift == 0) return;

		int_t word_shift = shift >> BIT_PER_WORD_LOG;

		if (word_shift)
		{
			for (int_t i = total_word - 1; i >= word_shift; i--)
				bit_vector[i] = bit_vector[i - word_shift];
			for (int_t i = 0; i < word_shift; i++)
				bit_vector[i] = 0;
		}
		else
		{
			Shift_Right_Logical(shift);
			return;
		}

		int_t word_ofs = shift & BIT_PER_WORD_MASK;
		if (word_ofs)
			Shift_Right_Logical(word_ofs);
	}

	inline void Shift_Left_Insert_1(int_t shift) // Shift left and insert (11...) 
	{
		if (shift == 0) return;

		int_t shift_inverse = BIT_PER_WORD - shift;
		for (int_t i = 0; i < total_word - 1; i++)
		{
			bit_vector[i] <<= shift;
			bit_vector[i] |= bit_vector[i + 1] >> shift_inverse;
		}
		bit_vector[total_word - 1] <<= shift;
		bit_vector[total_word - 1] |= POWER2(shift) - 1;
	}

	inline void Shift_Right_Insert_1(int_t shift) // Shift left and insert (11...) 
	{
		if (shift == 0) return;

		int_t shift_inverse = BIT_PER_WORD - shift;
		for (int_t i = total_word - 1; i > 0; i--)
		{
			bit_vector[i] >>= shift;
			bit_vector[i] |= bit_vector[i - 1] << shift_inverse;
		}
		bit_vector[0] >>= shift;
		bit_vector[0] |= (POWER2(shift) - 1) << (BIT_PER_WORD - shift);
	}

	inline void Long_Shift_Left_Insert_1(int_t shift) // Shift left and insert (11...) 
	{
		if (shift == 0) return;

		int_t word_shift = shift >> BIT_PER_WORD_LOG;

		if (word_shift)
		{
			for (int_t i = 0; i < total_word - word_shift; i++)
				bit_vector[i] = bit_vector[i + word_shift];
			for (int_t i = total_word - word_shift; i < num_word; i++)
				bit_vector[i] = ALL_F;
		}
		else
		{
			Shift_Left_Insert_1(shift);
			return;
		}

		int_t word_ofs = shift & BIT_PER_WORD_MASK;
		if (word_ofs)
			Shift_Left_Insert_1(word_ofs);
	}

	inline void Long_Shift_Right_Insert_1(int_t shift) // Shift left and insert (11...) 
	{
		if (shift == 0) return;

		int_t word_shift = shift >> BIT_PER_WORD_LOG;

		if (word_shift)
		{
			for (int_t i = total_word - 1; i >= word_shift; i--)
				bit_vector[i] = bit_vector[i - word_shift];
			for (int_t i = 0; i < word_shift; i++)
				bit_vector[i] = ALL_F;
		}
		else
		{
			Shift_Right_Insert_1(shift);
			return;
		}

		int_t word_ofs = shift & BIT_PER_WORD_MASK;
		if (word_ofs)
			Shift_Right_Insert_1(word_ofs);
	}

	inline SEQUENCE &operator=(SEQUENCE &opr)  // the Init function must be called before
	{
		// only data is copied the rest of the fields are set in Init()
		// there is no need to copy str value
		for (int_t i = sw_idx; i <= ew_idx; i++)
			bit_vector[i] = opr.Get_Word(i);
		return (*this);
	}

	inline void operator<<=(const int_t shift)
	{
		Long_Shift_Left_Logical(shift);
	}

	inline void operator>>=(const int_t shift)
	{
		Long_Shift_Right_Logical(shift);
	}

	inline void operator|=(SEQUENCE &opr)
	{
		for (int_t i = sw_idx; i <= ew_idx; i++)
		{
			bit_vector[i] |= opr.Get_Word(i);
		}
	}

	inline void operator^=(SEQUENCE &opr)
	{
		for (int_t i = sw_idx; i <= ew_idx; i++)
		{
			bit_vector[i] ^= opr.Get_Word(i);
		}
	}

	inline void operator&=(SEQUENCE &opr)
	{
		for (int_t i = sw_idx; i <= ew_idx; i++)
		{
			bit_vector[i] &= opr.Get_Word(i);
		}
	}

	inline void Compare(SEQUENCE &opr)
	{
		for (int_t i = sw_idx; i <= ew_idx; i++)
		{
			bit_vector[i] ^= opr.Get_Word(i);
			bit_vector[i] |= bit_vector[i] >> 1;
			bit_vector[i] &= ALL_5;
			bit_vector[i] |= bit_vector[i] << 1;
		}
		Mask_1();
	}

	inline void Mask_Short_MEM(int_t threshold) // remove every mem shorter than threshold
	{
		temp_seq[0] = *this; // F in the paper
		temp_seq[1] = *this;
		for (int_t i = 1; i < threshold; i++)
		{
			temp_seq[1].Shift_Left_Insert_1(BIT_PER_BASE);
			temp_seq[0] |= temp_seq[1];
		}

		*this = temp_seq[0];
		for (int_t i = 1; i < threshold; i++)
		{
			temp_seq[0].Shift_Right_Insert_1(BIT_PER_BASE);
			*this &= temp_seq[0];
		}
	}

	inline void Mark_Edge()
	{
		temp_seq[0] = *this;
		temp_seq[0].Shift_Right_Insert_1(1);
		*this ^= temp_seq[0];
	}

	inline void Extract_MEM(MEM *mems, int_t &idx, int_t offset, int_t offset_idx, int_t TM)
	{
		int_t bit_ofs = 0;
		int_t one_cnt = 0;
		for (int_t i = sw_idx; i <= ew_idx; i++)
		{
			if (bit_vector[i] == 0)
			{
				bit_ofs += BIT_PER_WORD;
			}
			else
			{
				int_t bit_left = BIT_PER_WORD;
				word_t temp = bit_vector[i];
				bool first_one = true;
				while (1)
				{
					int_t lzcnt = LZCNT(temp);

					temp <<= lzcnt + 1;

					int_t to_add = lzcnt;
					if (!first_one)
						to_add++;
					first_one = false;

					bit_left -= to_add;
					bit_ofs += to_add;

					if ((one_cnt & 1) == 0)
					{
						mems[idx].BQ = bit_ofs / 2;
					}
					else
					{
						mems[idx].EQ = (bit_ofs / 2) - 1;
						mems[idx].BT = mems[idx].BQ + offset;
						mems[idx].ET = mems[idx].EQ + offset;
						mems[idx].OFS = offset_idx;
						mems[idx].L = mems[idx].EQ - mems[idx].BQ + 1;
						idx++;
						if (idx > TM) // skip this to Smith-Waterman
							return;
					}
					one_cnt++;

					if (temp == 0)
					{
						bit_ofs += bit_left;
						break;
					}
				}
			}
		}
		// close last mem if not closed ( this will happen when the sequence len is aligned ot M_WWORD)
		if ((one_cnt & 1) != 0)
		{
			mems[idx].EQ = seq_len - 1;
			mems[idx].BT = mems[idx].BQ + offset;
			mems[idx].ET = mems[idx].EQ + offset;
			mems[idx].OFS = offset;
			mems[idx].L = mems[idx].EQ - mems[idx].BQ + 1;
			idx++;
			if (idx > TM) // skip this to Smith-Waterman
				return;
		}
	}

/*	inline int_t Count_Mismatches()
	{
		Mask_0();
		int_t cnt = 0;
		for (int_t i = sw_idx; i <= ew_idx; i++)
		{
			word_t t = bit_vector[i] & ALL_5;
			cnt += POPCNT(t);
		}
		Mask_1();
		return cnt;
	}
*/

	inline void Print_Bit(char *alphabet)
	{
		int_t idx = 0;
		for (int_t i = sw_idx; i <= ew_idx; i++)
		{
			word_t p = bit_vector[i];
			for (int_t j = 0; j < BIT_PER_WORD; j++)
			{
				if (idx == (seq_len * BIT_PER_BASE))
					Printf("\t");
				switch (p >> (BIT_PER_WORD - 1))
				{
				case 0:
					Printf("%c", alphabet[0]);
					break;
				case 1:
					Printf("%c", alphabet[1]);
					break;
				default:
					Printf("%c", alphabet[2]);
				}
				p <<= 1;
				idx++;
			}
		}
		Printf("\n");
	}
	inline void Print(char *alphabet)
	{
		int_t idx = 0;
		for (int_t i = sw_idx; i <= ew_idx; i++)
		{
			word_t p = bit_vector[i];
			for (int_t j = 0; j < BASE_PER_WORD; j++)
			{
				if (idx == seq_len)
					Printf("\t");
				switch (p >> (BIT_PER_WORD - BIT_PER_BASE))
				{
				case 0:
					Printf("%c", alphabet[0]);
					break;
				case 1:
					Printf("%c", alphabet[1]);
					break;
				case 2:
					Printf("%c", alphabet[2]);
					break;
				case 3:
					Printf("%c", alphabet[3]);
					break;
				default:
					Printf("%c", alphabet[4]);
				}
				p <<= BIT_PER_BASE;
				idx++;
			}
		}
		Printf("\n");
	}
	inline void Print_Seq() { Print("ACTG "); }
	inline void Print_Compare() { Print("-   ~"); }
	inline void Print_Edge() { Print("--^-~"); }
	inline void Print_Hex()
	{
		Printf("\n");
		for (int_t i = sw_idx; i < ew_idx; i++)
			Printf("%llX", bit_vector[i]);
	}
	inline void Print_Hex_All()
	{
		Printf("\n");
		for (int_t i = 0; i < total_word; i++)
			Printf("%016llX", bit_vector[i]);
	}
};
