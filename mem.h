#pragma once

#include "define.h"

struct MEM
{
public:
	int_t BQ;		// Begin position in Q
	int_t EQ;		// End position in Q
	int_t BT;		// Begin position in T
	int_t ET;		// End Position in T
	int_t OFS;		// Ofset (BT-BQ)
	int_t L;		// lengh (ET-BT+1)

	int_t Si;		// Alignment score for optimal alignment ending at this MEM
	int_t Wi;		// pervious MEM in optimal alignment
	int_t WRi;		// next MEM in optimal alignment

	int_t OLQ;		// overlapped len in query
	int_t OLT;		// overlapped len in target
	int_t MO;		// maximum overlap len MAX(OLQ, OLT)

	// These are computed in relation to previous MEM in optimal alignment
	int_t LQ;  // distance between this MEM and prevoious MEM in Q 
	int_t LT;  // distance between this MEM and prevoious MEM in T 
	int_t LD;  // LT-LQ
	int_t XBT; // BT excluding overlpse
	int_t XBQ; // BQ excluding overlapse
	
	// these are to compute penalty between this MEM and previous MEM 
	int_t gap_len;
	int_t insertion_len;
	int_t deletion_len;
	int_t mismatch_len;

	// if there is a gap between two MEM this parameter says to which MEM gap shoul be attached
	bool amended_match_stick_to_pervious_MEM;

	// find the number of matches mismatches and gaps between this MEM and M_j to computation of  P^j_i
	inline bool Extend_To(MEM &Mj)
	{
		// Check Mj ends and starts before Mi in both T and Q
		// Although sorted Mj.EQ can be equal to EQ
		if (Mj.BQ >= BQ || Mj.EQ >= EQ || Mj.BT >= BT || Mj.ET >= ET)
			return false;

		// find ovelapse length in T and Q
		OLQ = (Mj.EQ >= BQ) ? (Mj.EQ - BQ + 1) : 0;
		OLT = (Mj.ET >= BT) ? (Mj.ET - BT + 1) : 0;

		// find Max overlaped to be excluded from current MEM
		MO = (OLQ > OLT) ? OLQ : OLT;

		// exclude overlapse BT and BQ
		XBT = BT + MO;
		XBQ = BQ + MO;

		// find distance in Q and T
		LQ = XBQ - Mj.EQ - 1; //
		LT = XBT - Mj.ET - 1; //

		// find difference in distance
		LD = LT - LQ;

		// find the length of gap and missmatch
		insertion_len = 0;
		deletion_len = 0;
		mismatch_len = 0;
		if (LD == 0)
		{
			mismatch_len = LT; // or LQ
		}
		else if (LD < 0)
		{
			insertion_len = -LD;
			mismatch_len = LT;
		}
		else // if (LD > 0)
		{
			deletion_len = LD;
			mismatch_len = LQ;
		}
		gap_len = insertion_len + deletion_len;
		return true;
	}
	
	void Print(int_t index)
	{
		Printf("\n%5d>\t%5d\t%5d\t%5d\t%5d", index, BQ, L, OFS, Si);
	}
};