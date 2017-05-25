#include "define.h"
#include "seqence.h"
#include "mem.h"
#include "aligner.h"
#include "FastaFile.h"

SEQUENCE temp_seq[2];
ticks et[ET_SIZE];

int Test_StaticInput()
{
	char t[129] = "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACGGTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTTTGACTGACTGACTGACTGACTGACTGACTGACTGACTG";
	char q[129] = "GTCAGGACCGACTGATGGACTGACTGACTTTGACTGACTGACTGACTGACTCACTTACAGACTGACTGACAAGTACTGACTGACTGACTAGCTGACTGACTGACTGACTGACTGACTGCCTGATCAGT";
	ALIGNER aligner(128, 1000);
	aligner.Set_Penalties(2, 2, 1);
	aligner.Set_Parameter(5, 5, 20, 100);
	aligner.Align(t, q);
	aligner.Generate_Alignment();
	//printf("\n%d\n\n", aligner.oa_score);

	// example to show marginal gaps
	strcpy(t, "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACGGTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTTTGACTGACTGACTGACTGACTGACTGACTGACTGACTG");
	strcpy(q, "GTCAGGACCGACTGATGGACTGACTGACTTTGACTGACTGACTGACTGACTCACTTACAGACTGACTGACAAGTACTGACTGACTGACTAGCTGACTGACTGACTGACTGACTGACTGCTGATCAGTG");
	aligner.Align(t, q);
	aligner.Generate_Alignment();
	//printf("\n%d\n\n", aligner.oa_score);

	aligner.Set_Parameter(5, 3, 20, 100);

	strcpy(t, "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACGGTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTTTGACTGACTGACTGACTGACTGACTGACTGACTGACTG");
	strcpy(q, "GTCAGGACCGACTGATGGACTGACTGACTTTGACTGACTGACTGACTGACTCACTTACAGACTGACTGACAAGTACTGACTGACTGACTAGCTGACTGACTGACTGACTGACTGACTGCTGATCAGTG");
	aligner.Align(t, q);
	aligner.Generate_Alignment();
	//printf("\n%d\n\n", aligner.oa_score);

	// example to show gaps in the middle
	aligner.Set_Parameter(5, 7, 20, 100);
	strcpy(t, "ATCGATCGATATCATTATTCGGGATCGATATTCGTATAGCTATGCTATGCTTATTCGATTGCTATGCTTAGCTTATGCTATTTCGATTTACGTATATCGATTACTGATCGTATATCGATGACTGATCG");
	strcpy(q, "ATCGATCGATATCATTATTCGGGATCGATATTCGTATAGCTATGGTATGCATTATTAGATTGCTATGCTTAGCTTATGCTATTTCGATTTACGTATATCGATTACTGATCGTATATCGATGACTGATC");
	aligner.Align(t, q);
	aligner.Generate_Alignment();
	//printf("\n%d\n\n", aligner.oa_score);

	aligner.Set_Parameter(5, 3, 20, 100);
	strcpy(t, "ATCGATCGATATCATTATTCGGGATCGATATTCGTATAGCTATGCTATGCTTATTCGATTGCTATGCTTAGCTTATGCTATTTCGATTTACGTATATCGATTACTGATCGTATATCGATGACTGATCG");
	strcpy(q, "ATCGATCGATATCATTATTCGGGATCGATATTCGTATAGCTATGGTATGCATTATTAGATTGCTATGCTTAGCTTATGCTATTTCGATTTACGTATATCGATTACTGATCGTATATCGATGACTGATC");
	aligner.Align(t, q);
	aligner.Generate_Alignment();
	//printf("\n%d\n\n", aligner.oa_score);

	// example to show two gaps
	aligner.Set_Parameter(5, 7, 20, 100);
	strcpy(t, "ATCGATCGATATCATTATTCGGGATCGATATTCGTATAGCTATGCTATGCTTATTCGATTGCTATGCTTAGCTTATGCTATTTCGATTTACGTATATCGATTACTGATCGTATATCGATGACTGATCG");
	strcpy(q, "ATCGATCGATATCATTATTCGGGATCGATATTCGTATAGCTATGCTATGCTTCATTCGTTGCTATGCTTAGCTTATGCTATTTCGATTTACGTATATCGATTACTGATCGTATATCGATGACTGATCG");
	aligner.Align(t, q);
	aligner.Generate_Alignment();
	//printf("\n%d\n\n", aligner.oa_score);

	aligner.Set_Parameter(5, 3, 20, 100);
	strcpy(t, "ATCGATCGATATCATTATTCGGGATCGATATTCGTATAGCTATGCTATGCTTATTCGATTGCTATGCTTAGCTTATGCTATTTCGATTTACGTATATCGATTACTGATCGTATATCGATGACTGATCG");
	strcpy(q, "ATCGATCGATATCATTATTCGGGATCGATATTCGTATAGCTATGCTATGCTTCATTCGTTGCTATGCTTAGCTTATGCTATTTCGATTTACGTATATCGATTACTGATCGTATATCGATGACTGATCG");
	aligner.Align(t, q);
	aligner.Generate_Alignment();
	//printf("\n%d\n\n", aligner.oa_score);

	return 0;
}

int Test_Load(int argc, char *argv[])
{
	Printf("\n Start: %s", __FUNCTION__);
	if (argc < 4)
		Exit_Args;

	FASTA_FILE ft(argv[2], "rb");
	FASTA_FILE fq(argv[3], "rb");

	int_t len = 0;
	if (sscanf(argv[4], "%d", &len) != 1)
		Exit_Args;

	uint32 total_seq = 0;
	while (true)
	{
		int tidx = ft.GetSeq();
		int qidx = fq.GetSeq();
		if (tidx == -1 || qidx == -1)
			break;
		total_seq++;
	}

	Printf("\n %d sequence pairs are loaded\n", total_seq);
	return 0;
}

int Test_Convert(int argc, char *argv[])
{
	Printf("\n Start: %s", __FUNCTION__);
	if (argc < 5)
		Exit_Args;

	FASTA_FILE ft(argv[2], "rb");
	FASTA_FILE fq(argv[3], "rb");

	int_t len = 0;
	if (sscanf(argv[4], "%d", &len) != 1)
		Exit_Args;

	ALIGNER aligner(len, 10);
	aligner.Set_Penalties(4, 6, 1);
	aligner.Set_Parameter(20, 2, 20, 80);

	uint32 total_seq = 0;
	while (true)
	{
		int tidx = ft.GetSeq();
		int qidx = fq.GetSeq();
		if (tidx == -1 || qidx == -1)
			break;
		total_seq++;
		aligner.Align_Convert(ft.seqs[tidx], fq.seqs[qidx]);
	}

	Printf("\n %d sequence pairs are converted\n", total_seq);
	return 0;
}

int Test_Convert_Loop(int argc, char *argv[])
{
	Printf("\n Start: %s", __FUNCTION__);
	if (argc < 5)
		Exit_Args;

	FASTA_FILE ft(argv[2], "rb");
	FASTA_FILE fq(argv[3], "rb");

	int_t len = 0;
	if (sscanf(argv[4], "%d", &len) != 1)
		Exit_Args;

	ALIGNER aligner(len, 10);
	aligner.Set_Penalties(4, 6, 1);
	aligner.Set_Parameter(20, 2, 20, 80);

	uint32 total_seq = 0;
	while (true)
	{
		int tidx = ft.GetSeq();
		int qidx = fq.GetSeq();
		if (tidx == -1 || qidx == -1)
			break;
		total_seq++;
		aligner.Align_Convert_Loop(ft.seqs[tidx], fq.seqs[qidx]);
	}

	Printf("\n %d sequence pairs are converted\n", total_seq);
	return 0;
}

int Test_ExtractMEM(int argc, char *argv[])
{
	Printf("\n Start: %s", __FUNCTION__);
	if (argc < 8)
		Exit_Args;

	FASTA_FILE ft(argv[2], "rb");
	FASTA_FILE fq(argv[3], "rb");

	int_t len = 0;
	if (sscanf(argv[4], "%d", &len) != 1)
		Exit_Args;

	int_t gl = 0;
	if (sscanf(argv[5], "%d", &gl) != 1)
		Exit_Args;

	int_t sl = 0;
	if (sscanf(argv[6], "%d", &sl) != 1)
		Exit_Args;

	int_t TM = 0;
	if (sscanf(argv[7], "%d", &TM) != 1)
		Exit_Args;

	ALIGNER aligner(len, (TM + 10));
	aligner.Set_Penalties(4, 6, 1);
	aligner.Set_Parameter(gl, sl, 20, TM);

	uint32 total_seq = 0;
	double avg_mem = 0;
	while (true)
	{
		int tidx = ft.GetSeq();
		int qidx = fq.GetSeq();
		if (tidx == -1 || qidx == -1)
			break;
		total_seq++;
		aligner.Align_ExtractMEM(ft.seqs[tidx], fq.seqs[qidx]);
		if (aligner.num_mems <= aligner.TM)
			avg_mem += aligner.num_mems;
	}
	if(total_seq > aligner.bypass_SSW_by_TM)
		avg_mem /= (total_seq - aligner.bypass_SSW_by_TM);
	else
		avg_mem = 0;

	Printf("\n %d sequence pairs are processed\n", total_seq);
	Printf("\n %d sequence pairs are exceed TM\n", aligner.bypass_SSW_by_TM);
	Printf("\n %6.2f MEM exist per sequence pair on average\n", avg_mem);

	Printf("\n%d\t%d\t%6.2f\tAWKREP\n", total_seq, aligner.bypass_SSW_by_TM, avg_mem);
	return 0;
}

int Test_SortMEM(int argc, char *argv[])
{
	Printf("\n Start: %s", __FUNCTION__);
	if (argc < 8)
		Exit_Args;

	FASTA_FILE ft(argv[2], "rb");
	FASTA_FILE fq(argv[3], "rb");

	int_t len = 0;
	if (sscanf(argv[4], "%d", &len) != 1)
		Exit_Args;

	int_t gl = 0;
	if (sscanf(argv[5], "%d", &gl) != 1)
		Exit_Args;

	int_t sl = 0;
	if (sscanf(argv[6], "%d", &sl) != 1)
		Exit_Args;

	int_t TM = 0;
	if (sscanf(argv[7], "%d", &TM) != 1)
		Exit_Args;

	ALIGNER aligner(len, (TM + 10));
	aligner.Set_Penalties(4, 6, 1);
	aligner.Set_Parameter(gl, sl, 20, TM);

	uint32 total_seq = 0;
	while (true)
	{
		int tidx = ft.GetSeq();
		int qidx = fq.GetSeq();
		if (tidx == -1 || qidx == -1)
			break;
		total_seq++;
		aligner.Align_ExtractMEM(ft.seqs[tidx], fq.seqs[qidx]);
		if (aligner.num_mems <= aligner.TM)
			aligner.Sort_MEM();
	}

	Printf("\n %d sequence pairs are processed\n", total_seq);
	
	return 0;
}

int Test_QuickSortMEM(int argc, char *argv[])
{
	Printf("\n Start: %s", __FUNCTION__);
	if (argc < 8)
		Exit_Args;

	FASTA_FILE ft(argv[2], "rb");
	FASTA_FILE fq(argv[3], "rb");

	int_t len = 0;
	if (sscanf(argv[4], "%d", &len) != 1)
		Exit_Args;

	int_t gl = 0;
	if (sscanf(argv[5], "%d", &gl) != 1)
		Exit_Args;

	int_t sl = 0;
	if (sscanf(argv[6], "%d", &sl) != 1)
		Exit_Args;

	int_t TM = 0;
	if (sscanf(argv[7], "%d", &TM) != 1)
		Exit_Args;

	ALIGNER aligner(len, (TM + 10));
	aligner.Set_Penalties(4, 6, 1);
	aligner.Set_Parameter(gl, sl, 20, TM);

	uint32 total_seq = 0;
	while (true)
	{
		int tidx = ft.GetSeq();
		int qidx = fq.GetSeq();
		if (tidx == -1 || qidx == -1)
			break;
		total_seq++;
		aligner.Align_ExtractMEM(ft.seqs[tidx], fq.seqs[qidx]);
		if (aligner.num_mems <= aligner.TM)
			aligner.Quick_Sort_MEM();
	}

	Printf("\n %d sequence pairs are processed\n", total_seq);

	return 0;
}

int Test_Performance(int argc, char *argv[])
{
	Printf("\n Start: %s", __FUNCTION__);
	if (argc < 12)
		Exit_Args;

	FASTA_FILE ft(argv[2], "rb");
	FASTA_FILE fq(argv[3], "rb");

	int_t len = 0;
	if (sscanf(argv[4], "%d", &len) != 1)
		Exit_Args;

	int_t gl = 0;
	if (sscanf(argv[5], "%d", &gl) != 1)
		Exit_Args;

	int_t sl = 0;
	if (sscanf(argv[6], "%d", &sl) != 1)
		Exit_Args;

	int_t TM = 0;
	if (sscanf(argv[7], "%d", &TM) != 1)
		Exit_Args;

	int_t TD = 0;
	if (sscanf(argv[8], "%d", &TD) != 1)
		Exit_Args;

	int_t Px = 0;
	if (sscanf(argv[9], "%d", &Px) != 1)
		Exit_Args;

	int_t Po = 0;
	if (sscanf(argv[10], "%d", &Po) != 1)
		Exit_Args;

	int_t Pe = 0;
	if (sscanf(argv[11], "%d", &Pe) != 1)
		Exit_Args;

	OUT_FILE fas(argv[12], "wb");

	ALIGNER aligner(len, (TM + 10));
	aligner.Set_Penalties(Px, Po, Pe);
	aligner.Set_Parameter(gl, sl, TD, TM);

	uint32 total_seq = 0;
	while (true)
	{
		int tidx = ft.GetSeq();
		int qidx = fq.GetSeq();
		if (tidx == -1 || qidx == -1)
			break;
		total_seq++;
		aligner.Align(ft.seqs[tidx], fq.seqs[qidx]);
		fas.Write(aligner.oa_score);
	}

	Printf("\n%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\tAWKREP\n", total_seq, aligner.bypass_SSW_by_TM, aligner.num_SJI_normal, aligner.num_SJI_Omega, aligner.num_SJI_Omega_TD, aligner.numNeededSeqCmp, aligner.numSeqCmpAfterPreCompute);

	return 0;
}

int Test_Accuracy(int argc, char *argv[])
{
	Printf("\n Start Test Accuracy");
	Printf("\n The input should be interleaved T and Q sequences, all of the same length");
	Printf("\n Argument to this command: input_file_name Sequence_length gl sl TM TD Px Po Pe output_prefix input_golden_alignment_score_file");
	if (argc < 13)
		Exit_Args;

	FILE *fin = fopen(argv[2], "r");
	if (!fin)
		Exit_File;

	int_t len = 0;
	if (sscanf(argv[3], "%d", &len) != 1)
		Exit_Args;

	int_t gl = 0;
	if (sscanf(argv[4], "%d", &gl) != 1)
		Exit_Args;

	int_t sl = 0;
	if (sscanf(argv[5], "%d", &sl) != 1)
		Exit_Args;

	int_t TM = 0;
	if (sscanf(argv[6], "%d", &TM) != 1)
		Exit_Args;

	int_t TD = 0;
	if (sscanf(argv[7], "%d", &TD) != 1)
		Exit_Args;

	int_t Px = 0;
	if (sscanf(argv[8], "%d", &Px) != 1)
		Exit_Args;

	int_t Po = 0;
	if (sscanf(argv[9], "%d", &Po) != 1)
		Exit_Args;

	int_t Pe = 0;
	if (sscanf(argv[10], "%d", &Pe) != 1)
		Exit_Args;

	char fn[1000];
	sprintf(fn, "%s.AS", argv[11]);
	FILE *fout = fopen(fn, "w");
	if (!fout)
		Exit_File;
	sprintf(fn, "%s.T.fa", argv[11]);
	FILE *fT = fopen(fn, "w");
	if (!fout)
		Exit_File;
	sprintf(fn, "%s.Q.fa", argv[11]);
	FILE *fQ = fopen(fn, "w");
	if (!fout)
		Exit_File;

	FILE *fas = fopen(argv[12], "r");
	if (!fas)
		Exit_File;


	char *t = new char[len + 100];
	char *q = new char[len + 100];
	ALIGNER aligner(len, (TM + 10));
	aligner.Set_Penalties(Px, Po, Pe);
	aligner.Set_Parameter(gl, sl, TD, TM);
	//aligner.mmlioa = 4;

	int_t	sequence_pair_cnt = 0;
	int_t   golden_alignment_score;
	int_t   num_inacurate_AS=0;
	//int_t	num_unreliable = 0;
	while (fscanf(fin, "%s %s", t, q) == 2)
	{
		if (fscanf(fas, "%d", &golden_alignment_score) != 1)
			Exit("Error reading alignment score file");
		sequence_pair_cnt++;
		aligner.Align(t, q);
		fprintf(fout, "%d\n", aligner.oa_score);
		if (aligner.num_mems <= TM && aligner.oa_score != golden_alignment_score)
			num_inacurate_AS++;

		//if (aligner.not_reliable)
		//	num_unreliable++;
		
		if (aligner.num_mems > TM)
		{
			fprintf(fT, ">%d\n%s\n", sequence_pair_cnt, t);
			fprintf(fQ, ">%d\n%s\n", sequence_pair_cnt, q);
		}
	}

	fclose(fin);
	fclose(fout);
	fclose(fT);
	fclose(fQ);
	fclose(fas);
	//delete[]t;
	//delete[]q;
	aligner.num_SJI_normal;
	aligner.num_SJI_Omega;
	aligner.num_SJI_Omega_TD;
	aligner.numNeededSeqCmp;
	aligner.numSeqCmpAfterPreCompute;

	Printf("\nS1\t%llu", aligner.num_SJI_normal);
	Printf("\nS2\t%llu", aligner.num_SJI_Omega);
	Printf("\nS3\t%llu", aligner.num_SJI_Omega_TD);
	Printf("\nS4\t%llu", aligner.numNeededSeqCmp);
	Printf("\nS5\t%llu", aligner.numSeqCmpAfterPreCompute);
	//Printf("\nUntrustable\t%u", num_unreliable);
	Printf("\n %d sequence pairs processed. Stats: %d\t%d\n", sequence_pair_cnt, aligner.bypass_SSW_by_TM, num_inacurate_AS);

	return 0;
}

int Test_Alignemnt(int argc, char *argv[])
{
	Printf("\n Start Alignment Test");
	Printf("\n The input should be interleaved T and Q sequences, all of the same length");
	Printf("\n Argument to this command: input_file_name Sequence_length gl sl TM TD Px Po Pe");
	if (argc < 11)
		Exit_Args;

	FILE *fin = fopen(argv[2], "r");
	if (!fin)
		Exit_File;

	int_t len = 0;
	if (sscanf(argv[3], "%d", &len) != 1)
		Exit_Args;

	int_t gl = 0;
	if (sscanf(argv[4], "%d", &gl) != 1)
		Exit_Args;

	int_t sl = 0;
	if (sscanf(argv[5], "%d", &sl) != 1)
		Exit_Args;

	int_t TM = 0;
	if (sscanf(argv[6], "%d", &TM) != 1)
		Exit_Args;

	int_t TD = 0;
	if (sscanf(argv[7], "%d", &TD) != 1)
		Exit_Args;

	int_t Px = 0;
	if (sscanf(argv[8], "%d", &Px) != 1)
		Exit_Args;

	int_t Po = 0;
	if (sscanf(argv[9], "%d", &Po) != 1)
		Exit_Args;

	int_t Pe = 0;
	if (sscanf(argv[10], "%d", &Pe) != 1)
		Exit_Args;

	char *t = new char[len + 100];
	char *q = new char[len + 100];
	ALIGNER aligner(len, (TM + 10));
	aligner.Set_Penalties(Px, Po, Pe);
	aligner.Set_Parameter(gl, sl, TD, TM);

	int_t	sequence_pair_cnt = 0;
	//int_t   golden_alignment_score;
	int_t   num_inacurate_AS=0;
	Printf("\n\n");
	while (fscanf(fin, "%s %s", t, q) == 2)
	{
		sequence_pair_cnt++;
		aligner.Align(t, q);
		if (aligner.num_mems <= TM)
		{
			num_inacurate_AS++;
			aligner.Generate_Alignment();
			Printf("\n\nAlignemnt Score is: %d\n\n", aligner.oa_score);
		}
	}

	fclose(fin);
	//delete[]t;
	//delete[]q;
	Printf("\n %d sequence pairs processed");
	return 0;
}

int Test_Monitor(int argc, char *argv[])
{
	Printf("\n Start Test Extract MEM");
	Printf("\n The input should be interleaved T and Q sequences, all of the same length");
	Printf("\n Argument to this command: input_file_name Sequence_length gl sl TM");
	if (argc < 7)
		Exit_Args;

	FILE *fin = fopen(argv[2], "r");
	if (!fin)
		Exit_File;

	int_t len = 0;
	if (sscanf(argv[3], "%d", &len) != 1)
		Exit_Args;

	int_t gl = 0;
	if (sscanf(argv[4], "%d", &gl) != 1)
		Exit_Args;

	int_t sl = 0;
	if (sscanf(argv[5], "%d", &sl) != 1)
		Exit_Args;

	int_t TM = 0;
	if (sscanf(argv[6], "%d", &TM) != 1)
		Exit_Args;

	char *t = new char[len + 100];
	char *q = new char[len + 100];
	ALIGNER aligner(len, (TM + 10));
	aligner.Set_Penalties(4, 6, 1);
	aligner.Set_Parameter(gl, sl, 20, TM);

	int_t sequence_pair_cnt = 0;
	double avg_mem=0;
	while (fscanf(fin, "%s %s", t, q) == 2)
	{
		sequence_pair_cnt++;
		if (sequence_pair_cnt % 10000 == 1)
			Printf("\n%d", sequence_pair_cnt);
		aligner.Align_ExtractMEM(t, q);
		if (aligner.num_mems <= aligner.TM)
			avg_mem += aligner.num_mems;
	}
	avg_mem /= (sequence_pair_cnt - aligner.bypass_SSW_by_TM);

	fclose(fin);
	//delete[]t;
	//delete[]q;
	Printf("\n %d sequence pairs processed. %d pairs excceds TM. %6.3f average MEM per pair\n", sequence_pair_cnt, aligner.bypass_SSW_by_TM, avg_mem);
	return 0;
}

int Test_ET(int argc, char *argv[])
{
	Printf("\n Start: %s", __FUNCTION__);
	if (argc < 12)
		Exit_Args;

	FASTA_FILE ft(argv[2], "rb");
	FASTA_FILE fq(argv[3], "rb");

	int_t len = 0;
	if (sscanf(argv[4], "%d", &len) != 1)
		Exit_Args;

	int_t gl = 0;
	if (sscanf(argv[5], "%d", &gl) != 1)
		Exit_Args;

	int_t sl = 0;
	if (sscanf(argv[6], "%d", &sl) != 1)
		Exit_Args;

	int_t TM = 0;
	if (sscanf(argv[7], "%d", &TM) != 1)
		Exit_Args;

	int_t TD = 0;
	if (sscanf(argv[8], "%d", &TD) != 1)
		Exit_Args;

	int_t TS = 0;
	if (sscanf(argv[9], "%d", &TS) != 1)
		Exit_Args;

	int_t Px = 0;
	if (sscanf(argv[10], "%d", &Px) != 1)
		Exit_Args;

	int_t Po = 0;
	if (sscanf(argv[11], "%d", &Po) != 1)
		Exit_Args;

	int_t Pe = 0;
	if (sscanf(argv[12], "%d", &Pe) != 1)
		Exit_Args;

	OUT_FILE fas(argv[13], "wb");

	ALIGNER aligner(len, (TM + 10));
	aligner.Set_Penalties(Px, Po, Pe);
	aligner.Set_Parameter(gl, sl, TD, TM, TS);

	uint32 total_seq = 0;
	double avg_mem=0;
	while (true)
	{
		int tidx = ft.GetSeq();
		int qidx = fq.GetSeq();
		if (tidx == -1 || qidx == -1)
			break;
		total_seq++;
		add_ticks(ET_IO);
		aligner.Align(ft.seqs[tidx], fq.seqs[qidx]);
		fas.Write(aligner.oa_score); 
		if (aligner.processed_by == MEM_ALIGN)
			avg_mem += aligner.num_mems;
	}
	
	if (total_seq > aligner.bypass_SSW_by_TM)
		avg_mem /= (total_seq - aligner.bypass_SSW_total);
	else
		avg_mem = 0;

	aligner.bypass_SSW_by_noMEM = aligner.bypass_SSW_total - aligner.bypass_SSW_by_TM - aligner.bypass_SSW_by_TS;
	
	double SJI_reduced_by_Omega         = aligner.num_SJI_normal ? ((double)(aligner.num_SJI_normal - aligner.num_SJI_Omega   ) * 100) / aligner.num_SJI_normal : 0;
	double SJI_reduced_by_TD_over_Omega = aligner.num_SJI_Omega ? ((double)(aligner.num_SJI_Omega - aligner.num_SJI_Omega_TD) * 100) / aligner.num_SJI_Omega : 0;
	double SeqCmp_optimisation = aligner.numNeededSeqCmp ? ((double)(aligner.numNeededSeqCmp - aligner.numSeqCmpAfterPreCompute) * 100) / aligner.numNeededSeqCmp : 0;

	Printf("\nAWKREP");
	Printf("\t%llu", total_seq);
	Printf("\t%10.6f", avg_mem);
	Printf("\t%llu", aligner.bypass_SSW_by_TM);
	Printf("\t%llu", aligner.bypass_SSW_by_TS);
	Printf("\t%llu", aligner.bypass_SSW_by_noMEM);
	Printf("\t%llu", aligner.bypass_SSW_total);
	Printf("\t%10.6f", (double)(aligner.bypass_SSW_by_TM * 100) / total_seq);
	Printf("\t%10.6f", (double)(aligner.bypass_SSW_by_TS * 100) / total_seq);
	Printf("\t%10.6f", (double)(aligner.bypass_SSW_by_noMEM * 100) / total_seq);
	Printf("\t%llu", aligner.num_SJI_normal);
	Printf("\t%llu", aligner.num_SJI_Omega);
	Printf("\t%llu", aligner.num_SJI_Omega_TD);
	Printf("\t%10.6f", SJI_reduced_by_Omega);
	Printf("\t%10.6f", SJI_reduced_by_TD_over_Omega);
	Printf("\t%llu", aligner.numNeededSeqCmp);
	Printf("\t%llu", aligner.numSeqCmpAfterPreCompute);
	Printf("\t%10.6f", SeqCmp_optimisation);
	Printf("\t%llu", et[ET_IO]);
	Printf("\t%llu", et[ET_CONVERT]);
	Printf("\t%llu", et[ET_EXTRACT]);
	Printf("\t%llu", et[ET_SORT]);
	Printf("\t%llu", et[ET_ALN]);
	Printf("\t%llu", et[ET_BACKTRACK]);
	Printf("\t%llu", et[ET_SW]);
	Printf("\n");

	return 0;
}

int Test_ET_TEXT_REPORT(int argc, char *argv[])
{
	Printf("\n Start: %s", __FUNCTION__);
	if (argc < 12)
		Exit_Args;

	FASTA_FILE ft(argv[2], "rb");
	FASTA_FILE fq(argv[3], "rb");

	int_t len = 0;
	if (sscanf(argv[4], "%d", &len) != 1)
		Exit_Args;

	int_t gl = 0;
	if (sscanf(argv[5], "%d", &gl) != 1)
		Exit_Args;

	int_t sl = 0;
	if (sscanf(argv[6], "%d", &sl) != 1)
		Exit_Args;

	int_t TM = 0;
	if (sscanf(argv[7], "%d", &TM) != 1)
		Exit_Args;

	int_t TD = 0;
	if (sscanf(argv[8], "%d", &TD) != 1)
		Exit_Args;

	int_t TS = 0;
	if (sscanf(argv[9], "%d", &TS) != 1)
		Exit_Args;

	int_t Px = 0;
	if (sscanf(argv[10], "%d", &Px) != 1)
		Exit_Args;

	int_t Po = 0;
	if (sscanf(argv[11], "%d", &Po) != 1)
		Exit_Args;

	int_t Pe = 0;
	if (sscanf(argv[12], "%d", &Pe) != 1)
		Exit_Args;

	OUT_FILE fas(argv[13], "wb");

	ALIGNER aligner(len, (TM + 10));
	aligner.Set_Penalties(Px, Po, Pe);
	aligner.Set_Parameter(gl, sl, TD, TM, TS);

	uint32 total_seq = 0;
	double avg_mem = 0;
	while (true)
	{
		int tidx = ft.GetSeq();
		int qidx = fq.GetSeq();
		if (tidx == -1 || qidx == -1)
			break;
		total_seq++;
		add_ticks(ET_IO);
		aligner.Align(ft.seqs[tidx], fq.seqs[qidx]);
		fas.Write(aligner.oa_score);
		if (aligner.processed_by == MEM_ALIGN)
			avg_mem += aligner.num_mems;
	}

	if (total_seq > aligner.bypass_SSW_by_TM)
		avg_mem /= (total_seq - aligner.bypass_SSW_total);
	else
		avg_mem = 0;

	aligner.bypass_SSW_by_noMEM = aligner.bypass_SSW_total - aligner.bypass_SSW_by_TM - aligner.bypass_SSW_by_TS;

	double SJI_reduced_by_Omega = aligner.num_SJI_normal ? ((double)(aligner.num_SJI_normal - aligner.num_SJI_Omega) * 100) / aligner.num_SJI_normal : 0;
	double SJI_reduced_by_TD_over_Omega = aligner.num_SJI_Omega ? ((double)(aligner.num_SJI_Omega - aligner.num_SJI_Omega_TD) * 100) / aligner.num_SJI_Omega : 0;
	double SeqCmp_optimisation = aligner.numNeededSeqCmp ? ((double)(aligner.numNeededSeqCmp - aligner.numSeqCmpAfterPreCompute) * 100) / aligner.numNeededSeqCmp : 0;

	Printf("\n");
	Printf("\n%50s:\t%llu", "Total number of input sequence pair", total_seq);
	Printf("\n%50s:\t%10.6f", "Average Number of Extracted MEM", avg_mem);
	Printf("\n%50s:\t%llu", "Bypassed by TM", aligner.bypass_SSW_by_TM);
	Printf("\n%50s:\t%llu", "Bypassed by TS", aligner.bypass_SSW_by_TS);
	Printf("\n%50s:\t%llu", "Bypassed by no MEM", aligner.bypass_SSW_by_noMEM);
	Printf("\n%50s:\t%llu", "Bypassed in total", aligner.bypass_SSW_total);
	Printf("\n%50s:\t%10.6f", "Percentage of Bypassed by TM", (double)(aligner.bypass_SSW_by_TM * 100) / total_seq);
	Printf("\n%50s:\t%10.6f", "Percentage of Bypassed by TS", (double)(aligner.bypass_SSW_by_TS * 100) / total_seq);
	Printf("\n%50s:\t%10.6f", "Percentage of Bypassed by no MEM", (double)(aligner.bypass_SSW_by_noMEM * 100) / total_seq);
	Printf("\n%50s:\t%llu", "Number of expected alignment extention", aligner.num_SJI_normal);
	Printf("\n%50s:\t%llu", "Number of expected alignment extention after Omega", aligner.num_SJI_Omega);
	Printf("\n%50s:\t%llu", "Number of executed alignment extention after Omega and TD", aligner.num_SJI_Omega_TD);
	Printf("\n%50s:\t%10.6f", "Percentage of Alignemnt extension avoided by Omega", SJI_reduced_by_Omega);
	Printf("\n%50s:\t%10.6f", "Percentage of Alignemnt extension avoided by TD", SJI_reduced_by_TD_over_Omega);
	Printf("\n%50s:\t%llu", "Number of expected String Compare operation", aligner.numNeededSeqCmp);
	Printf("\n%50s:\t%llu", "Number of Executed String Compare optimisation", aligner.numSeqCmpAfterPreCompute);
	Printf("\n%50s:\t%10.6f", "percentage of avoided String Compare", SeqCmp_optimisation);
	Printf("\n%50s:\t%llu", "Number of CPU cycle spent on IO", et[ET_IO]);
	Printf("\n%50s:\t%llu", "Number of CPU cycle spent on Convert", et[ET_CONVERT]);
	Printf("\n%50s:\t%llu", "Number of CPU cycle spent on Extract", et[ET_EXTRACT]);
	Printf("\n%50s:\t%llu", "Number of CPU cycle spent on Sort", et[ET_SORT]);
	Printf("\n%50s:\t%llu", "Number of CPU cycle spent on Alignment", et[ET_ALN]);
	Printf("\n%50s:\t%llu", "Number of CPU cycle spent on Backtracking", et[ET_BACKTRACK]);
	Printf("\n%50s:\t%llu", "Number of CPU cycle spent on SSW", et[ET_SW]);
	Printf("\n");

	return 0;
}

int Test_ET_CS_A(int argc, char *argv[])
{
	Printf("\n Start: %s", __FUNCTION__);
	if (argc < 12)
		Exit_Args;

	FASTA_FILE ft(argv[2], "rb");
	FASTA_FILE fq(argv[3], "rb");

	int_t len = 0;
	if (sscanf(argv[4], "%d", &len) != 1)
		Exit_Args;

	int_t gl = 0;
	if (sscanf(argv[5], "%d", &gl) != 1)
		Exit_Args;

	int_t sl = 0;
	if (sscanf(argv[6], "%d", &sl) != 1)
		Exit_Args;

	int_t TM = 0;
	if (sscanf(argv[7], "%d", &TM) != 1)
		Exit_Args;

	int_t TD = 0;
	if (sscanf(argv[8], "%d", &TD) != 1)
		Exit_Args;

	int_t TS = 0;
	if (sscanf(argv[9], "%d", &TS) != 1)
		Exit_Args;

	int_t Px = 0;
	if (sscanf(argv[10], "%d", &Px) != 1)
		Exit_Args;

	int_t Po = 0;
	if (sscanf(argv[11], "%d", &Po) != 1)
		Exit_Args;

	int_t Pe = 0;
	if (sscanf(argv[12], "%d", &Pe) != 1)
		Exit_Args;

	OUT_FILE fas(argv[13], "wb");

	ALIGNER aligner(len, (TM + 10));
	aligner.Set_Penalties(Px, Po, Pe);
	aligner.Set_Parameter(gl, sl, TD, TM, TS);

	uint32 total_seq = 0;
	double avg_mem = 0;
	while (true)
	{
		int tidx = ft.GetSeq();
		int qidx = fq.GetSeq();
		if (tidx == -1 || qidx == -1)
			break;
		total_seq++;
		add_ticks(ET_IO);
		aligner.Align_CS_A(ft.seqs[tidx], fq.seqs[qidx]);
		//fas.Write(aligner.oa_score);
		//if (aligner.processed_by == MEM_ALIGN)
			avg_mem += aligner.num_mems;
	}

	if (total_seq > aligner.bypass_SSW_by_TM)
		avg_mem /= (total_seq - aligner.bypass_SSW_total);
	else
		avg_mem = 0;

	aligner.bypass_SSW_by_noMEM = aligner.bypass_SSW_total - aligner.bypass_SSW_by_TM - aligner.bypass_SSW_by_TS;

	double SJI_reduced_by_Omega = aligner.num_SJI_normal ? ((double)(aligner.num_SJI_normal - aligner.num_SJI_Omega) * 100) / aligner.num_SJI_normal : 0;
	double SJI_reduced_by_TD_over_Omega = aligner.num_SJI_Omega ? ((double)(aligner.num_SJI_Omega - aligner.num_SJI_Omega_TD) * 100) / aligner.num_SJI_Omega : 0;
	double SeqCmp_optimisation = aligner.numNeededSeqCmp ? ((double)(aligner.numNeededSeqCmp - aligner.numSeqCmpAfterPreCompute) * 100) / aligner.numNeededSeqCmp : 0;

	Printf("\nAWKREP");
	Printf("\t%llu", total_seq);
	Printf("\t%10.6f", avg_mem);
	Printf("\t%llu", aligner.bypass_SSW_by_TM);
	Printf("\t%llu", aligner.bypass_SSW_by_TS);
	Printf("\t%llu", aligner.bypass_SSW_by_noMEM);
	Printf("\t%llu", aligner.bypass_SSW_total);
	Printf("\t%10.6f", (double)(aligner.bypass_SSW_by_TM * 100) / total_seq);
	Printf("\t%10.6f", (double)(aligner.bypass_SSW_by_TS * 100) / total_seq);
	Printf("\t%10.6f", (double)(aligner.bypass_SSW_by_noMEM * 100) / total_seq);
	Printf("\t%llu", aligner.num_SJI_normal);
	Printf("\t%llu", aligner.num_SJI_Omega);
	Printf("\t%llu", aligner.num_SJI_Omega_TD);
	Printf("\t%10.6f", SJI_reduced_by_Omega);
	Printf("\t%10.6f", SJI_reduced_by_TD_over_Omega);
	Printf("\t%llu", aligner.numNeededSeqCmp);
	Printf("\t%llu", aligner.numSeqCmpAfterPreCompute);
	Printf("\t%10.6f", SeqCmp_optimisation);
	Printf("\t%llu", et[ET_IO]);
	Printf("\t%llu", et[ET_CONVERT]);
	Printf("\t%llu", et[ET_EXTRACT]);
	Printf("\t%llu", et[ET_SORT]);
	Printf("\t%llu", et[ET_ALN]);
	Printf("\t%llu", et[ET_BACKTRACK]);
	Printf("\t%llu", et[ET_SW]);
	Printf("\n");

	return 0;
}

int Test_ET_CS_B(int argc, char *argv[])
{
	Printf("\n Start: %s", __FUNCTION__);
	if (argc < 12)
		Exit_Args;

	FASTA_FILE ft(argv[2], "rb");
	FASTA_FILE fq(argv[3], "rb");

	int_t len = 0;
	if (sscanf(argv[4], "%d", &len) != 1)
		Exit_Args;

	int_t gl = 0;
	if (sscanf(argv[5], "%d", &gl) != 1)
		Exit_Args;

	int_t sl = 0;
	if (sscanf(argv[6], "%d", &sl) != 1)
		Exit_Args;

	int_t TM = 0;
	if (sscanf(argv[7], "%d", &TM) != 1)
		Exit_Args;

	int_t TD = 0;
	if (sscanf(argv[8], "%d", &TD) != 1)
		Exit_Args;

	int_t TS = 0;
	if (sscanf(argv[9], "%d", &TS) != 1)
		Exit_Args;

	int_t Px = 0;
	if (sscanf(argv[10], "%d", &Px) != 1)
		Exit_Args;

	int_t Po = 0;
	if (sscanf(argv[11], "%d", &Po) != 1)
		Exit_Args;

	int_t Pe = 0;
	if (sscanf(argv[12], "%d", &Pe) != 1)
		Exit_Args;

	OUT_FILE fas(argv[13], "wb");

	ALIGNER aligner(len, (TM + 10));
	aligner.Set_Penalties(Px, Po, Pe);
	aligner.Set_Parameter(gl, sl, TD, TM, TS);

	uint32 total_seq = 0;
	double avg_mem = 0;
	while (true)
	{
		int tidx = ft.GetSeq();
		int qidx = fq.GetSeq();
		if (tidx == -1 || qidx == -1)
			break;
		total_seq++;
		add_ticks(ET_IO);
		aligner.Align_CS_B(ft.seqs[tidx], fq.seqs[qidx]);
		//fas.Write(aligner.oa_score);
		//if (aligner.processed_by == MEM_ALIGN)
			avg_mem += aligner.num_mems;
	}

	if (total_seq > aligner.bypass_SSW_by_TM)
		avg_mem /= (total_seq - aligner.bypass_SSW_total);
	else
		avg_mem = 0;

	aligner.bypass_SSW_by_noMEM = aligner.bypass_SSW_total - aligner.bypass_SSW_by_TM - aligner.bypass_SSW_by_TS;

	double SJI_reduced_by_Omega = aligner.num_SJI_normal ? ((double)(aligner.num_SJI_normal - aligner.num_SJI_Omega) * 100) / aligner.num_SJI_normal : 0;
	double SJI_reduced_by_TD_over_Omega = aligner.num_SJI_Omega ? ((double)(aligner.num_SJI_Omega - aligner.num_SJI_Omega_TD) * 100) / aligner.num_SJI_Omega : 0;
	double SeqCmp_optimisation = aligner.numNeededSeqCmp ? ((double)(aligner.numNeededSeqCmp - aligner.numSeqCmpAfterPreCompute) * 100) / aligner.numNeededSeqCmp : 0;

	Printf("\nAWKREP");
	Printf("\t%llu", total_seq);
	Printf("\t%10.6f", avg_mem);
	Printf("\t%llu", aligner.bypass_SSW_by_TM);
	Printf("\t%llu", aligner.bypass_SSW_by_TS);
	Printf("\t%llu", aligner.bypass_SSW_by_noMEM);
	Printf("\t%llu", aligner.bypass_SSW_total);
	Printf("\t%10.6f", (double)(aligner.bypass_SSW_by_TM * 100) / total_seq);
	Printf("\t%10.6f", (double)(aligner.bypass_SSW_by_TS * 100) / total_seq);
	Printf("\t%10.6f", (double)(aligner.bypass_SSW_by_noMEM * 100) / total_seq);
	Printf("\t%llu", aligner.num_SJI_normal);
	Printf("\t%llu", aligner.num_SJI_Omega);
	Printf("\t%llu", aligner.num_SJI_Omega_TD);
	Printf("\t%10.6f", SJI_reduced_by_Omega);
	Printf("\t%10.6f", SJI_reduced_by_TD_over_Omega);
	Printf("\t%llu", aligner.numNeededSeqCmp);
	Printf("\t%llu", aligner.numSeqCmpAfterPreCompute);
	Printf("\t%10.6f", SeqCmp_optimisation);
	Printf("\t%llu", et[ET_IO]);
	Printf("\t%llu", et[ET_CONVERT]);
	Printf("\t%llu", et[ET_EXTRACT]);
	Printf("\t%llu", et[ET_SORT]);
	Printf("\t%llu", et[ET_ALN]);
	Printf("\t%llu", et[ET_BACKTRACK]);
	Printf("\t%llu", et[ET_SW]);
	Printf("\n");

	return 0;
}

int Test_Normal(int argc, char *argv[])
{
	Printf("\n Start: %s", __FUNCTION__);
	if (argc < 12)
		Exit_Args;

	FASTA_FILE ft(argv[2], "rb");
	FASTA_FILE fq(argv[3], "rb");

	int_t len = 0;
	if (sscanf(argv[4], "%d", &len) != 1)
		Exit_Args;

	int_t gl = 0;
	if (sscanf(argv[5], "%d", &gl) != 1)
		Exit_Args;

	int_t sl = 0;
	if (sscanf(argv[6], "%d", &sl) != 1)
		Exit_Args;

	int_t TM = 0;
	if (sscanf(argv[7], "%d", &TM) != 1)
		Exit_Args;

	int_t TD = 0;
	if (sscanf(argv[8], "%d", &TD) != 1)
		Exit_Args;

	int_t TS = 0;
	if (sscanf(argv[9], "%d", &TS) != 1)
		Exit_Args;

	int_t Px = 0;
	if (sscanf(argv[10], "%d", &Px) != 1)
		Exit_Args;

	int_t Po = 0;
	if (sscanf(argv[11], "%d", &Po) != 1)
		Exit_Args;

	int_t Pe = 0;
	if (sscanf(argv[12], "%d", &Pe) != 1)
		Exit_Args;

	OUT_FILE fas(argv[13], "wb");

	ALIGNER aligner(len, (TM + 10));
	aligner.Set_Penalties(Px, Po, Pe);
	aligner.Set_Parameter(gl, sl, TD, TM, TS);

	uint32 total_seq = 0;
	double avg_mem = 0;
	while (true)
	{
		int tidx = ft.GetSeq();
		int qidx = fq.GetSeq();
		if (tidx == -1 || qidx == -1)
			break;
		total_seq++;
		add_ticks(ET_IO);
		aligner.Align(ft.seqs[tidx], fq.seqs[qidx]);
		fas.Write(aligner.oa_score);
		if (aligner.processed_by == MEM_ALIGN)
			avg_mem += aligner.num_mems;
		aligner.Generate_Alignment();
	}

	if (total_seq > aligner.bypass_SSW_by_TM)
		avg_mem /= (total_seq - aligner.bypass_SSW_total);
	else
		avg_mem = 0;

	aligner.bypass_SSW_by_noMEM = aligner.bypass_SSW_total - aligner.bypass_SSW_by_TM - aligner.bypass_SSW_by_TS;

	double SJI_reduced_by_Omega = aligner.num_SJI_normal ? ((double)(aligner.num_SJI_normal - aligner.num_SJI_Omega) * 100) / aligner.num_SJI_normal : 0;
	double SJI_reduced_by_TD_over_Omega = aligner.num_SJI_Omega ? ((double)(aligner.num_SJI_Omega - aligner.num_SJI_Omega_TD) * 100) / aligner.num_SJI_Omega : 0;
	double SeqCmp_optimisation = aligner.numNeededSeqCmp ? ((double)(aligner.numNeededSeqCmp - aligner.numSeqCmpAfterPreCompute) * 100) / aligner.numNeededSeqCmp : 0;

	Printf("\nAWKREP");
	Printf("\t%llu", total_seq);
	Printf("\t%10.6f", avg_mem);
	Printf("\t%llu", aligner.bypass_SSW_by_TM);
	Printf("\t%llu", aligner.bypass_SSW_by_TS);
	Printf("\t%llu", aligner.bypass_SSW_by_noMEM);
	Printf("\t%llu", aligner.bypass_SSW_total);
	Printf("\t%10.6f", (double)(aligner.bypass_SSW_by_TM * 100) / total_seq);
	Printf("\t%10.6f", (double)(aligner.bypass_SSW_by_TS * 100) / total_seq);
	Printf("\t%10.6f", (double)(aligner.bypass_SSW_by_noMEM * 100) / total_seq);
	Printf("\t%llu", aligner.num_SJI_normal);
	Printf("\t%llu", aligner.num_SJI_Omega);
	Printf("\t%llu", aligner.num_SJI_Omega_TD);
	Printf("\t%10.6f", SJI_reduced_by_Omega);
	Printf("\t%10.6f", SJI_reduced_by_TD_over_Omega);
	Printf("\t%llu", aligner.numNeededSeqCmp);
	Printf("\t%llu", aligner.numSeqCmpAfterPreCompute);
	Printf("\t%10.6f", SeqCmp_optimisation);
	Printf("\t%llu", et[ET_IO]);
	Printf("\t%llu", et[ET_CONVERT]);
	Printf("\t%llu", et[ET_EXTRACT]);
	Printf("\t%llu", et[ET_SORT]);
	Printf("\t%llu", et[ET_ALN]);
	Printf("\t%llu", et[ET_BACKTRACK]);
	Printf("\t%llu", et[ET_SW]);
	Printf("\n");

	return 0;
}

int main(int argc, char *argv[])
{
	Printf("\nUsage: %s CMD Target_File Query_File Seq_Len gl sl TM TD TS Px Po Pe OutputFile(Alignment_Scores)", argv[0]);
	Printf("\n>>CMD=K process fixed sequences with fixed parametes. No need for additional argumnts.(used in alignment generation example of the paper) \n");
	Printf("\n>>CMD=X Only store alignment scores (used for performance and accuracy test in the paper) \n");
	Printf("\n>>CMD=L Print alignments in color along with cigar string and alignment score  (alignment generation is not optimised and should not be used for performance test) \n");
	Printf("\n>>CMD=W As same as X with human redable report \n");
	Printf("\n>>Other CMD are used for test purpose. You may want to read C source code for them. \n");
	Printf("\n===========================================================\n");

	init_ticks();
	
	if (argc < 2)
		return 0;

	switch (argv[1][0])
	{
	case 'A':
		Test_Load(argc, argv);
		break;
	case 'B':
		Test_Convert(argc, argv);
		break;	
	case 'C':
		Test_Convert_Loop(argc, argv);
		break;
	case 'D':
		Test_ExtractMEM(argc, argv);
		break;
	case 'E':
		Test_SortMEM(argc, argv);
		break;
	case 'F':
		Test_QuickSortMEM(argc, argv);
		break;
	case 'G':
		Test_Performance(argc, argv);
		break;
	case 'H':
		Test_Accuracy(argc, argv);
		break;
	case 'I':
		Test_Alignemnt(argc, argv);
		break;
	case 'J':
		Test_Monitor(argc, argv);
		break;
	case 'K':
		Test_StaticInput();
		break;
	case 'L':
		Test_Normal(argc, argv);
		break;
	case 'W':
		Test_ET_TEXT_REPORT(argc, argv);
		break;
	case 'X':
		Test_ET(argc, argv);
		break;
	case 'Y':
		Test_ET_CS_A(argc, argv);
		break;
	case 'Z':
		Test_ET_CS_B(argc, argv);
		break;
	default:
		break;
	}
	return 0;
}

#ifdef USE_LZCNT
	#ifdef __linux
		inline int_t LZCNT(word_t X)
		{
			return(int_t)__builtin_clzl(X);
		}
	#else
		inline int_t LZCNT(word_t X)
		{
			return(int_t)__lzcnt64(X);
		}
	#endif
#else
	inline int_t LZCNT(word_t X)
	{
		int_t cnt = 0;
		while (!(X & 0X8000000000000000llu))
		{
			if (cnt == 63)
				return -1;
			X <<= 1;
			cnt++;
		}
		return cnt;
	}
#endif

#ifdef USE_POPCNT
	#ifdef __linux
		inline int_t POPCNT(word_t X)
		{
			return ((int_t)__builtin_popcountl(X));
		}
	#else
		inline int_t POPCNT(word_t X)
		{
			return ((int_t)__popcnt64(X));
		}
	#endif
#else
	inline int_t POPCNT(word_t X) { return POPCNT_MUL(X); }

	inline int_t POPCNT_MUL(word_t X)
	{
		X |= X >> 1;
		X &= ALL_5;
		X = (X & ALL_3) + (X >> 2 & ALL_3);
		X = (X & ALL_0F) + ((X >> 4) &  ALL_0F);
		X = ((X * ALL_01) >> 56) & 0xFF;
		return X;
	}

	inline int_t POPCNT_ADD(word_t X)
	{
		X = (X & 0x3333333333333333llu) + (X >> 2 & 0x3333333333333333llu);
		X = ((X + (X >> 4)) & 0x0F0F0F0F0F0F0F0Fllu);
		X = ((X + (X >> 8)) & 0x00FF00FF00FF00FFllu);
		X = ((X + (X >> 16)) & 0x0000FFFF0000FFFFllu);
		X = ((X + (X >> 32)) & 0x00000000FFFFFFFFllu);
		return X;
	}
#endif