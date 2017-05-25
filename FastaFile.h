#pragma once
#include "define.h"

#define bufsize (16*1024*1024)

inline void reverse(char str[], int length)
{
	int start = 0;
	int end = length - 1;
	while (start < end)
	{
		char temp = str[start];
		str[start] = str[end];
		str[end] = temp;
		start++;
		end--;
	}
}

inline int my_itoa(int num, char* str)
{
	#define base (10)
	int i = 0;
	bool isNegative = false;

	if (num == 0)
	{
		str[i++] = '0';
		str[i] = '\n';
		return 2;
	}

	if (num < 0 && base == 10)
	{
		isNegative = true;
		num = -num;
	}

	while (num != 0)
	{
		int rem = num % base;
		str[i++] = rem + '0';
		num = num / base;
	}

	// If number is negative, append '-'
	if (isNegative)
		str[i++] = '-';

	str[i] = '\n'; // Append string terminator

				   // Reverse the string
	reverse(str, i);

	return (i+1);
}

class FASTA_FILE
{
public:
	bool file_end;
	uint64 file_pos;
	FILE* file_pointer;
	char* buf;
	char** seqs;
	char** ids;
	uint32 seq_cnt;
	uint32 idx;


	FASTA_FILE(char *file_name, char *mode)
	{
		file_pointer = fopen(file_name, mode);
		if (!file_pointer)
		{
			printf("\n **** Error cannot open file: %s ", file_name);
			exit(0);
		}
		file_end = false;
		file_pos = 0;
		buf = new char[bufsize];
		seqs = new char *[bufsize];
		ids = new char *[bufsize];
		seq_cnt = 0;
		idx = 0;
	}
	~FASTA_FILE()
	{
		delete[]buf;
		delete[]seqs;
		delete[]ids;
		fclose(file_pointer);
	};

	bool Fasta_Read()
	{
		if (file_end)
			return false;
		rewind(file_pointer);
		fseek(file_pointer, file_pos, SEEK_SET);
		uint32 read_cnt = fread(buf, 1, bufsize, file_pointer);
		if (read_cnt < bufsize)
			file_end = true;
		uint32 line_cnt = 0;
		seq_cnt = 0;
		idx = 0;
		ids[seq_cnt] = buf;
		uint64 pre_seq_pos = file_pos;
		for (uint32 i = 0; i < read_cnt; i++)
		{
			if (buf[i] == '\n')
			{
				buf[i] = 0;
				line_cnt++;
				if ((line_cnt & 0x1))
				{
					seqs[seq_cnt] = &buf[i+1];
				}
				else
				{
					seq_cnt++;
					ids[seq_cnt] = &buf[i+1];
					pre_seq_pos = file_pos;
					//cur_seq_pos = file_pos;
				}
			}
			file_pos++;
		}
		file_pos = pre_seq_pos+1;
		return true;
	};
	
	int GetSeq()
	{
		int ret = -1;
		if (seq_cnt == 0)
		{
			if (!Fasta_Read())
				return ret;
		}
		if (idx < seq_cnt)
		{
			ret = idx;
			idx++;
		}
		else
		{
			if (Fasta_Read())
			{
				ret = 0;
				idx = 1;
			}
		}
		return ret;
	}
	//uint Write_Text_Block(char *buf, uint len);
};

class OUT_FILE
{
public:
	FILE* file_pointer;
	char* buf;
	uint32 idx;

	OUT_FILE(char *file_name, char *mode)
	{
		file_pointer = fopen(file_name, mode);
		if (!file_pointer)
		{
			printf("\n **** Error cannot open file: %s ", file_name);
			exit(0);
		}
		buf = new char[bufsize];
		idx = 0;
	}
	~OUT_FILE()
	{
		fwrite(buf, 1, idx, file_pointer);
		delete[]buf;
		fclose(file_pointer);
	};

	bool Write(int32 score)
	{
		idx += my_itoa(score, &buf[idx]);

		if (idx > bufsize - 100)
		{
			fwrite(buf, 1, idx, file_pointer);
			idx = 0;
		}

		return true;
	};

	//uint Write_Text_Block(char *buf, uint len);
};