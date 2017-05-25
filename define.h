#pragma once

#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#define ET_IO			1
#define ET_CONVERT		2
#define ET_EXTRACT		3
#define ET_SORT			4
#define ET_ALN			5
#define ET_SW			6
#define ET_BACKTRACK	7
#define	ET_SIZE			8


typedef unsigned long long ticks;
extern ticks et[ET_SIZE];

#ifdef TICKPROF
#include <sys/time.h>
#include <time.h>

static __inline__ void init_ticks(void)
{
	for (int i = 0; i < ET_SIZE; i++)
		et[i] = 0;

	unsigned a, d;
	
	//asm("cpuid");
	asm volatile("rdtsc" : "=a" (a), "=d" (d));
		
	et[0] = (((ticks)a) | (((ticks)d) << 32));
}

static __inline__ void add_ticks(int idx)
{
	unsigned a, d;

	//asm("cpuid");
	asm volatile("rdtsc" : "=a" (a), "=d" (d));
	
	ticks t = (((ticks)a) | (((ticks)d) << 32));
	et[idx] += t - et[0];
	et[0] = t;
}
#else

void init_ticks(void)
{
}

void add_ticks(int idx)
{
}
#endif

#define POWER2(X) ((unsigned long long int)1<<X)
#define IS_POWER_OF_2(X) ((X!=0)&&(!(X&(X-1))))

#define Printf printf
#define Exit(STR) {Printf("\n***Error: %s\n\n", STR); exit(1);}
#define Exit_Args {Printf("\n***Cannot parse arguments\n\n"); exit(1);}
#define Exit_File {Printf("\n***Cannot Open File\n\n"); exit(1);}

#define BIT_PER_BASE_LOG		1
#define BIT_PER_WORD_LOG		6	// 64-bit
#define BASE_PER_WORD_LOG		(BIT_PER_WORD_LOG-BIT_PER_BASE_LOG)
#define BIT_PER_BASE			POWER2(BIT_PER_BASE_LOG)
#define BIT_PER_WORD			POWER2(BIT_PER_WORD_LOG)
#define BASE_PER_WORD			POWER2(BASE_PER_WORD_LOG)
#define BIT_PER_WORD_MASK		(BIT_PER_WORD-1)
#define BYTE_PER_WORD			(BIT_PER_WORD>>3)

#define ALL_F  (0xFFFFFFFFFFFFFFFFllu)
#define ALL_5  (0x5555555555555555llu)
#define ALL_C  (0xCCCCCCCCCCCCCCCCllu)
#define ALL_A  (0xAAAAAAAAAAAAAAAAllu)
#define ALL_3  (0x3333333333333333llu)
#define ALL_0F (0x0F0F0F0F0F0F0F0Fllu)
#define ALL_F0 (0xF0F0F0F0F0F0F0F0llu)
#define ALL_01 (0x0101010101010101llu)
#define ALL_06 (0x0606060606060606llu)

typedef unsigned long long int uint64;
typedef unsigned int uint32;
typedef int int32;

// default integer type
typedef signed int int_t;

// machine word type
typedef unsigned long long int word_t;


#define USE_LZCNT
#define USE_POPCNT

#ifdef USE_LZCNT
#ifdef __linux

#else
	#include <intrin.h>  
	using namespace std;
	#pragma intrinsic(_BitScanReverse)  
	#pragma intrinsic(_BitScanReverse64)  
	#pragma intrinsic(_BitScanForward)  
	#pragma intrinsic(_BitScanForward64)  
#endif
#endif

#ifdef USE_POPCNT
#ifdef __linux
	#include <nmmintrin.h>
#else
	#include <intrin.h>  
	using namespace std;
	#pragma intrinsic(__popcnt64)  
	#pragma intrinsic(__popcnt)  
	#endif
#endif

inline int_t LZCNT(word_t X);
inline int_t POPCNT(word_t X);


#include "iostream"
#include "color.h"

#define Print_M(T,R) {if(R&1){std::cout<<termcolor::grey<<termcolor::on_cyan<<T<<termcolor::reset;}else{std::cout<<termcolor::grey<<termcolor::on_yellow<<T<<termcolor::reset;}}
#define Print_X(T) {std::cout<<termcolor::red<<T<<termcolor::reset;}
#define Print_RM(T) {std::cout<<termcolor::green<<T<<termcolor::reset;}
#define Print_C(T) {std::cout<<termcolor::red<<termcolor::on_white<<T<<termcolor::reset;}
#define Print_G(T) {std::cout<<T<<termcolor::reset;}

class INT128
{
public:
	uint64 l;
	uint64 h;
	INT128()
	{
		l = h = 0;
	}
	void inc()
	{
		//l++;
		//if (l == 0)
		//	h++;
	}
	void sh()
	{
		l >>= 7;
		l |= h << (64 - 7);
	}

};