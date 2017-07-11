Video Description on https://youtu.be/_SnPC793QMM

Description

MEM-Align is a fast, near optimal nucleotide sequence alignment algorithm using Maximal Exact Matches (MEMs)

Developed by: Arash Bayat (a.bayat@unsw.edu.au)

MEM-Align is under publication

Title: "Fast and accurate nucleotide sequence alignment algorithm suitable for DNA read mapping"

Authors: Arash Bayat, Bruno Gaeta, Aleksandar Ignjatovic and Sri Parameswaran.

University of New South Wales - School of Computer Science and Engineering

A technical report on MEM-Align is available on: ftp://ftp.cse.unsw.edu.au/pub/doc/papers/UNSW/201701.pdf

However, we recomend you to wait for the final publication which include latest changes.

website: https://sites.google.com/site/memalignv1/

Code Ocean: https://codeocean.com/2017/05/26/mem-align-colon-a-fast-near-optimal-nucleotide-sequence-alignment-algorithm/ ***(One-click simple run is available)***


*** In order to use CodeOcean please run the algorithm using the provided interface in Interface tab.

Parameters are as follow:

Program: (choosing what you want to run)
 1. MEM-Align (Test Mode): *** Produce colourful alignment for sample sequences (those used in the paper)
 2. MEM-Align (Generate Colourful Alignment): *** Produce colourful alignment for the given dataset and parameters (output file would be really large if large dataset is chosen we recommend to only execute this with "TinyQ" dataset)
 3. MEM-Align (Performance Mode): Produce alignment score for the given dataset and parameters.
 4. SSW (SIMD Smith-Waterman): Produce alignment score for the given dataset and scoring value using Smith-Waterman (SSW:https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
 5. MEM-Align and SSW (Compare results): Produce alignment score using both MEM-Align and Smith-Waterman and compare produced alignment scores and report the comparison results at the end.
 
*** Note that printing colourful text is not supported in Code Ocean Console. You can see colours if you run the program in Linux shell with a black background. In order to do so download the code from link in which you can use Makefile to compile the code.

Dataset: (refer to the paper)
 1. DSL: 1 million sequences of length 125 with low edit rate
 2. DSH: 1 million sequences of length 125 with low edit rate
 3. DLL: 1 million sequences of length 500 with low edit rate
 4. DLH: 1 million sequences of length 500 with low edit rate
 5. DRQ: 1 million sequences of length 250 with natural edit rate
 6. TinyQ: the first 1000 sequences of DRQ

gl: Band width (refer to the paper).

sl: Short MEM length threshold (refer to the paper).

TM: Maximum MEM threshold (refer to the paper).

TD: Maximum distance threshold (refer to the paper).

TS: Minimum score threshold (refer to the paper).
