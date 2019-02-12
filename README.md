KmerEstimate: A Streaming Algorithm for Estimating k-mer Counts with Optimal Space Usage - Parallel Implementation
------------------------------------------------------------------------------------------------------------------

Counting distinct number of k-mers (substrings of length k in a DNA/RNA sequence), and computing the frequency distribution of k-mers (k-mer abundance histogram), in a genome data is central component in many bioinformatics applications. Some example applications are sequence assembly, read error correction, genome size prediction and estimation of its characteristics, changes in copy number of highly repetitive sequences, digital normalization, and parameter tuning of k-mer analysis based tools. People have tried several approaches for computing counts of all distinct k-mers present in sequence data like sorting, suffix-array, filter and sketch data structure, and parallel disk-based partitioning. These exact count algorithms did not scale well for large datasets because of time and space complexity. Instead approximate count algorithms give better results. Streaming algorithms have been proved to be efficient for approximating the k-mer abundance histogram and related counts because of their memory and space usage. The tools which use streaming algorithms approach for k-mer counting problem are KmerGenie, KmerStream, ntCard, Kmerlight. KmerEstimate is yet another streaming algorithm that approximates the number of k-mers with a given frequency in a genomic data set.


Hardware Requirements
---------------------

OS - Windows 10

RAM - 8GB

Software Requirements
---------------------

C++11

Cygwin

Sublime Text Editor

Installing
----------

zlib is designed to be a free, general-purpose, legally unencumbered -- that is, not covered by any patents -- lossless data-compression library for use on virtually any computer hardware and operating system. The link for downloading zlib is http://gnuwin32.sourceforge.net/packages/zlib.htm



Compile and Run Serial Implementation
------------------------------


Compile:


		g++ kmerCountEstimate_serial.cpp -lpthread -lz -std=c++11
		
		
Run:


		./a.exe -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -o <out.txt>
  
 Compile and Run Parallel Implementation
------------------------------ 

Compile:


		g++ kmerCountEstimate_parallel.cpp FastxParser.cpp -lpthread -lz -std=c++11
		
		
Run:


		./a.exe -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -o <out.txt> -np <NoOfProducers> -nt <NoOfConsumers> 
  
