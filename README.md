# Masters Dissertation Project #
Queen Mary University of London
MSc Mathematical Finance

Purpose of this project to optimize numerical solutions of parabolic PDEs by testing high performance computing techniques and comparing compilers/os/32bit/64bit.
The idea of this project is to study how to take advantage of this parallelism and explore how much faster we can make these calculations.

# TODO #

delta s + random 10^-20 prevent undesired optimization 
log spot bs bak compare errors

error için precision arttır
Generalize ADI implement multi assset BS

farkli tridiagonal solverlara bak thomas algo cyclic reduction double sweep
Intel numerical library for tridiagonal systems

32 bit vs 64  bit
Intel compiler vs VSCompiler vs gcc
Optimization switches

OpenMP
Multithreading (HPC slidesa bak)
AVX/Intrinsics  (internal to the CPU, probably only through a Visual Studio compilation switch, see SIMD Registers section in attached slides) 

Timinge gore heatmap hazirla

cash-flushes: When you ask for a piece of data from memory, as getting stuff from memory is slow but can be vectorised, a whole load of consecutive values are copied from RAM into the cache. This means that if you are doing a matrix operation that uses the first column of values, then the second etc. It is much better to store these values in memory as column 1, then column 2 etc. If you save it as Row1, Row 2 etc, then every time you access a new element the CPU will have to flush the cache: delete all its values and again load a set of consecutive values from memory. The differences can be very big. So you should have a look at how you implement the ADI method and decide whether the way the calculations are done, it is much better to have the arrays stored by a list of columns or of rows. And it could be (but I have not checked this, so this is just an optimistic dream) that ArrayOld and ArrayNew need different storage. If this were to be so then this would be a quite interesting discovery.

# FURTHER WORK #
Cloud applications
More type of options BSPde bdegistir callputflag duffy c++ 677 
More models
rannacher trick (168 foreign exchange pricing rannacher ve stencil)
knuth premature optimization


Numerical Solution of Partial Differential Equations: Finite Difference Methods, Gordon D. Smith

Numerical Partial Differential Equations: Finite Difference Methods, Thomas, J.W.

Modern X86 Assembly Language Programming, Daniel Kusswurm

Inner Loops: A Sourcebook for Fast 32-bit Software Development Paperback, Rick Booth  


