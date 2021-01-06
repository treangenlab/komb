AM_CXXFLAGS = -fopenmp -O3 -std=c++17 -g -I/home/software/miniconda3/include/igraph -I/home/software/miniconda3/include/bifrost -I$(top_srcdir)/external
komb_LDADD = -lz -ligraph -lbifrost
komb_LDFLAGS = -L/home/software/miniconda3/lib -Wl,-rpath=/home/software/miniconda3/lib -fopenmp -pthread
bin_PROGRAMS = komb
komb_SOURCES = gfa.cpp  graph.cpp komb.cpp
bindir = $(top_srcdir)/bin

