AM_CXXFLAGS = -fopenmp -O3 -std=c++17 -g -I${CONDA_PREFIX}/include/igraph -I$(top_srcdir)/external
komb2_LDADD = -lz -ligraph
komb2_LDFLAGS = -L${CONDA_PREFIX}/lib -Wl,-rpath, ${CONDA_PREFIX}/lib -fopenmp -pthread
bin_PROGRAMS = komb2
komb2_SOURCES = gfa.cpp graph.cpp komb2.cpp
bindir = $(top_srcdir)/bin