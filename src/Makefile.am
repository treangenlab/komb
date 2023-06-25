AM_CXXFLAGS = -fopenmp -O3 -std=c++17 -g -I${CONDA_PREFIX}/include/igraph -I${CONDA_PREFIX}/include/bifrost -I$(top_srcdir)/external
komb_LDADD = -lz -ligraph -lbifrost
komb_LDFLAGS = -L${CONDA_PREFIX}/lib -Wl,-rpath,${CONDA_PREFIX}/lib -fopenmp -pthread
komb2_LDADD = -lz -ligraph -lbifrost
komb2_LDFLAGS = -L${CONDA_PREFIX}/lib -Wl,-rpath,${CONDA_PREFIX}/lib -fopenmp -pthread
bin_PROGRAMS = komb2
komb_SOURCES = gfa.cpp  graph.cpp komb.cpp
komb2_SOURCES = gfa.cpp  graph.cpp komb2.cpp
bindir = $(top_srcdir)
