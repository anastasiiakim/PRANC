.SUFFIXES: .c .o
.SUFFIXES: .cpp .o

BINDIR = ../BIN
INCDIR = ../INCLUDE/
CC = g++

CFLAGS = -O2 -g -std=c++11 -I$(INCDIR)
CINC   =
LDFLAGS= -lm

LD=$(CC)

LBFGSB  = lbfgsb.o linesearch.o subalgorithms.o print.o

LINPACK = linpack.o

BLAS 	= miniCBLAS.o
#POSTFLAGS += -D_USE_OPTIMIZED_BLAS -lblas

TIMER   = timer.o

SRC = $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(F2CFILES)

OBJS = main.o probs_calculation.o string_manipulations.o queue.o ranking_unr_gt.o min_ancient_coal.o write_ranked_tree.o maxlike.o maxlike_brent.o symbolic.o get_unranked_topology.o get_ranked_topology.o get_greedy_consensus.o rank_distance.o branch_distance.o $(SRC)
$(BINDIR)/pranc: $(OBJS) $(SRC)
	$(LD) $(OBJS) $(LDFLAGS) -o $@
.c.o:
	$(CC) -c $(CFLAGS) $(CINC) $<
.cpp.o:
	$(CC) -c $(CFLAGS) $(CINC) $<




# The dependencies between modules
main.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)min_ancient_coal.h $(INCDIR)write_ranked_tree.h $(INCDIR)maxlike.h $(INCDIR)symbolic.h $(INCDIR)get_unranked_topology.h $(INCDIR)get_ranked_topology.h $(INCDIR)get_greedy_consensus.h $(INCDIR)rank_distance.h $(INCDIR)branch_distance.h
string_manipulations.cpp: $(INCDIR)string_manipulations.h
probs_calculation.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h
ranking_unr_gt.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h
symbolic.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)symbolic.h
maxlike.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)maxlike.h 
maxlike_brent.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)maxlike.h $(INCDIR)maxlike_brent.h
min_ancient_coal.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)min_ancient_coal.h
write_ranked_tree.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)write_ranked_tree.h
get_unranked_topology.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)ranking_unr_gt.h $(INCDIR)probs_calculation.h $(INCDIR)get_unranked_topology.h
get_ranked_topology.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)probs_calculation.h $(INCDIR)get_ranked_topology.h
get_greedy_consensus.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)ranking_unr_gt.h $(INCDIR)probs_calculation.h $(INCDIR)get_unranked_topology.h $(INCDIR)write_ranked_tree.h $(INCDIR)get_greedy_consensus.h
rank_distance.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)probs_calculation.h $(INCDIR)rank_distance.h
branch_distance.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)probs_calculation.h $(INCDIR)branch_distance.h

clean:
	rm -f *.o $(BINDIR)*.x
