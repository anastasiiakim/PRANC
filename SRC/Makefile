.SUFFIXES: .c .o
.SUFFIXES: .cpp .o

BINDIR = ../BIN
INCDIR = ../INCLUDE/
CC = g++

CFLAGS = -O2 -g -std=c++11 -I$(INCDIR)
CINC   =
LDFLAGS= -lm

LD=$(CC)

OBJS = main.o probs_calculation.o string_manipulations.o queue.o ranking_unr_gt.o min_ancient_coal.o write_ranked_tree.o maxlike.o symbolic.o get_unranked_topology.o get_ranked_topology.o
$(BINDIR)/pranc: $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $@
.c.o:
	$(CC) -c $(CFLAGS) $(CINC) $<
.cpp.o:
	$(CC) -c $(CFLAGS) $(CINC) $<

# The dependencies between modules
main.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)min_ancient_coal.h $(INCDIR)write_ranked_tree.h $(INCDIR)maxlike.h $(INCDIR)symbolic.h $(INCDIR)get_unranked_topology.h $(INCDIR)get_ranked_topology.h
string_manipulations.cpp: $(INCDIR)string_manipulations.h
probs_calculation.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h
ranking_unr_gt.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h
symbolic.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)symbolic.h
maxlike.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)maxlike.h
min_ancient_coal.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)min_ancient_coal.h
write_ranked_tree.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)write_ranked_tree.h
get_unranked_topology.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)ranking_unr_gt.h $(INCDIR)probs_calculation.h $(INCDIR)get_unranked_topology.h
get_ranked_topology.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)probs_calculation.h $(INCDIR)get_ranked_topology.h

clean:
	rm -f *.o $(BINDIR)*.x
