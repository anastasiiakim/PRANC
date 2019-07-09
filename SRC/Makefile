.SUFFIXES: .c .o
.SUFFIXES: .cpp .o

BINDIR = ../BIN
INCDIR = ../INCLUDE/
CC = g++

CFLAGS = -O2 -g -std=c++11 -I$(INCDIR)
CINC   =
LDFLAGS= -lm

LD=$(CC)

OBJS = main.o probs_calculation.o string_manipulations.o queue.o ranking_unr_gt.o min_ancient_coal.o write_ranked_tree.o maxlike.o
$(BINDIR)/pranc: $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $@
.c.o:
	$(CC) -c $(CFLAGS) $(CINC) $<
.cpp.o:
	$(CC) -c $(CFLAGS) $(CINC) $<

# The dependencies between modules
main.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)min_ancient_coal.h $(INCDIR)write_ranked_tree.h $(INCDIR)maxlike.h
probs_calculation.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h
ranking_unr_gt.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h
maxlike.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)maxlike.h
string_manipulations.cpp: $(INCDIR)string_manipulations.h
min_ancient_coal.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)min_ancient_coal.h
write_ranked_tree.cpp: $(INCDIR)node.h $(INCDIR)string_manipulations.h $(INCDIR)queue.h $(INCDIR)probs_calculation.h $(INCDIR)ranking_unr_gt.h $(INCDIR)write_ranked_tree.h
clean:
	rm -f *.o $(BINDIR)*.x
