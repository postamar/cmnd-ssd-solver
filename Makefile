# (C) Marius Posta, 2012
# Check LICENSE.txt for the legal blah-blah.

# variables

CFLAGS = -Wall -std=c99 -O3 
GRBROOT  = /Library/gurobi501/mac64
GRBFLAGS = -I$(GRBROOT)/include/ $(CFLAGS)
COMPILER = clang


# build rules

ssd: main.o cmnd.o shortestpath.o colgen.o
	$(COMPILER) $(GRBFLAGS) main.o cmnd.o shortestpath.o colgen.o -o ssd -L$(GRBROOT)/lib/ -lgurobi50 -lm -lpthread

main.o: main.c
	$(COMPILER) -c $(GRBFLAGS) main.c -o main.o 

cmnd.o: cmnd.c
	$(COMPILER) -c $(CFLAGS) cmnd.c -o cmnd.o 

shortestpath.o: shortestpath.c
	$(COMPILER) -c $(CFLAGS) shortestpath.c -o shortestpath.o 

colgen.o: colgen.c
	$(COMPILER) -c $(GRBFLAGS) colgen.c -o colgen.o


# delete temporary files

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log


