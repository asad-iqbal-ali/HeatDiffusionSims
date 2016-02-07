OUTPUT=heat_test
MAINFILE=driver.c
HEADERS=class-func.h heatsolvers.h
DATAFILES=
OBJFILES=hs.o nrutil.o
CC=gcc

CFLAGS =-g -lm

all: $(OUTPUT)

$(OUTPUT): driver.c class-func.c hs.o nrutil.o
	$(CC) -o $(OUTPUT) hs.o $(CFLAGS) driver.c class-func.c nrutil.o -lrt

hs.o: heatsolvers.c heatsolvers.h
	$(CC) -o hs.o -c heatsolvers.c $(CFLAGS)

nrutil.o: nrutil.c nrutil.h
	$(CC) -o nrutil.o -c nrutil.c
clean:
	-@rm -f $(OUTPUT) $(DATAFILES) $(OBJFILES)
