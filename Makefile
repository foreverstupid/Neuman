CFLAGS = -O3 -DNDEBUG -DSHOUT -DN_COUNT=10001 -DITERS=500

all: pnm.x
	true

pnm.x: p_neuman.c
	gcc p_neuman.c -o pnm.x -lfftw3 -lm $(CFLAGS)

clean:
	rm -f pnm.x
