CC=gcc
FLAGS=-Wall -O3 -lm -lfftw3 -DSHOUT

SRC=vecs.c
OBJ=$(SRC:%.c=%.o)

NAME=neuman

%.o: %.c %.h
	$(CC) $(FLAGS) -c $< -o $@

$(NAME): main.c $(OBJ)
	$(CC) $(FLAGS) $^ -o $@

ifneq (clean, $(MAKECMDGOALS))
-include deps.mk
endif

deps.mk: $(SRC)
	$(CC) -MM $^ > $@

clean:
	rm -f deps.mk
	rm -f *.o
	rm -f $(NAME)
