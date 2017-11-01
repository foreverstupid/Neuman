CC=gcc
FLAGS=-Wall -O3 -lm -lfftw3 -DSHOUT

SRC=vecs.c
OBJ=$(SRC:%.c=%.o)

NAME=neuman
VALUES_PATH=./values
IMG_PATH=./images

%.o: %.c %.h
	$(CC) $(FLAGS) -c $< -o $@

$(NAME): main.c $(OBJ)
	$(CC) $(FLAGS) $^ -o $@

ifneq (clean, $(MAKECMDGOALS))
-include deps.mk
endif

deps.mk: $(SRC)
	$(CC) -MM $^ > $@

calc:
	if [ ! -e $(NAME) ]; then make; fi
	if [ ! -d $(VALUES_PATH) ]; then mkdir $(VALUES_PATH); fi
	./calculate.sh ./$(NAME) $(VALUES_PATH)

draw:
	if [ ! -d $(VALUES_PATH) ]; then make calc; fi
	if [ -z `ls $(VALUES_PATH) | head -n 1` ]; then make calc; fi
	if [ ! -d $(IMG_PATH) ]; then mkdir $(IMG_PATH); fi
	./draw.sh $(VALUES_PATH) $(IMG_PATH)

clean:
	rm -f deps.mk
	rm -f *.o
	rm -f $(NAME)
	rm -rf $(IMG_PATH)
	rm -rf $(VALUES_PATH)

