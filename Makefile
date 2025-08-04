CC=gcc
CFLAGS=-Wall -Wextra -std=c99
OBJECTS=main.c dna.c dna.h
NAME=encoder

program: $(OBJECTS)
	$(CC) -o $(NAME) $(OBJECTS)

clean:
	rm -f *.o $(NAME)
