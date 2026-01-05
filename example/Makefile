CC = gcc
CFLAGS = -O3 -Wall -Wextra -std=c99 -pedantic
LDFLAGS = -lm

MOWAO_PATH = ../mowao/
MOWAO = mowao
CEXMPL = example3.c

all: c-example

$(MOWAO):
	$(CC) $(CFLAGS) -c $(MOWAO_PATH)$@.c -o $@.o

c-example: $(MOWAO)
	$(CC) $(CFLAGS) $^.o $(CEXMPL) -o example $(LDFLAGS)
	
clean: 
	rm -vf example *.o *.so
