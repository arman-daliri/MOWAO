CC = gcc
CFLAGS = -O3 -Wall -Wextra -std=c99 -pedantic
LDFLAGS = -lm

MOWAO = mowao

all: c-example py-example

$(MOWAO):
	$(CC) $(CFLAGS) -c $@.c -o $@.o

c-example: $(MOWAO).o
	$(CC) $(CFLAGS) $^ example.c -o example $(LDFLAGS)
	
py-example: fpy cpy
cpy: $(MOWAO)
	$(CC) $(CFLAGS) $(MOWAO).o -o $(MOWAO).so $(LDFLAGS)
fpy:
	$(eval CFLAGS += -fPIC -shared)

clean: 
	rm -vf example *.o *.so


