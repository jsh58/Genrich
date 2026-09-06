Genrich: Genrich.c Genrich.h
	$(CC) -g -Wall -std=gnu99 $(LDFLAGS) $(CFLAGS) -O2 -o Genrich Genrich.c -lz -lm
