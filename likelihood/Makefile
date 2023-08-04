.POSIX:
.SUFFIXES:
CFLAGS = -O2 -g
S = main.c
likelihood: $S
	$(CC) $S -o $@ $(CFLAGS) -lm $(LDFLAGS)
clean:; rm -f likelihood
