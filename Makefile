BIN=fp

$(BIN) : fp.c
	$(CC) -o $@ $^ -lm

