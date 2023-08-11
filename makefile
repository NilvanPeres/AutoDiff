CC := gcc
CFLAGS := -O4 -I/home/junior/Documents/TCC/spg/ADOL-C/ADOL-C
LDFLAGS := -L/home/junior/Documents/TCC/spg/ADOL-C/ADOL-C/.libs -ladolc -lm

SRC := spg.c toyprob.c spgma.c
OBJ := $(SRC:.c=.o)
EXECUTABLE := spgma

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXECUTABLE) $(OBJ)
