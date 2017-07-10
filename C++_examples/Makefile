CC = g++
CFLAGS = -Weverything
LIBS = -std=c++11 -stdlib=libc++

PROGS = lambda
OBJS = lambda.o

.PHONY: ALL clean

ALL: $(PROGS)

$(PROGS): %: %.o
	$(CC) $(LIBS) $(CFLAGS) $< -o $@

%.o: %.cpp
	$(CC) -c $(LIBS) $(CFLAGS) $< -o $@

clean:
	rm $(PROGS) $(PROGS).o
