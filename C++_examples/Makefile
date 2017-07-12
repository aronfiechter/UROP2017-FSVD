CC = g++
CFLAGS = -Weverything -Wno-c++98-compat
LIBS = -std=c++11 -stdlib=libc++

PROGS = lambda reference
OBJS = lambda.o reference.o

.PHONY: ALL clean

ALL: $(PROGS)

$(PROGS): %: %.o
	$(CC) $(LIBS) $(CFLAGS) $< -o $@

%.o: %.cpp
	$(CC) -c $(LIBS) $(CFLAGS) $< -o $@

clean:
	rm $(PROGS) $(OBJS)
