# Makefile for duplex_analyzer
#
# Requires htslib.  If htslib is not in the default system paths, set:
#   make HTSLIB_PREFIX=/path/to/htslib

CC        = gcc
CFLAGS    = -O2 -Wall -Wextra -std=c99
TARGET    = duplex_analyzer
SRCS      = duplex_analyzer.c kmer.c
OBJS      = duplex_analyzer.o kmer.o

# Optional: point to a local htslib install
# HTSLIB_PREFIX = /usr/local
ifdef HTSLIB_PREFIX
	CFLAGS  += -I$(HTSLIB_PREFIX)/include
	LDFLAGS += -L$(HTSLIB_PREFIX)/lib -Wl,-rpath,$(HTSLIB_PREFIX)/lib
endif

LDLIBS = -lhts -lz -lm -lpthread

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS) tiny.h kmer.h
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)

duplex_analyzer.o: duplex_analyzer.c tiny.h kmer.h
	$(CC) $(CFLAGS) -c -o $@ duplex_analyzer.c

kmer.o: kmer.c kmer.h
	$(CC) $(CFLAGS) -c -o $@ kmer.c

clean:
	rm -f $(TARGET) $(OBJS)
