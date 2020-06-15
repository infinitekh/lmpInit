SHELL        = /bin/bash
CC           = gcc
FLAGS        =           # -std=gnu99
CFLAGS       = #-fPIC  -g # -pedantic -Wall -Wextra -march=native -ggdb3
LDFLAGS      = -lm -llmpio -lmymath  #-shared


TARGET       = lmpInit.out
SOURCES      = main.c
HEADERS      = 
OBJECTS      = $(SOURCES:.c=.o)

PREFIX     = $(DESTDIR)/usr/local
BINDIR     = $(PREFIX)/bin
LIBDIR     = $(PREFIX)/lib
INCLUDEDIR = $(PREFIX)/include

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS)  $(OBJECTS) -o $@ $(LDFLAGS)

debug: $(OBJECTS)
	$(CC) $(CFLAGS) -ggdb  $(SOURCES) -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJECTS) $(TARGET)

#.PHONY: install
#install: $(TARGET)
#	mkdir -p $(LIBDIR)
#	cp $< $(LIBDIR)/
#	cp $(HEADERS) $(INCLUDEDIR)/
#
#.PHONY: uninstall
#uninstall:
#	rm -f $(LIBDIR)/$(TARGET)
#	rm -f $(INCLUDEDIR)/$(HEADERS) 
#
	
