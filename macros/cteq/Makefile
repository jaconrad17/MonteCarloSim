.PHONY : clean

SHELL = /bin/sh
CC    = gcc
CFLAGS = -fPIC -pipe -pthread -O2 
LDFLAGS= -shared

SOURCES = $(shell echo *.c)
HEADERS = $(shell echo *.h)
OBJECTS = $(SOURCES:.c=.o)

TARGET  = libcteq.so

all: $(TARGET)

clean: 
	rm -f $(OBJECTS) $(TARGET)

$(TARGET): $(OBJECTS) 
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)
