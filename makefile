# Makefile for fluffy-eureka
# Based on the exceptional boilerplate from http://stackoverflow.com/a/1484873/3280431

TARGET = fluffy-eureka
LIBS = -lm
CC = gcc
CFLAGS = -g -std=gnu11 -Wall -Wextra -Wpedantic -O3 # Let's make GCC very strict.
#CFLAGS = -std=c99 

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)