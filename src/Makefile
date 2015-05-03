CC      = gcc
CFLAGS  = -Wall -Wextra -O3 -march=native -std=c99
CFLAGS += -pedantic -g
#FFLAGS += $(shell pkg-config --cflags plplotd-f95)
#LDFLAGS = $(shell pkg-config --libs plplotd-f95)
LIBS    = -lm

COMPILE = $(CC) $(CFLAGS)
LINK = $(CC) $(LDFLAGS)

TARGET = reacnet.exe       # Name of final executable to produce
OBJS = reacnet.o nucleus.o reaction.o  # List of object dependencies

$(TARGET): $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o:%.c
	$(COMPILE) -c $<

clean: 
	rm $(OBJS)
