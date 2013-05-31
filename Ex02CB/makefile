CC = g++
CFLAGS = -O3 -Wall -Wshadow -Winline -ansi -pedantic 
LIBS = 
INCLUDES = -I./Colsamm
TARGET = waveguide

SRC = $(wildcard *.cpp)
OBJS = $(patsubst %.c, %.o, $(SRC))

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LIBS)

.PHONY : clean depend

clean: 
	@/bin/rm -f $(OBJS) 
	@/bin/rm -f $(TARGET)

depend: 
	@makedepend -- $(CFLAGS) -- $(SRC)