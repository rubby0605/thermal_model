CC = g++
LD = g++
SRCS = $(wildcard *.cpp)
OBJS = $(wildcard *.cpp)

INCLUDE = -I/usr/local/include/eigen3/
LIB = -L.
CFLAGS = -O0 -g
TARGET= main
.PHONY:all clean

all: $(TARGET)
$(TARGET): $(OBJS)
	$(LD) -o $@ $^ $(INCLUDE) -I. $(LIB) $(CFLAGS)

%.o:%.cpp
	$(CC) -c $^ $(INCLUDE) -I. $(CFLAGS)

clean:
	rm -rf $(TARGET)

