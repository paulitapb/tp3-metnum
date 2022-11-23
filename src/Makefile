CC = g++
CFLAGS = -pedantic -O3 -lm -Wall -Wextra -g
CFLAGS2 := -O3 -Wall -std=c++11

TARGET = tp3
INCLUDES = -I ../src/  -I /usr/include/python3.8

all: $(TARGET)

valgrind: 
	valgrind --leak-check=full  --undef-value-errors=yes ./tp2

$(TARGET): main_cpp.o eigen_cpp.o utils_cpp.o
	$(CC) $(CFLAGS) $^  -o $@ $(LDLIBS)

main_cpp.o: main.cpp eigen.h types.h utils.h
	$(CC) $(CFLAGS) -c $(INCLUDES) $< -o $@

eigen_cpp.o: eigen.cpp eigen.h types.h
	$(CC) $(CFLAGS) -c $(INCLUDES) $< -o $@

utils_cpp.o: utils.cpp utils.h types.h
	$(CC) $(CFLAGS) -c $(INCLUDES) $< -o $@

.PHONY: clean

clean:
	rm -f *.o
	rm -f tp3