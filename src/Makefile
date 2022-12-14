CC = g++
CFLAGS = -O3 -std=c++11

TARGET = tp3
INCLUDES = -I ../src/  

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