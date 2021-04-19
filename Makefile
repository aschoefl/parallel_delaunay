.PHONY = all clean

CC = mpicxx
TARGET = delauney
OBJECTS = main.o polygon.o


# all: $(OBJECTS)
# 	$(CC) $(INCLUDES) $(CFLAGS) -o $(OBJECT_OUT) $(OBJECTS)

all: main.o polygon.o point.o
	$(CC) -o $(TARGET) main.o polygon.o point.o
	rm -rf *.o

main.o: main.cpp
	$(CC) -c main.cpp -I ./

polygon.o: polygon.cpp
	$(CC) -c polygon.cpp -I ./

point.o: point.cpp
	$(CC) -c point.cpp -I ./

clean:
	rm -rf $(TARGET)

