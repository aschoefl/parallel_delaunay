.PHONY = all clean

CC = mpicxx
TARGET = delauney
OBJECTS = main.o polygon.o


# all: $(OBJECTS)
# 	$(CC) $(INCLUDES) $(CFLAGS) -o $(OBJECT_OUT) $(OBJECTS)

all: main.o p.o
	$(CC) -o $(TARGET) main.o polygon.o
	rm -rf *.o

main.o: main.cpp
	$(CC) -c main.cpp -I ./

p.o: polygon.cpp
	$(CC) -c polygon.cpp -I ./


clean:
	rm -rf $(TARGET)

