.PHONY = all clean

CC = mpicxx
TARGET = delauney
OBJECTS = main.o polygon.o


# all: $(OBJECTS)
# 	$(CC) $(INCLUDES) $(CFLAGS) -o $(OBJECT_OUT) $(OBJECTS)

all: main.o point.o bucket.o
	$(CC) -o $(TARGET) main.o point.o bucket.o
	rm -rf *.o

main.o: main.cpp
	$(CC) -c main.cpp -I ./

point.o: point.cpp
	$(CC) -c point.cpp -I ./

bucket.o: bucket.cpp
	$(CC) -c bucket.cpp -I ./

clean:
	rm -rf $(TARGET)

