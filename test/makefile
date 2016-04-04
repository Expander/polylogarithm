CPPFLAGS := -Wall -Wextra -pedantic
CXXFLAGS := -I.. $(shell gsl-config --cflags)
LIBS     := $(shell gsl-config --libs)
DILOGLIB := ../dilog.o

all: test.x

clean:
	-rm -f test.x test.o

.PHONY: all clean

test.x: test.o $(DILOGLIB)
	$(CXX) -o $@ $^ $(LIBS)
