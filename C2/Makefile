CC ?= gcc # c-compiler
CFLAGS = # compiler flags such as -O3 -W
LIBS = -lgsl -lm # libraries such as -lm -lgsl -lgslcblas
SOURCES = linalg.c gsl_example.c
OBJECTS = add.o main.o # object files to produce
PROGRAM = program # name of executable to be built

LOCAL = ${HOME}/.local
CFLAGS += \
	-Wl,-rpath,$(LOCAL)/lib \
	-I$(LOCAL)/include \
	-L$(LOCAL)/lib \


# build .o files
%add.o: gsl_example.c
	$(CC) -c -o $@ $< $(CFLAGS)
# compile program
all: $(PROGRAM)
$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
