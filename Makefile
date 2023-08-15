.SUFFIXES:
.SUFFIXES: .c
.SUFFIXES: .cpp
.SUFFIXES: .o

CXXFLAGS = -O3
CFLAGS = -O3
DEFS = -D_BLOCKSIZE_=8

OBJECTS =\
main.o\
lib.o\
MRAGBoundaryBlockInfo.o\
MRAGWavelets_StaticData.o\

brain: $(OBJECTS); $(CXX) $(OBJECTS) -o $@ $(LDFLAGS)
.cpp.o:
	$(CXX) $< -o $@ -c $(DEFS) $(CXXFLAGS)
.c.o:
	$(CC) $< -o $@ -c $(DEFS) $(CFLAGS)
clean:; rm -rf $(OBJECTS) brain
main.o: lib.h
lib.o: lib.h
