.SUFFIXES:
.SUFFIXES: .cpp
.SUFFIXES: .o

CXXFLAGS = -O3
DEFS = -D_BLOCKSIZE_=8
VPATH = MRAG MRAG/MRAGcore

OBJECTS =\
main.o\
lib.o\
MRAGBoundaryBlockInfo.o\
MRAGWavelets_StaticData.o\

brain: $(OBJECTS); $(CXX) $(OBJECTS) -o $@ $(LDFLAGS)
.cpp.o:
	$(CXX) $< -o $@ -c $(DEFS) -IMRAG $(CXXFLAGS)
clean:; rm -rf $(OBJECTS) brain
main.o: write.h lib.h
lib.o: lib.h
