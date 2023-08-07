.SUFFIXES:
.SUFFIXES: .cpp
.SUFFIXES: .o

CXXFLAGS = -O3
DEFS = -D_BLOCKSIZE_=16\
	-D_BPD_=16\
	-D_FMMSILENT\
	-D_MAXLEVEL_=4\
	-D_RESJUMP_=1\

VPATH = MRAG MRAG/MRAGcore

OBJECTS =\
main.o\
MRAGBoundaryBlockInfo.o\
MRAGWavelets_StaticData.o\

brain: $(OBJECTS); $(CXX) $(OBJECTS) -o $@ $(LDFLAGS)
.cpp.o:
	$(CXX) $< -o $@ -c $(DEFS) -Wno-deprecated -IMRAG $(CXXFLAGS)
clean:; rm -rf $(OBJECTS) brain
