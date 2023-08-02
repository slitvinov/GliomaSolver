.SUFFIXES:
.SUFFIXES: .cpp
.SUFFIXES: .o

CXXFLAGS = -O3
CXXFLAGS_VTK = -IVTK/include/vtk-5.2
LIBS_VTK = -LVTK/lib/vtk-5.2
DEFS = -D_BLOCKSIZE_=16\
	-D_BPD_=8\
	-D_DIM=3\
	-D_FMMSILENT\
	-D_MAXLEVEL_=4\
	-DNDEBUG\
	-D_RESJUMP_=1\

LIBS = \
	$(LIBS_VTK)\
	-lvtkCommon \
	-lvtkDICOMParser\
	-lvtkexpat \
	-lvtkFiltering\
	-lvtkftgl\
	-lvtkIO\
	-lvtkjpeg\
	-lvtkmetaio\
	-lvtkNetCDF\
	-lvtkpng\
	-lvtksqlite\
	-lvtksys\
	-lvtktiff\
	-lvtkverdict\
	-lvtkzlib\

VPATH = MRAG MRAG/MRAGcore

OBJECTS = \
main.o \
MRAGBoundaryBlockInfo.o \
MRAGProfiler.o \
MRAGWavelets_StaticData.o \

brain: $(OBJECTS); $(CXX) $(OBJECTS) -o $@ $(LIBS)
.cpp.o:
	$(CXX) $< -o $@ -c $(DEFS) -Wno-deprecated $(CXXFLAGS_VTK) -IMRAG $(CXXFLAGS)
clean:; rm -rf $(OBJECTS) brain
