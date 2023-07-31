.POSIX:
CXXFLAGS = -O2 -g
S = main.cpp
likelihood: $S ArgumentParser.h  Matrix.h
	$(CXX) $(CXXFLAGS) $S -o $@
clean:; rm -f likelihood
