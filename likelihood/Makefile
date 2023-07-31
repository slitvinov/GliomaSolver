.POSIX:
CXXFLAGS = -O2 -g
S = Likelihood_main.cpp HGG_Likelihood.cpp
likelihood: $S ArgumentParser.h  HGG_Likelihood.h  Matrix.h
	$(CXX) $(CXXFLAGS) $S -o $@
clean:; rm -f likelihood
