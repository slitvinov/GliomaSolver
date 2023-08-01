#include "Glioma_ReactionDiffusion.h"
using namespace MRAG;
int main(int argc,const char ** argv)
{
    ArgumentParser parser(argc, argv);
    Environment::setup(max(1, parser("-nthreads").asInt()));
    Glioma_ReactionDiffusion s(argc, (const char **)argv);
    s.run();
}
