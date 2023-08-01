/*
 *  mainGlioma.cpp
 *  Untitled
 *
 *  Created by Lipkova on 9/14/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
#include <xmmintrin.h>

#include "Tests/Test.h"
#include "Glioma_ReactionDiffusion.h"
#include "dat2VP.h"

using namespace std;
using namespace MRAG;

int main(int argc,const char ** argv)
{
    
    printf("\n\n       MRAG Launched         \n\n");
    ArgumentParser parser(argc, argv);
    Environment::setup(max(1, parser("-nthreads").asInt()));
    Glioma * s = NULL;
    s = new Glioma_ReactionDiffusion(argc, (const char **)argv);
    tbb::tick_count t1,t0;
    {
        t0=tbb::tick_count::now();
        s->run();
        t1=tbb::tick_count::now();
    }
    delete s;
    printf("we spent: %2.2f sec \n",(t1-t0).seconds());
    std::cout << std::endl << "MRAG Terminated" << std::endl;
}



