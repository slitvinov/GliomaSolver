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
#include "Glioma_UQ_DataPreprocessing.h"
#include "Glioma_ComputePFF_CahnHilliard.h"

#include "dat2VP.h"

#ifdef HYPRE
#include <mpi.h>
#include "Glioma_BrainDeformation.h"
#include "../Tests/HelmholtzTest.h"
#endif


using namespace std;
using namespace MRAG;

int main(int argc,const char ** argv)
{
    
#ifdef HYPRE
    MPI_Init(&argc, (char ***)&argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if(rank==0){
        printf("\n\n       MRAG Launched         \n\n");
        printf("Running with %d MPI processes: \n", size);
    }
    
#else
    printf("\n\n       MRAG Launched         \n\n");
#endif
    
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

#ifdef HYPRE
    MPI_Finalize();
    delete s;
    
    if(rank==0) printf("we spent: %2.2f sec \n",(t1-t0).seconds());
    if(rank==0) std::cout << std::endl << "MRAG Terminated" << std::endl;

#else
    delete s;
    printf("we spent: %2.2f sec \n",(t1-t0).seconds());
    std::cout << std::endl << "MRAG Terminated" << std::endl;
#endif
    
}



