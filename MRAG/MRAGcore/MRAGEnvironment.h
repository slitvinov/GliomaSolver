/*
 *  MRAGEnvironment.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/21/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
namespace MRAG
{
	
//OS STUFF
#define _MRAG_OS_WINDOWS _WIN32
#define _MRAG_OS_APPLE __APPLE__

#ifdef _WIN32
#define _MRAG_OS _MRAG_OS_WINDOWS
#else 
#define _MRAG_OS _MRAG_OS_APPLE
#endif

}


#pragma once
namespace MRAG
{
	//ENVIRONMENT: RT SETUP
	namespace Environment
	{
        /**
         * General setup of MRAG Environment. Should be called before doing stuff.
         */
		inline void setup(int threads=-1)
		{
		}
	}
}
