/*
 *  MRAGProfiler.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/13/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */


#include "MRAGProfiler.h"
	#include <time.h>
	void MRAG::ProfileAgent::_getTime(clock_t& time)
	{
		time = clock();
	}

	float MRAG::ProfileAgent::_getElapsedTime(const clock_t& tS, const clock_t& tE)
	{
		return (tE - tS)/(double)CLOCKS_PER_SEC;
	}

	
