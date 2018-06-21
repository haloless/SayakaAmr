
#include <cmath>

#include <iostream>

#include "log.h"
#include "vector3d.h"

#include "SayakaCommons.h"
#include "Sayaka.h"
#include "SayakaTree.h"
#include "SayakaTreeData.h"
#include "SayakaBoundaryPatch.h"
#include "SayakaFillPatch.h"
#include "SayakaVarLocInterpolator.h"
#include "SayakaMacSolver.h"

//
//#include "Test4.h"
//#include "Test5.h"
//#include "Test6.h"
//#include "Test7.h"
#include "Test8.h"

using namespace sayaka;


int main(int argc, char *argv[] ) {

	sayaka::Initialize(argc, argv);

	//Test2 t2;
	//t2.init(argc, argv);

	//Test3 test3;
	//test3.init(argc, argv);

	//Test4 t4;
	//t4.init(argc, argv);

	//Test5 t5;
	//t5.init(argc, argv);

	//Test6 t6;
	//t6.init(argc, argv);

	//Test7 t7;
	//t7.init(argc, argv);

	Test8 t;
	t.init(argc, argv);

	return 0;
}

