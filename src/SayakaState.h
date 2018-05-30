#pragma once




enum SayakaCellStateComponents {
	SAYAK_CELLSTATE_COMP_U,
	SAYAK_CELLSTATE_COMP_V,
	SAYAK_CELLSTATE_COMP_W,


	SAYAKA_CELLSTATE_NCOMP,
};




struct SayakaFaceStateVector {
	double un;
	double var;
	//double value[];
};

struct SayakaStateVector {
	/// face state
	SayakaFaceStateVector f[SAYAKA_SPACEDIM*2];
	
	// cell state
	double data[SAYAKA_CELLSTATE_NCOMP];
};




