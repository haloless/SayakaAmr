

#include "SayakaFaceIndex.h"



namespace sayaka 
{


//
FaceIndex FaceIndex::opposite() const {
	const static int opp[] = {
		1, 0, 
		3, 2, 
		5, 4,
	};
	return opp[m_index]; 
}

} // namespace_sayaka




