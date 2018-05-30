#include "SayakaSurroundingIndex.h"

SAYAKA_NS_BEGIN;

// surrounding connected to given face
SurroundingIndex SurroundingIndex::FaceSurrounding(const FaceIndex & face)
{
	int dir = face.dir();
	int side = face.side();

	int pos[MAX_DIM] = { 0 };
	pos[dir] = 2 * side - 1;
	assert(pos[dir] == -1 || pos[dir] == 1);

	SurroundingIndex face_surr(pos[0], pos[1], pos[2]);
	return face_surr;
}


SAYAKA_NS_END;

