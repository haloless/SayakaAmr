
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <vector>
#include <algorithm>

#include "log.h"
#include "vector3d.h"

#include "SayakaCommons.h"
#include "SayakaTree.h"
#include "SayakaTreeData.h"


SAYAKA_NS_BEGIN;


void TreeStateData::swapTimeStates() {
	assert(new_data!=NULL && old_data!=NULL);
	new_data->swapValue(*old_data);
}
void TreeStateData::overwritePrevTimeState() {
	assert(new_data!=NULL && old_data!=NULL);
	old_data->copyValue(*new_data);
}


/*
 *
 */

void AmrTree::init() {
	//if (0 != config.load_from_file("input_sayaka.ini")) {
	//	LOGPRINTF("Failed to open config file\n");
	//	exit(1);
	//}

	//

	//
	const int maxblocks = 40000;

	this->maxBlockNum = maxblocks;

	/*
	 * tree structure
	 */

	const int maxtreeblocks = maxblocks * 1;
	this->maxTreeBlockNum = maxtreeblocks;
	//
	this->numBlocks = 0;
	this->gridChanged = 0;
	//
	this->blocks = new AmrTreeNode[maxtreeblocks];
	this->refineFlag = new int[maxtreeblocks];
	this->coarsenFlag = new int[maxtreeblocks];
	this->newChildFlag = new int[maxtreeblocks];
	this->stayFlag = new int[maxtreeblocks];
	this->surrBlocks = new SurroundingBlocks[maxtreeblocks];
	// 
	for (int i=0; i<maxtreeblocks; i++) {
		blocks[i].setUnused();

		refineFlag[i] = 0;
		coarsenFlag[i] = 0;
		newChildFlag[i] = 0;
		stayFlag[i] = 1;

		surrBlocks[i].setUnused();
	}


	/*
	 * block structure
	 */
	// TODO assert(NBX%2 == 0);


	/*
	 * 
	 */
	this->levels.resize(maxAmrLevel+1);
	for (int ilevel=0; ilevel<=maxAmrLevel; ilevel++) {
		levels[ilevel].level = ilevel;
		levels[ilevel].clear();
	}

} // amrtree_init

void AmrTree::initUniformRootLevel(
	const Vector3i &ndiv,
	const RealBox &probbox,
	const int bc_type[MAX_FACE],
	const bool is_periodic[MAX_FACE])
{
	const int root_level = 1; // 

	const int n0x = ndiv(0);
	const int n0y = ndiv(1);
	const int n0z = ndiv(2);
	const int nroot = n0x * n0y * n0z;

	// define indexing
	auto ind = [=](int i, int j, int k) {
		return i + j * n0x + k * n0x * n0y;
	};

	this->numBlocks = nroot;

	for (int k = 0; k<n0z; k++) {
		for (int j = 0; j<n0y; j++) {
			for (int i = 0; i<n0x; i++) {
				// grab a node
				int iroot = ind(i, j, k);
				AmrTreeNode &root = this->blocks[iroot];

				// set physical size
				Vector3d block_lo = probbox.lo();
				Vector3d block_hi = probbox.lo();
				block_lo.x += probbox.length(0) / n0x * i;
				block_hi.x += probbox.length(0) / n0x * (i + 1);
				block_lo.y += probbox.length(1) / n0y * j;
				block_hi.y += probbox.length(1) / n0y * (j + 1);
				if (NDIM == 3) {
					block_lo.z += probbox.length(2) / n0z * k;
					block_hi.z += probbox.length(2) / n0z * (k + 1);
				}

				RealBox rootbox(block_lo, block_hi);
				root.boundBox = rootbox;
				root.blockCenter = rootbox.center();
				root.blockLength = rootbox.length();

				// set info
				root.nodeType = TreeNodeType_IsLeaf;
				root.setLevel(root_level);

				// set connections and BC
				for (FaceIndex iface = 0; iface<FaceIndex::NumFace; iface++) {
					int ineigh = -1;

					if (iface == FaceIndex::FACE_XMINUS) {
						if (i == 0) {
							ineigh = is_periodic[iface] ? ind(n0x - 1, j, k) : bc_type[iface];
						}
						else {
							ineigh = ind(i - 1, j, k);
						}
					}
					else if (iface == FaceIndex::FACE_XPLUS) {
						if (i == n0x - 1) {
							ineigh = is_periodic[iface] ? ind(0, j, k) : bc_type[iface];
						}
						else {
							ineigh = ind(i + 1, j, k);
						}
					}
					else if (iface == FaceIndex::FACE_YMINUS) {
						if (j == 0) {
							ineigh = is_periodic[iface] ? ind(i, n0y - 1, k) : bc_type[iface];
						}
						else {
							ineigh = ind(i, j - 1, k);
						}
					}
					else if (iface == FaceIndex::FACE_YPLUS) {
						if (j == n0y - 1) {
							ineigh = is_periodic[iface] ? ind(i, 0, k) : bc_type[iface];
						}
						else {
							ineigh = ind(i, j + 1, k);
						}
					}
					else if (iface == FaceIndex::FACE_ZMINUS) {
						if (k == 0) {
							ineigh = is_periodic[iface] ? ind(i, j, n0z - 1) : bc_type[iface];
						}
						else {
							ineigh = ind(i, j, k - 1);
						}
					}
					else if (iface == FaceIndex::FACE_ZPLUS) {
						if (k == n0z - 1) {
							ineigh = is_periodic[iface] ? ind(i, j, 0) : bc_type[iface];
						}
						else {
							ineigh = ind(i, j, k + 1);
						}
					}

					//assert(0<=ineigh && ineigh<n0x*n0y*n0z);
					root.neighbor[iface] = ineigh;

				} // end loop face
			}
		}
	}
}


void AmrTree::initRefineToMinLevel() {
	// before this is called, we should have setup the root level(=1)

	//{
	//	// build surrounding list
	//	const int check_corner = 1;
	//	buildSurroundingBlocks(check_corner);

	//	// build level cache for root level
	//	const int root_level = 1;
	//	cacheTreeLevel(root_level);
	//}

	//// initial static refine up to min level
	//for (int level=1; level<minAmrLevel; level++) {
	//	for (int i=0; i<numBlocks; i++) {
	//		refineFlag[i] = 1;
	//	}

	//	regrid();
	//}

	// directly perform regrid without set refine flags
	// this should exactly do the work for us
	for (int level=1; level<=minAmrLevel; level++) {
		regrid();
	}

	assert(currentFinestLevel() == minAmrLevel);
}




int AmrTree::writeTreeLeafState(const char *filename, int step, double time) const {
	if (NDIM != 2) {
		LOGPRINTF("%s: 2D only\n", __FUNCTION__);
		exit(1);
	}

	FILE *fp = fopen(filename, "w");
	if (!fp) {
		LOGPRINTF("%s: failed to open %s\n", __FUNCTION__, filename);
		return 1;
	}

	//
	fprintf(fp, "TITLE=\"Amr Leaf\"\n");
	fprintf(fp, "VARIABLES=\"X\" \"Y\" \"level\"\n");

	for (int i=0; i<numBlocks; i++) {
		if (blocks[i].isLeaf()) {
			fprintf(fp, "ZONE DATAPACKING=BLOCK, I=2, J=2, VARLOCATION=([3]=CELLCENTERED)\n");
			fprintf(fp, "STRANDID=%d, SOLUTIONTIME=%e\n", step, time);
			// x
			for (int jj=0; jj<2; jj++) {
				for (int ii=0; ii<2; ii++) {
					double xx = blocks[i].boundBox.vlo.x + blocks[i].blockLength.x * ii;
					double yy = blocks[i].boundBox.vlo.y + blocks[i].blockLength.y * jj;
					fprintf(fp, "%lf ", xx);
				}
			}
			fprintf(fp, "\n");
			// y
			for (int jj=0; jj<2; jj++) {
				for (int ii=0; ii<2; ii++) {
					double xx = blocks[i].boundBox.vlo.x + blocks[i].blockLength.x * ii;
					double yy = blocks[i].boundBox.vlo.y + blocks[i].blockLength.y * jj;
					fprintf(fp, "%lf ", yy);
				}
			}
			fprintf(fp, "\n");
			// level
			fprintf(fp, "%d", blocks[i].getLevel());
			fprintf(fp, "\n");
		}
	}

	fclose(fp);

	LOGPRINTF("%s: write %s\n", __FUNCTION__, filename);

	return 0;
}

int AmrTree::writeTreeLeafBlockGrid(const char *filename, int step, double time) const {
	if (NDIM != 2) {
		LOGPRINTF("%s: 2D only\n", __FUNCTION__);
		exit(1);
	}


	FILE *fp = fopen(filename, "w");
	if (!fp) {
		LOGPRINTF("%s: failed to open %s\n", __FUNCTION__, filename);
		return 1;
	}

	//
	fprintf(fp, "TITLE=\"Amr Grid\"\n");
	//
	fprintf(fp, "FILETYPE=GRID\n");
	fprintf(fp, "VARIABLES=\"X\" \"Y\"\n");

	const IndexBox &validBox = validBlockCellBox();
	const Vector3i vnb = validBox.size();

	for (int i=0; i<numBlocks; i++) {
		if (blocks[i].isLeaf()) {
			double xlo = blocks[i].boundBox.vlo.x;
			double ylo = blocks[i].boundBox.vlo.y;
			double zlo = blocks[i].boundBox.vlo.z;
			double dx = blocks[i].blockLength.x / vnb.x;
			double dy = blocks[i].blockLength.y / vnb.y;
			double dz = blocks[i].blockLength.z / vnb.z;

			fprintf(fp, "ZONE\n");
			fprintf(fp, "I=%d, J=%d\n", vnb.x+1, vnb.y+1);
			fprintf(fp, "DATAPACKING=BLOCK\n");
			fprintf(fp, "STRANDID=%d, SOLUTIONTIME=%e\n", step, time);
			// x
			for (int jj=0; jj<=vnb.y; jj++) {
				for (int ii=0; ii<=vnb.x; ii++) {
					double xx = xlo + dx * ii;
					double yy = ylo + dy * jj;
					fprintf(fp, "%lf ", xx);
				}
			}
			fprintf(fp, "\n");
			// y
			for (int jj=0; jj<=vnb.y; jj++) {
				for (int ii=0; ii<=vnb.x; ii++) {
					double xx = xlo + dx * ii;
					double yy = ylo + dy * jj;
					fprintf(fp, "%lf ", yy);
				}
			}
			fprintf(fp, "\n");
		}
	}

	fclose(fp);

	LOGPRINTF("%s: write %s\n", __FUNCTION__, filename);

	return 0;
}








SAYAKA_NS_END;

