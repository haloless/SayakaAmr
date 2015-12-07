
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


namespace sayaka 
{


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








} // namespace sayaka
