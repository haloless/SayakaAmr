#pragma once

#include <vector>
#include <map>
#include <algorithm>

#include "config.h"
#include "log.h"
#include "vector3d.h"

#include "SayakaCommons.h"
#include "SayakaFaceIndex.h"
#include "SayakaChildIndex.h"
#include "SayakaSurroundingIndex.h"
#include "SayakaBox.h"


namespace sayaka
{



//
const int MAX_FLAG = 1;
const int MAX_LEVEL = 16;



enum TreeNodeType {
	TreeNodeType_Invalid = 0,
	TreeNodeType_IsLeaf = 1, 
	TreeNodeType_HasLeafChild = 2,
	TreeNodeType_NoLeafChild = 3,
};


/**
 * Physical boundary <= this value
 */
enum NeighborType {
	NeighborType_Boundary = -20,
};


// forward declare
struct TreeData;
struct SurroundingBlocks;
struct TreeStateData;
struct InterpPatch;
struct FillPatch;
struct BoundaryConditionPatch;


struct AmrTreeNode
{
	//
	int refineLevel;

	//
	int neighbor[MAX_FACE]; 

	//
	int child[MAX_CHILD];
	int parent;
	int whichChild; // child index in its parent

	// see TreeNodeType
	int nodeType; 
	//
	int nodeType2[MAX_FACE];

	//int empty; 
	int blockFlags[MAX_FLAG];

	//int newChild;
	//int refineFlag;
	//int coarsenFlag;
	//int stay;

	// physical size
	Vector3d blockCenter;
	Vector3d blockLength;
	RealBox boundBox;


public:
	int getLevel() const { return refineLevel; }
	void setLevel(int level) { refineLevel = level; }

	int getNodeType() const { return nodeType; }
	void setNodeType(int type) { nodeType = type; }
	bool isLeaf() const { return nodeType == TreeNodeType_IsLeaf; }
	bool isNotLeaf() const { return !isLeaf(); }
	bool isParentOfLeafChild() const { return nodeType == TreeNodeType_HasLeafChild; }

	//
	bool hasNeighbor(int iface) const { return neighbor[iface] >= 0; }
	void setNeighbor(int iface, int ineigh) { neighbor[iface] = ineigh; }
	int getNeighbor(int iface) const { return neighbor[iface]; }
	const int* getNeighbors() const { return neighbor; }
	int getChild(int ichild) const { return child[ichild]; }
	const int* getChildren() const { return child; }
	int getParent() const { return parent; }
	void setParent(int iparent) { parent = iparent; }

	void fillBlockFlags(int flag) {
		std::fill(blockFlags, blockFlags+MAX_FLAG, flag);
	}

	//
	int hasChildren() const {
		// test the first child will be sufficient
		return child[0] >= 0;
	}
	int hasChild(int ichild) const {
		return child[ichild] >= 0;
	}

	int hasParent() const { return parent >= 0; }

	//
	void setUnused() {
		refineLevel = -1;
		
		std::fill(neighbor, neighbor+MAX_FACE, -1);

		std::fill(child, child+MAX_CHILD, -1);
		parent = -1;
		whichChild = -1;

		nodeType = -1;
		std::fill(nodeType2, nodeType2+MAX_FACE, -1);

		fillBlockFlags(-1);

		blockCenter.setZero();
		blockLength.setZero();
		boundBox.vlo.setZero();
		boundBox.vhi.setZero();
	}

	void setBoundBox(const RealBox &bbox) {
		this->boundBox = bbox;
		this->blockCenter = bbox.center();
		this->blockLength = bbox.length();
	}
}; // struct_amrtreenode

struct AmrTreeLevel
{
	int level;

	//Vector3d cellSize;

	int numLevelBlock;
	std::vector<int> levelBlock;

public:
	//
	//int numLevelBlock() const { return levelBlock.size(); }
	bool isEmptyLevel() const { 
		assert(numLevelBlock>=0);
		return numLevelBlock == 0; 
	}

	typedef std::vector<int>::const_iterator level_block_iter;

	const int& operator[](int igrid) const { 
		assert(0<=igrid && igrid<numLevelBlock); 
		return levelBlock[igrid];
	}

	void addBlock(int iblock) {
		levelBlock.push_back(iblock);
		numLevelBlock += 1;
	}

	void clear() {
		numLevelBlock = 0;
		levelBlock.clear();
	}
}; // struct_amrtreelevel


/**
 *
 */
struct AmrTree
{

	int maxBlockNum;
	int maxTreeBlockNum;

	int numBlocks;
	int numBlocksNew;

	//RealBox probBound;
	//Vector3i probDivNum;

	// cell-centered index
	IndexBox blockValidIndex; // cell-centered, valid cells
	IndexBox blockGrownIndex; // cell-centered, valid & ghost cells
	int blockNumGrow;

	//double levelCellSize[MAX_LEVEL][MAX_DIM];

	int maxAmrLevel;
	int minAmrLevel;

	//
	AmrTreeNode *blocks;

	//
	int *refineFlag;
	int *coarsenFlag;
	int *newChildFlag;
	int *stayFlag;

	//
	SurroundingBlocks *surrBlocks;

	int gridChanged;

	//
	std::vector<AmrTreeLevel> levels;

	// references to data that will be dynamically processed 
	// during AMR regrid
	std::vector<TreeStateData> treeStateData;

public:
	
	void init();

	// call this after root level is setup
	void initRefineToMinLevel();

	// 
	void regrid();
	void checkCoarsen();
	void refineBlocks();
	int coarsenBlocks(int oldBlockNum);
	void setNodeType();
	void buildSurroundingBlocks(int buildCorner);
	void initNewBlockDataAfterRegrid();


	//
	int writeTreeLeafState(const char *filename, int step, double time) const;
	int writeTreeLeafBlockGrid(const char *filename, int step, double time) const;

	//
	//const AmrTreeNode& getBlock(int iblock) const { return blocks[iblock]; }
	//AmrTreeNode& getBlock(int iblock) { return blocks[iblock]; }
	const AmrTreeNode& operator[](int iblock) const { return blocks[iblock]; }
	AmrTreeNode& operator[](int iblock) { return blocks[iblock]; }

	const IndexBox& validBlockCellBox(int iblock=0) const { return blockValidIndex; }
	const IndexBox& grownBlockCellBox(int iblock=0) const { return blockGrownIndex; }
	//const int

	// TODO move these functions to each level
	Vector3d getBlockCellSize(int iblock) const {
		const IndexBox &vbox = this->validBlockCellBox();
		Vector3d size = blocks[iblock].blockLength;
		for (int dir=0; dir<NDIM; dir++) {
			size(dir) /= vbox.size(dir);
		}
		return size;
	}
	double getBlockCellVolume(int iblock) const {
		Vector3d dx = getBlockCellSize(iblock);
		double vol = 1.0;
		for (int dir=0; dir<NDIM; dir++) {
			vol *= dx(dir);
		}
		return vol;
	}
	Vector3d getBlockCellArea(int iblock) const {
		Vector3d dx = getBlockCellSize(iblock);
		double vol = getBlockCellVolume(iblock);
		Vector3d ds = Vector3d::VecMake(1.0, 1.0, 1.0);
		for (int dir=0; dir<NDIM; dir++) {
			ds(dir) = vol / dx(dir);
		}
		return ds;
	}

	int maxRefineLevel() const { return maxAmrLevel; }
	int minRefineLevel() const { return minAmrLevel; }
	int currentFinestLevel() const {
		int finest_level = 0;
		for (int i=0; i<numBlocks; i++) {
			if (blocks[i].refineLevel > finest_level) {
				finest_level = blocks[i].refineLevel;
			}
		}
		assert(minAmrLevel<=finest_level && finest_level<=maxAmrLevel);
		return finest_level;
	}

	const AmrTreeLevel& getTreeLevel(int ilevel) const { return levels[ilevel]; }
	AmrTreeLevel& getTreeLevel(int ilevel) { return levels[ilevel]; }

	void cacheTreeLevels();
	void cacheTreeLevel(int ilevel);

	//
	void fillFlags(int *flag, int value) const {
		for (int i=0; i<maxTreeBlockNum; i++) {
			flag[i] = value;
		}
	}


	//

}; // class_amrtree


struct SurroundingBlocks
{
	int blocks[3][3][3];
	int types[3][3][3];
	int bcrefs[3][3][3];


public:
	// index [-1, 0, +1]
	const int& operator()(int ii, int jj, int kk) const {
		assert(-1<=ii && ii<=1);
		assert(-1<=jj && jj<=1);
		assert(-1<=kk && kk<=1);
		return blocks[ii+1][jj+1][kk+1];
	}
	int& operator()(int ii, int jj, int kk) {
		assert(-1<=ii && ii<=1);
		assert(-1<=jj && jj<=1);
		assert(-1<=kk && kk<=1);
		return blocks[ii+1][jj+1][kk+1];
	}

	const int& type(int ii, int jj, int kk) const {
		assert(-1<=ii && ii<=1);
		assert(-1<=jj && jj<=1);
		assert(-1<=kk && kk<=1);
		return types[ii+1][jj+1][kk+1];
	}
	int& type(int ii, int jj, int kk) {
		assert(-1<=ii && ii<=1);
		assert(-1<=jj && jj<=1);
		assert(-1<=kk && kk<=1);
		return types[ii+1][jj+1][kk+1];
	}

	const int& bcref(int ii, int jj, int kk) const {
		return bcrefs[ii+1][jj+1][kk+1];
	}
	int& bcref(int ii, int jj, int kk) {
		return bcrefs[ii+1][jj+1][kk+1];
	}

	int isCoarseLevelSurrounding(int ii, int jj, int kk) const {
		int surr = (*this)(ii,jj,kk);
		return surr<0 && surr>NeighborType_Boundary;
	}
	int isExternalBoundarySurrounding(int ii, int jj, int kk) const {
		int surr = (*this)(ii,jj,kk);
		return surr <= NeighborType_Boundary;
	}
	int isSameLevelSurrounding(int ii, int jj, int kk) const {
		int surr = (*this)(ii,jj,kk);
		return surr >= 0;
	}
	int hasCoarseLevelSurrounding() const {
		for (int kk=-ZDIM; kk<=ZDIM; kk++) {
			for (int jj=-YDIM; jj<=YDIM; jj++) {
				for (int ii=-XDIM; ii<=XDIM; ii++) {
					if (isCoarseLevelSurrounding(ii,jj,kk)) {
						assert((*this)(ii,jj,kk) == -1);
						return 1;
					}
				}
			}
		}
		return 0;
	}
	int hasExternalBoundarySurrounding() const {
		for (int kk=-ZDIM; kk<=ZDIM; kk++) {
			for (int jj=-YDIM; jj<=YDIM; jj++) {
				for (int ii=-XDIM; ii<=XDIM; ii++) {
					if (isExternalBoundarySurrounding(ii,jj,kk)) {
						return 1;
					}
				}
			}
		}
		return 0;
	}

	void setUnused() {
		for (int kk=0; kk<3; kk++) {
			for (int jj=0; jj<3; jj++) {
				for (int ii=0; ii<3; ii++) {
					blocks[ii][jj][kk] = -1;
					types[ii][jj][kk] = -1;
					bcrefs[ii][jj][kk] = -1;
				}
			}
		}
	}

	

	int collectNeighborOrFillBoundary(
		int cneigh[MAX_FACE],
		const AmrTree &tree, int iblock) const 
	{
		std::fill(cneigh, cneigh+MAX_FACE, -1);

		if (iblock >= 0) { // valid block, collect its neighbors
			const int *neigh = tree.blocks[iblock].neighbor;
			std::copy(neigh, neigh+MAX_FACE, cneigh);
			return 1;
		} else if (iblock <= NeighborType_Boundary) { 
			// hits boundary, fill with this boundary type
			std::fill(cneigh, cneigh+MAX_FACE, iblock);
			return 1;
		} else {
			return 0;
		}
	}

	// return (corner,bcref)
	std::pair<int,int> collectCorner(int ic, int jc, int kc, const AmrTree &tree) const {
		assert(ic==-1 || ic==1);
		assert(jc==-1 || jc==1);
		assert(kc==-1 || kc==1);

		if (NDIM != 3) {
			LOGPRINTF("%s: only for NDIM=3\n", __FUNCTION__);
			exit(1);
		}

		int iface = ic==-1 ? 0 : 1;
		int jface = jc==-1 ? 2 : 3;
		int kface = kc==-1 ? 4 : 5;

		//int corner = -1;
		//int lb;
		//int cneigh[MAX_FACE];

		//if (corner == -1) {
		//	lb = (*this)(0,jc,kc);
		//	if (collectNeighborOrFillBoundary(cneigh, tree, lb)) {
		//		corner = cneigh[iface];
		//	}
		//}
		//if (corner == -1) {
		//	lb = (*this)(ic,0,kc);
		//	if (collectNeighborOrFillBoundary(cneigh, tree, lb)) {
		//		corner = cneigh[jface];
		//	}
		//}
		//if (corner == -1) {
		//	lb = (*this)(ic,jc,0);
		//	if (collectNeighborOrFillBoundary(cneigh, tree, lb)) {
		//		corner = cneigh[kface];
		//	}
		//}

		const int lbref_boundary = 0;
		const int lbref_samelevel = 1;
		const int lbref_crselevel = 2;

		std::pair<int,int> corner = std::make_pair(-1,-1);
		int lbref = -1;
		for (int dir=0; dir<NDIM; dir++) {
			int ii = dir==0 ? 0 : ic;
			int jj = dir==1 ? 0 : jc;
			int kk = dir==2 ? 0 : kc;
			int ff = dir==0 ? iface : dir==1 ? jface : kface;

			int lb = (*this)(ii,jj,kk);
			if (lb >= 0) { // in-level true cell
				// reset to that position 
				lbref = lbref_samelevel;
				corner.first = tree[lb].neighbor[ff];
				corner.second = SurroundingIndex(ii,jj,kk);
			} else if (lb <= NeighborType_Boundary) { // boundary
				// check no cell position yet
				if (lbref == -1) {
					lbref = lbref_boundary;
					corner.first = lb;
					corner.second = SurroundingIndex(ii,jj,kk);
				}
			} else { // coarse-fine cell
				// set to coarse-fine if no in-level reference yet
				if (lbref != lbref_samelevel) {
					lbref = lbref_crselevel;
					corner.first = -1;
					corner.second = SurroundingIndex(ii,jj,kk);
				}
			}
		}

		return corner;
	}


}; // class_surroundingblocks

struct TreeStateData
{
	// const AmrTree *m_tree;

	TreeData *old_data;
	TreeData *new_data;


	//InterpPatch *interp_patch;
	BoundaryConditionPatch *bc_patch;
	FillPatch *fill_new_patch; // bind with new_data

public:
	//TreeStateData(const AmrTree &tree) 
	//	: m_tree(&tree),
	//	old_data(NULL), new_data(NULL)
	//{}
	TreeStateData()
		: old_data(NULL), new_data(NULL), 
		fill_new_patch(NULL), 
		bc_patch(NULL)
	{}
	~TreeStateData() {
		//if (old_data) delete old_data;
		//if (new_data) delete new_data;
	}

	bool hasCurrData() const { return new_data != NULL; }
	bool hasPrevData() const { return old_data != NULL; }
	TreeData& currData() { assert(new_data); return *new_data; }
	TreeData& prevData() { assert(old_data); return *old_data; }
	const TreeData& currData() const { assert(new_data); return *new_data; }
	const TreeData& prevData() const { assert(old_data); return *old_data; }

	//
	void swapTimeStates();
	void overwritePrevTimeState();
}; // struct_treestatedata



} // namespace sayaka



