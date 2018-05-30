#pragma once

#include "SayakaCommons.h"
#include "SayakaTree.h"
#include "SayakaTreeData.h"


SAYAKA_NS_BEGIN;


/**
 * The type of block face.
 * 
 */
enum TowerBlockFaceType {
	FACE_FINE_FINE = 0,
	FACE_FINE_BC,
	FACE_FINE_CRSE,
	FACE_CRSE_FINE,
};

/**
 * If face_type is fine_crse, 
 * then face_neigh holds the coarse neighbor block.
 *
 * If face_type is crse_fine,
 * face_neigh holds the coarse neighbor at same coarse level;
 * crse_fine_neigh holds the fine neighbors at fine level.
 */
struct BlockFaceRegister
{
	int face_type;
	int face_neigh;

	Vector3i fine_crse_offset;

	int crse_fine_neigh[ChildIndex::NumChildOnFace];
	Vector3i crse_fine_offset[ChildIndex::NumChildOnFace];
	IndexBox crse_fine_subbox[ChildIndex::NumChildOnFace];

	bool isFaceBC() const { return face_type == FACE_FINE_BC; }
	bool isFaceFineFine() const { return face_type == FACE_FINE_FINE; }
	bool isFaceFineCrse() const { return face_type == FACE_FINE_CRSE; }
	bool isFaceCrseFine() const { return face_type == FACE_CRSE_FINE; }
}; // class_blockfaceregister

typedef std::vector<BlockFaceRegister> BlockFaceRegArray;


/**
 * A simple container class to hold the blocks on a MG level.
 */
struct MGLevel : public AmrTreeLevel
{
	typedef AmrTreeLevel super_type;


	// [begin, end)
	std::vector<int> level_begin; // closed range
	std::vector<int> level_end; // open range

public:
	void define(int mg_level, int maxlevel) {
		assert(1<=mg_level && mg_level<=maxlevel);

		this->level = mg_level;
		super_type::clear();

		//
		level_begin.resize(maxlevel+1);
		level_end.resize(maxlevel+1);
		std::fill(level_begin.begin(), level_begin.end(), 0);
		std::fill(level_end.begin(), level_end.end(), 0);
	}

	void clearMGLevel() {
		// clear blocks cached
		super_type::clear();

		// reset begin/end ranges
		std::fill(level_begin.begin(), level_begin.end(), 0);
		std::fill(level_end.begin(), level_end.end(), 0);
	}

	int getMGLevel() const { return this->level; }
	void setMGLevel(int mg_level) { this->level = mg_level; }

	int coversTreeLevel(int ilevel) const { 
		assert(ilevel>=1); 
		return ilevel<=this->level; 
	}
	int hasTreeLevel(int ilevel) const {
		assert(ilevel>=1); 
		return coversTreeLevel(ilevel) && size(ilevel)>0;
	}

	int begin(int ilevel) const {
		assert(1<=ilevel && ilevel<=this->level);
		return level_begin[ilevel];
	}
	int end(int ilevel) const {
		assert(1<=ilevel && ilevel<=this->level);
		return level_end[ilevel];
	}
	int size(int ilevel) const {
		assert(1<=ilevel && ilevel<=this->level);
		return end(ilevel) - begin(ilevel);
	}

	void setBegin(int ilevel, int ibegin) {
		assert(1<=ilevel && ilevel<=this->level);
		level_begin[ilevel] = ibegin;
	}
	void setEnd(int ilevel, int iend) {
		assert(1<=ilevel && ilevel<=this->level);
		level_end[ilevel] = iend;
	}
};


/**
 * The tower is composed of blocks that satisfy
 * (a) block is leaf and its level<mg_level;
 * (b) block level==mg_level
 *
 * Apparently, on the finest level of tree,
 * the tower is equal to the union of all leaf blocks.
 */
class MGLevelTower : public MGLevel
{
	typedef MGLevel super_type;

public:
	enum TowerCellMask {
		//MASK_INVALID = -1,
		//MASK_INTERNAL = 0,
		//MASK_FINEFINE = (1 << FACE_FINE_FINE),
		MASK_INVALID = 0,
		MASK_INTERNAL = (1 << FACE_FINE_FINE),
		MASK_PHYSBC = (1 << FACE_FINE_BC),
		MASK_FINECRSE = (1 << FACE_FINE_CRSE),
		MASK_CRSEFINE = (1 << FACE_CRSE_FINE),
	};
	enum TowerFluxSyncMode {
		FLUX_SYNC_AVG,
		FLUX_SYNC_SUM,
	};

	// null construction
	MGLevelTower()
		: m_tree(NULL)
	{}
	// construction at given mg_level
	MGLevelTower(const AmrTree &tree, int mg_level)
		: m_tree(NULL)
	{
		define(tree, mg_level);
	}
	// 
	virtual ~MGLevelTower() 
	{}

	// attach to tree, allocate memory, etc.
	void define(const AmrTree &tree, int mg_level);

	//
	void update(bool change_to_current_finest=false);

	//
	void maskTowerCell(TreeData &mask, int maskcomp, int nlayer=1) const;

	//
	void fillTowerBndry(TreeData &inout, int dcomp, int ncomp, int nlayer,
		int fillPhysBC=1, int fillCorner=1, int usePiecewiseConstant=1);

	//
	void syncTowerCrseFineFlux(int dir, TreeData &flux, int fluxcomp,
		TowerFluxSyncMode sync_mode) const;

	// 
	void syncSubLevelCrseFineFlux(int fine_ilevel, int crse_ilevel,
		int dir, TreeData &flux, int fluxcomp,
		TowerFluxSyncMode sync_mode) const;

	// set value for all blocks in the data
	void setTowerDataValue(TreeData &data, int scomp, int ncomp, int ngrow, 
		double value) const;
	// set value for blocks on the specified level
	void setSubLevelDataValue(int ilevel, 
		TreeData &data, int scomp, int ncomp, int ngrow, 
		double value) const;


protected:
	// test and add a block
	// the block face is tested and set.
	bool registerBlock(int mg_level, int ilevel, int iblock);

public:
	const BlockFaceRegister& getFaceReg(int iblock, int iface) const {
		return face_reg[iface][iblock];
	}
	BlockFaceRegister& getFaceReg(int iblock, int iface) {
		return face_reg[iface][iblock];
	}
	void getBlockFaceRegs(int iblock, const BlockFaceRegister* blockFaceRegs[]) const {
		for (FaceIndex iface=0; iface<FaceIndex::NumFace; iface++) {
			blockFaceRegs[iface] = &getFaceReg(iblock, iface);
		}
	}

	int isSubLevelEmpty(int ilevel) const {
		assert(1<=ilevel && ilevel<=this->getMGLevel());
		return (this->size(ilevel) <= 0);
	}

	const AmrTree& getTree() const { assert(m_tree); return *m_tree; }

protected:
	const AmrTree *m_tree;

	BlockFaceRegArray face_reg[FaceIndex::NumFace];

};

typedef std::vector<MGLevelTower> MGLevelTowerArray;



SAYAKA_NS_END;



