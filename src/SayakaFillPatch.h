#pragma once


#include "SayakaCommons.h"
#include "SayakaBox.h"
#include "SayakaTree.h"
#include "SayakaTreeData.h"
// TODO remove this dependency??
#include "SayakaInterpPatch.h"

namespace sayaka
{


/**
 * Treatment of ghost cells enclosing blocks, including
 * (a) internal in-level boundary 
 * (b) fine/coarse boundary
 * (c) physical boundary
 */
struct FillPatch : public InterpPatch
{
	/*
	 * coord_flag: 
	 * -1 for all directions;
	 */
	enum CoordFillFlag {
		COORD_ALL = -1,
		COORD_XDIR = 0,
		COORD_YDIR,
		COORD_ZDIR,
	};
	enum DiagFillFlag {
		DIAG_IGNORE = 0,
		DIAG_FILL = 1,
	};

	struct FillPatchRange {
		int nlayer;

		// cell-centered target ghost range
		// can be used to locate range
		IndexBox dstCellRanges[SurroundingIndex::NumSurr];
		// source valid range corresponding to target range
		IndexBox srcCellRanges[SurroundingIndex::NumSurr];

		// target ghost range
		// have the same index type as the data (e.g. face-centered)
		// May have overlap with valid range (e.g. for face-centered)
		IndexBox dstRanges[SurroundingIndex::NumSurr];
		// source valid range corresponding to target range
		IndexBox srcRanges[SurroundingIndex::NumSurr];

		// target ghost range
		// have the same index type as the data
		// Do NOT have overlap with the valid range 
		// (e.g. for UMAC data, boundary parts are shrinked in x-axis)
		IndexBox dstNonOverlapRanges[SurroundingIndex::NumSurr];
		IndexBox srcNonOverlapRanges[SurroundingIndex::NumSurr];


	};

	typedef InterpPatch super_type;

	BoundaryConditionPatch *bc_patch;

	// pre-build boundary range
	std::vector<FillPatchRange> range_cache;

	// pre-build child index offset
	Vector3i child_offset_cache[ChildIndex::NumChild];

public:
	
	// deprecated
	FillPatch(AmrTree &tree_in, 
		TreeData &data_in, 
		VariableLocation varloc_in)
		: super_type(tree_in, data_in, varloc_in), 
		bc_patch(NULL)
	{
		LOGPRINTF("%s: deprecated\n", __FUNCTION__);
		exit(1);
	}

	FillPatch(AmrTree &tree_in, TreeData &data_in)
		: super_type(tree_in, data_in, data_in.validBox().type()),
		bc_patch(NULL)
	{
		// possible fill-patch band width
		const int nlayer_max = data_in.numGrow();
		assert(nlayer_max >= 0);

		// cache usable stuff for fill-patch
		cacheOffsets();
		cacheRanges(nlayer_max);
	}



	// Fill ghost cells for
	// in-level, physical
	void fillBoundary(int scomp, int ncomp, int nlayer);

	// Fill blocks in specified level only
	// NOTE this must be called from coarse level to fine level
	void fillLevelBoundary(int level, int scomp, int ncomp, int nlayer,
		int onlySameLevel=0, int applyPhysBC=1, int fillCorner=1);

	// Fill given block
	void fillBlockBoundary(int iblock, int scomp, int ncomp, int nlayer,
		int onlySameLevel=0, int applyPhysBC=1, int fillCorner=1);


	// Directly set value in boundary ghost points
	void setBlockBoundaryValue(
		int iblock, /*DoubleBlockData &out, */
		int scomp, int ncomp, int nlayer, 
		double bndryValue);
	void setBlockBoundaryValue(
		int iblock, /*DoubleBlockData &out, */
		int scomp, int ncomp, int nlayer, 
		double bndryValue, const SurroundingIndex &surr_idx);

	//
	void restrictAll(int scomp, int ncomp);

	// target level is the coarse level
	// average (level+1) down to its level
	void restrictLevel(int crse_level, int scomp, int ncomp);

	// Target block must be a parent block with children.
	// Take data from its children and perform restriction on target block.
	// It may requires proper fill-patch of fine blocks.
	// And it is the caller's responsibility to do that.
	void restrictBlock(int iblock, int scomp, int ncomp);

	//
	//void prolongAllNewBlock(int scomp, int ncomp, int nlayer);

	void prolongLevel(int fine_level, int scomp, int ncomp);

	// Target block is a child with parent.
	// Take data from its parent and perform prolongation on target child block.
	// It may requires proper fill-patch of fine blocks.
	// And it is the caller's responsibility to do that.
	void prolongBlock(int iblock, int scomp, int ncomp);

	// Only for face-centered block
	// if the data is staggered in face direction,
	// set the face value from its neighbor
	void correctFaceBlockData(int iface, 
		int iblock, DoubleBlockData &dstdata,
		int jblock, const DoubleBlockData &srcdata,
		int dstcomp, int srccomp, int ncomp);

	// Only for face-centered block
	// if the data is staggered in face direction,
	// set the face value to be the average of the two adjacent blocks
	// NOTE the input face refers to the IBLOCK,
	// which is opposite in terms of the JBLOCK
	void syncFaceBlockData(int iface, 
		int iblock, DoubleBlockData &idata,
		int jblock, DoubleBlockData &jdata,
		int scomp, int ncomp);

protected:


	// for boundaries connecting in-level neighbors
	void fillBlockBndryAtSameLevel(
		int iblock, DoubleBlockData &blockdata,
		int scomp, int ncomp, int nlayer, 
		int coordflag, DiagFillFlag cornerflag);

	// boundaries on coarse/fine interface
	void fillBlockBndryFromCoarseLevel(
		int iblock, DoubleBlockData &blockdata,
		int iparent, DoubleBlockData &parentdata,
		int scomp, int ncomp, int nlayer, 
		int coordflag, DiagFillFlag cornerflag);

	// boundaries on physical boundary
	void fillBlockBndryCond(
		int iblock, DoubleBlockData &data,
		int scomp, int ncomp, int nlayer,
		int coordflag, DiagFillFlag cornerflag);


	// 
	void cacheRanges(int nlayer_max);
	void cacheOffsets();

	const FillPatchRange& getRangeCache(int nlayer) const {
		assert(0<nlayer && nlayer<=tree_data.numGrow());
		return range_cache[nlayer-1];
	}
	const Vector3i& getChildOffsetCache(int ichild) const {
		assert(0<=ichild && ichild<ChildIndex::NumChild);
		return child_offset_cache[ichild];
	}

}; // class_fillpatch






} // namespace_sayaka


