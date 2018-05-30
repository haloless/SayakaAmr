#pragma once

#include <cassert>
#include <cstdlib>

#include <vector>
#include <algorithm>

#include "SayakaCommons.h"
#include "SayakaBox.h"
#include "SayakaBlockData.h"
#include "SayakaTree.h"


SAYAKA_NS_BEGIN;



/**
 * The mapping between coarse/fine grid index
 */
struct IndexMapping {
	enum {
		// Used to shift negative fine index to be positive.
		// Its value should be larger than max ghost layer number,
		// also it must be even number!
		// So 64 shall be enough.
		LargeIndex = 64, 
		HalfLargeIndex = LargeIndex/2,
	};


	// 1D, coarse -> fine, i_f = 
	static inline int CrseToFine(int ic, int if_ind, int ic_off) {
		return (ic-ic_off)*2 + if_ind;
	}

	// 1D, fine -> coarse, found by the following cases
	// if + large = (ic-ioff)*2 + large
	// if + large = (ic-ioff)*2 + 1 + large
	static inline int FineToCrse(int ifine, int ioff) {
		assert(ifine+LargeIndex >= 0);
		return (ifine+LargeIndex)/2 - LargeIndex/2 + ioff;
	}

	// fine side
	// if + large = (ic-ioff)*2 + offset + large
	static inline int FineSide(int ifine) {
		assert(ifine+LargeIndex >= 0);
		return (ifine+LargeIndex)%2;
	}

	//
	template<int DIR>
	static inline int DirCrseToFine(int ic, int if_ind, int ic_off) {
		if (DIR >= NDIM) {
			return ic;
		} else {
			return CrseToFine(ic, if_ind, ic_off);
		}
	}
	template<int DIR>
	static inline int DirFineToCrse(int ifine, int ioff) {
		if (DIR >= NDIM) {
			return ifine;
		} else {
			return FineToCrse(ifine, ioff);
		}
	}
	template<int DIR>
	static inline int DirFineSide(int ifine) {
		if (DIR >= NDIM) {
			return 0;
		} else {
			return FineSide(ifine);
		}
	}
}; // class_indexmapping


/**
 * Hierarchical data defined on a tree.
 * It holds its own block data,
 * while the tree structure must be found from the tree.
 */
class TreeData
{
protected:

	AmrTree *m_tree; // pointer to tree

	IndexBox m_validBox;
	IndexBox m_indexBox;

	int m_numComp; // number of components
	int m_numGrow; // number of ghosts

	//
	using value_type = double;

	GridData<value_type> *m_data;
	std::vector<value_type> m_buf;

public:

	// construct empty
	TreeData();

	// construct by TREE and BOX
	TreeData(AmrTree &tree, 
		const IndexBox &validBox, 
		const IndexBox &indexBox, 
		int ngrow,
		int ncomp = 1);

	// copy
	TreeData(const TreeData &src);
	
	// move
	TreeData(TreeData &&src);

	// destroy
	~TreeData();

	// define empty data
	void define(AmrTree &tree,
		const IndexBox &validBox,
		const IndexBox &indexBox,
		int ngrow, int ncomp);

	// TODO copy assign
	TreeData& operator=(const TreeData &rhs) = delete;

	// move assign
	TreeData& operator=(TreeData &&rhs);


protected:
	
	// allocate data blocks
	// actual data buffer is allocated as a whole chunk
	// and divided to each block
	void allocData(int maxblock);
	
	void freeData();

public:
	//
	const AmrTree& tree() const { return *m_tree; }
	AmrTree& tree() { return *m_tree; }

	int maxBlockNum() const { return m_tree->maxBlockNum; }
	int blockNum() const { return m_tree->numBlocks; }

	const IndexBox& indexBox() const { return m_indexBox; }
	const IndexBox& validBox() const { return m_validBox; }

	int numComp() const { return m_numComp; }
	int numGrow() const { return m_numGrow; }

	VariableLocation indexType() const { return m_validBox.type(); }

	int isCellData() const { return m_validBox.isCellBox(); }
	int isFaceData() const { return m_validBox.isFaceBox(); }
	int isFaceData(int dir) const { return m_validBox.isFaceBox(dir); }
	int isNodeData() const { return m_validBox.isNodeBox(); }

	int getFaceDataDir() const {
		int dir = validBox().type().getFaceVarDir();
		assert(0<=dir && dir<NDIM);
		assert(isFaceData(dir));
		return dir;
	}

	//
	const DoubleBlockData& operator[](int iblock) const { return m_data[iblock]; }
	DoubleBlockData& operator[](int iblock) { return m_data[iblock]; }
	//
	const AmrTreeNode& block(int i) const { return tree()[i]; }
	AmrTreeNode& block(int i) { return tree()[i]; }

	//
	void setValue(double val) {
		const int nblock = blockNum();
		for (int i=0; i<nblock; i++) {
			m_data[i].setValue(val);
		}
	}
	void fillValue(double val) {
		const int nblock = maxBlockNum();
		for (int i=0; i<nblock; i++) {
			m_data[i].setValue(val);
		}
	}

	void copyValue(const TreeData &src) {
		assert(this->indexBox() == src.indexBox());
		assert(this->blockNum() == src.blockNum());
		assert(this->numComp() == src.numComp());

		const int nblock = this->blockNum();
		for (int i=0; i<nblock; i++) {
			m_data[i].setValue(src.m_data[i]);
		}
	}

	// Simply swap internal data pointer
	void swapValue(TreeData &rhs) {
		assert(this->indexBox() == rhs.indexBox());
		assert(this->blockNum() == rhs.blockNum());
		assert(this->numComp() == rhs.numComp());
		assert(this->numGrow() == rhs.numGrow());

		std::swap(this->m_data, rhs.m_data);
		std::swap(this->m_buf, rhs.m_buf);
	}

	//
	//void fillBoundary(int nlayer);

public:

	/*
	 * Easy creation functions.
	 * You MUST delete the data by yourself!!!
	 */
	static TreeData* CreateCellData(const AmrTree &tree, int ncomp, int ngrow);
	static TreeData* CreateFaceData(const AmrTree &tree, int dir, int ncomp, int ngrow);
	static TreeData* CreateNodeData(const AmrTree &tree, int ncomp, int ngrow);
	//
	static std::vector<TreeData*> CreateFaceDataPArray(const AmrTree &tree, int ncomp, int ngrow);
	static void ReleaseDataPArray(std::vector<TreeData*> &parray, bool clear_array=true);

	//
	// move-based version
	//
	static TreeData MakeCellData(const AmrTree &tree, int ncomp, int ngrow);
	static TreeData MakeFaceData(const AmrTree &tree, int dir, int ncomp, int ngrow);
	static TreeData MakeNodeData(const AmrTree &tree, int ncomp, int ngrow);


	//
	static void Copy(TreeData &dst, const TreeData& src,
		int dcomp, int scomp, int ncomp, int ngrow);
	static void SetValue(TreeData &dst, double val,
		int dcomp, int ncomp, int ngrow);

	/*
	 * The following functions act on a list of blocks.
	 * The given range is close/open, i.e. [ibegin, iend)
	 */
	
	//
	static void Copy(TreeData &dst, const TreeData &src, 
		int dcomp, int scomp, int ncomp, int ngrow,
		const int iblocks[], const int ibegin, const int iend);
	static inline void Copy(TreeData &dst, const TreeData &src, 
		int dcomp, int scomp, int ncomp, int ngrow,
		const std::vector<int> &iblocks, const int ibegin, const int iend)
	{
		Copy(dst, src, dcomp, scomp, ncomp, ngrow, 
			&iblocks[0], ibegin, iend);
	}
	static inline void Copy(TreeData &dst, const TreeData &src, 
		int dcomp, int scomp, int ncomp, int ngrow,
		const std::vector<int> &iblocks)
	{
		Copy(dst, src, dcomp, scomp, ncomp, ngrow,
			iblocks, 0, iblocks.size());
	}

	//
	static void SetValue(TreeData &dst, double val,
		int dcomp, int ncomp, int ngrow,
		const int iblocks[], const int ibegin, const int iend);
	static inline void SetValue(TreeData &dst, double val,
		int dcomp, int ncomp, int ngrow,
		const std::vector<int> &iblocks, const int ibegin, const int iend)
	{
		SetValue(dst, val, dcomp, ncomp, ngrow, 
			&iblocks[0], ibegin, iend);
	}

	//
	static void AddEqual(TreeData &dst, const TreeData &src, double coef,
		int dcomp, int scomp, int ncomp, int ngrow,
		const int iblocks[], const int ibegin, const int iend);
	static inline void AddEqual(TreeData &dst, const TreeData &src, double coef,
		int dcomp, int scomp, int ncomp, int ngrow,
		const std::vector<int> &iblocks)
	{
		AddEqual(dst, src, coef, dcomp, scomp, ncomp, ngrow,
			&iblocks[0], 0, iblocks.size());
	}

public:
	// write all leaf data
	int writeTreeLeafBlockData(const char *filename, int step, double time) const;
	// write specified block data, with ghost cells
	int writeBlockData(int iblock, const char *filename, int step, double time) const;
	// 
	int writeTreeLevelBlockData(int ilevel, const char *filename, int step, double time) const;
protected:
	void writeBlockZone(FILE *fp, int iblock, int step, double time, 
		int writeGrid) const;
}; // class_treedata



SAYAKA_NS_END;


