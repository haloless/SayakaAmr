#pragma once

#include <cstdlib>
#include <cassert>

#include <vector>
#include <map>
#include <algorithm>
#include <functional>

#include "log.h"

#include "SayakaCommons.h"
#include "SayakaBox.h"
#include "SayakaSurroundingIndex.h"
#include "SayakaTree.h"
#include "SayakaTreeData.h"

SAYAKA_NS_BEGIN;

/**
 *
 */
class TreeDataAlgorithm
{
protected:
	AmrTree &tree;
	TreeData &tree_data;
	VariableLocation var_loc;

public:
	TreeDataAlgorithm(
		AmrTree &tree_in, 
		TreeData &data_in, 
		VariableLocation varloc_in)
		: tree(tree_in),
		tree_data(data_in),
		var_loc(varloc_in)
	{}

	virtual ~TreeDataAlgorithm() {}


}; // 


//
// BCType is the numerical type when treating BC for the variable.
// - Neumann
// - Dirichlet
// - SimpleFill
// 
enum BCType {

	Neumann = 0, 
	Dirichlet,
	SimpleFill, // directly set value in ghost cell
};

//
// BCRegister is used to keep 
//
struct BCRegister {
	typedef std::map<int,int> bcmap_type;
	typedef std::map<int,double> bcval_type;

	using bcfun = std::function<double(*)(const double x[], const double t)>;

	bcmap_type bcmaps[SurroundingIndex::NumSurr];
	bcval_type bcvals[SurroundingIndex::NumSurr];
	

	//void fillBCMap(int phys_bc, int math_bc)

	void setBCMap(int isurr, int jsurr, int ksurr, 
		int phys_bc, int math_bc)
	{
		SurroundingIndex ind(isurr, jsurr, ksurr);
		bcmaps[ind][phys_bc] = math_bc;
		
		// init BC value to homogeneous
		bcvals[ind][phys_bc] = 0.0;
	}
	void setBCVal(int isurr, int jsurr, int ksurr,
		int phys_bc, double bc_val)
	{
		SurroundingIndex ind(isurr, jsurr, ksurr);
		bcvals[ind][phys_bc] = bc_val;
	}

	void setBCMap(const SurroundingIndex &isurr, int phys_bc, int math_bc) {
		bcmaps[isurr][phys_bc] = math_bc;
		bcvals[isurr][phys_bc] = 0.0;
	}

	void setBCVal(const SurroundingIndex &isurr, int phys_bc, double bc_val) {
		bcvals[isurr][phys_bc] = bc_val;
	}

	int getBCMap(int isurr, int jsurr, int ksurr, 
		int phys_bc) const
	{
		// TODO
		SurroundingIndex ind(isurr, jsurr, ksurr);
		auto it = bcmaps[ind].find(phys_bc);
		if (it == bcmaps[ind].end()) {
			LOGPRINTF("Unknown BC mapping for phys_bc=%d\n", phys_bc);
			exit(1);
		} 
		return it->second;
	}
	double getBCVal(int isurr, int jsurr, int ksurr, 
		int phys_bc) const
	{
		SurroundingIndex ind(isurr, jsurr, ksurr);
		bcval_type::const_iterator it = bcvals[ind].find(phys_bc);
		if (it == bcvals[ind].end()) {
			LOGPRINTF("Unknown BC mapping for phys_bc=%d\n", phys_bc);
			exit(1);
		} 
		return it->second;
	}
};


struct BlockBCRecord
{
protected:

	int bc_type[FaceIndex::NumFace];
	double bc_val[FaceIndex::NumFace];



public:

	void setInteriorBC(int face) {
		assert(0<=face && face<FaceIndex::NumFace);
		bc_type[face] = -1;
		bc_val[face] = 0;
	}
	void setInteriorBC() {
		for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
			setInteriorBC(face);
		}
	}

	void setBlock(const AmrTreeNode &block, const BCRegister &bcreg) {
		for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
			SurroundingIndex surr = SurroundingIndex::FaceSurrounding(face);
			int isurr = surr.ipos();
			int jsurr = surr.jpos();
			int ksurr = surr.kpos();
			
			int physbc = block.neighbor[face];
			if (physbc <= NeighborType::PhysBndry) {
				// block face is physical BC
				int mathbc = bcreg.getBCMap(isurr,jsurr,ksurr, physbc);
				double bcval = bcreg.getBCVal(isurr,jsurr,ksurr, physbc);
				bc_type[face] = mathbc;
				bc_val[face] = bcval;
			} else {
				// interior block face
				setInteriorBC(face);
			}
		}
	}

	//
	const int& type(int face) const { return bc_type[face]; }
	int& type(int face) { return bc_type[face]; }
	//
	const double& value(int face) const { return bc_val[face]; }
	double& value(int face) { return bc_val[face]; }
};


// 
// BCPatch is used to fill physical boundary.
// 
// 
class BCPatch
{
protected:

	AmrTree * m_tree; // the underlying tree

	VariableLocation m_varloc; // cell/face/node type

	int m_ncomp;
	int m_ngrow;

	IndexBox bndryBoxes[3][3][3];

	std::vector<BCRegister> bndryRegs;

public:

	// empty, undefined BC
	BCPatch();

	__declspec(deprecated)
	BCPatch(
		AmrTree &tree_in, 
		TreeData &data_in, 
		VariableLocation varloc_in);

	void define(
		AmrTree &tree, 
		VariableLocation varloc, 
		int ncomp, int ngrow);

	void defineByData(TreeData &data) {
		define(data.tree(), data.indexType(), data.numComp(), data.numGrow());
	}

	//----------------------------------------

	// 
	const BCRegister& boundaryRegister(int comp) const { 
		//assert(0<=comp && comp<tree_data.numComp());
		assert(0 <= comp && comp < m_ncomp);
		return bndryRegs[comp];
	}
	BCRegister& boundaryRegister(int comp) { 
		//assert(0<=comp && comp<tree_data.numComp());
		assert(0 <= comp && comp < m_ncomp);
		return bndryRegs[comp];
	}

	//
	// Fill physical boundary for given block.
	// 
	//
	void fillBlockBC(int physbc, 
		int isweep, int sweepdir,
		int xbndry, int ybndry, int zbndry,
		int iblock, DoubleBlockData &data, 
		int scomp, int ncomp) const;

	//void fillCellBlockBC(int bctype, int xbndry, int ybndry, int zbndry,
	//	int iblock, DoubleBlockData &blockdata) const;

protected:

	// cell-centered
	void fillBasicCellBlockBC(
		int bctype, double bcval,
		int isweep, int sweepdir,
		int xbndry, int ybndry, int zbndry,
		int iblock, DoubleBlockData &blockdata, int comp) const;
	
	// 
	void fillBasicFacexBlockBC(
		int bctype, double bcval,
		int xbndry, int ybndry, int zbndry,
		int iblock, DoubleBlockData &blockdata, int comp) const;
	void fillBasicFaceyBlockBC(
		int bctype, double bcval,
		int xbndry, int ybndry, int zbndry,
		int iblock, DoubleBlockData &blockdata, int comp) const;
	void fillBasicFacezBlockBC(
		int bctype, double bcval,
		int xbndry, int ybndry, int zbndry,
		int iblock, DoubleBlockData &blockdata, int comp) const;
public:

	const IndexBox& bndryBox(int ibndry, int jbndry, int kbndry) const { 
		assert(-XDIM<=ibndry && ibndry<=XDIM);
		assert(-YDIM<=jbndry && jbndry<=YDIM);
		assert(-ZDIM<=kbndry && kbndry<=ZDIM);
		return bndryBoxes[ibndry+1][jbndry+1][kbndry+1];
	}

	void calcRanges();
}; // class_boundaryconditionpatch



SAYAKA_NS_END;


