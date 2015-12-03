#pragma once

#include <cstdlib>
#include <cassert>

#include <vector>
#include <map>
#include <algorithm>

#include "log.h"

#include "SayakaCommons.h"
#include "SayakaBox.h"
#include "SayakaSurroundingIndex.h"
#include "SayakaTree.h"
#include "SayakaTreeData.h"

namespace sayaka 
{


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


enum BndryCondType {

	BCType_Neumann, 
	BCType_Dirichlet,
	BCType_SimpleFill, // directly set value in ghost cell
};

struct BCRegister {
	typedef std::map<int,int> bcmap_type;
	bcmap_type bcmaps[SurroundingIndex::NumSurr];

	//void fillBCMap(int phys_bc, int math_bc)

	void setBCMap(int isurr, int jsurr, int ksurr, 
		int phys_bc, int math_bc)
	{
		SurroundingIndex ind(isurr, jsurr, ksurr);
		bcmaps[ind][phys_bc] = math_bc;
	}

	int getBCMap(int isurr, int jsurr, int ksurr, 
		int phys_bc) const
	{
		// TODO
		SurroundingIndex ind(isurr, jsurr, ksurr);
		bcmap_type::const_iterator it = bcmaps[ind].find(phys_bc);
		if (it == bcmaps[ind].end()) {
			LOGPRINTF("Unknown BC mapping for phys_bc=%d\n", phys_bc);
			exit(1);
		} 
		return it->second;
	}
};

/**
 *
 */
struct BoundaryConditionPatch : public TreeDataAlgorithm
{

	typedef TreeDataAlgorithm super_type;

	IndexBox bndryBoxes[3][3][3];

	std::vector<BCRegister> bndryRegs;

public:
	BoundaryConditionPatch(
		AmrTree &tree_in, 
		TreeData &data_in, 
		VariableLocation varloc_in)
		: super_type(tree_in, data_in, varloc_in),
		bndryRegs(data_in.numComp())
	{
		calcRanges();
	}

	const BCRegister& boundaryRegister(int comp) const { 
		assert(0<=comp && comp<tree_data.numComp());
		return bndryRegs[comp];
	}
	BCRegister& boundaryRegister(int comp) { 
		assert(0<=comp && comp<tree_data.numComp());
		return bndryRegs[comp];
	}

	void fillBlockBC(int physbc, 
		int isweep, int sweepdir,
		int xbndry, int ybndry, int zbndry,
		int iblock, DoubleBlockData &data, 
		int scomp, int ncomp) const;

	//void fillCellBlockBC(int bctype, int xbndry, int ybndry, int zbndry,
	//	int iblock, DoubleBlockData &blockdata) const;

protected:
	//
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


} // namespace_sayaka

