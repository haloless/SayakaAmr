#pragma once

#include <cstdlib>
#include <cassert>

#include <vector>
#include <map>
#include <algorithm>

#include "vector3d.h"
#include "log.h"

#include "SayakaCommons.h"
#include "SayakaBox.h"
#include "SayakaTree.h"
#include "SayakaTreeData.h"
// TODO remove this dependency
#include "SayakaBoundaryPatch.h"


namespace sayaka 
{

/**
 * Interpolation for blocks from different levels
 */
struct InterpPatch : public TreeDataAlgorithm
{
	/**
	 * Coarse -> Fine
	 */
	enum ProlongMethod {
		PROLONG_IGNORE = -1,
		PROLONG_INJECTION = 0,
		PROLONG_CENTERED_LINEAR, 
		PROLONG_CENTERED_LINEAR_LIMITED,
	};

	/**
	 * Fine -> Coarse
	 */
	enum RestrictMethod {
		RESTRICT_IGNORE = -1,
		RESTRICT_AVERAGE = 0,
		RESTRICT_SUM,
	};

	typedef TreeDataAlgorithm super_type;

	std::vector<int> prolong_flag;
	std::vector<int> restrict_flag;

	//
	//const AmrTree *m_ptree;

public:

	InterpPatch(AmrTree &tree_in, TreeData &data_in, VariableLocation varloc_in)
		: super_type(tree_in, data_in, varloc_in),
		prolong_flag(data_in.numComp(), PROLONG_INJECTION),
		restrict_flag(data_in.numComp(), RESTRICT_AVERAGE)
	{
		assert(varloc_in == data_in.validBox().type());
	}

	void setProlongation(int comp, int prolong) {
		assert(0<=comp && comp<tree_data.numComp());
		prolong_flag[comp] = prolong;
	}
	void setRestriction(int comp, int restrict) {
		assert(0<=comp && comp<tree_data.numComp());
		restrict_flag[comp] = restrict;
	}
	const std::vector<int>& prolongFlag() const { return prolong_flag; }
	const std::vector<int>& restrictFlag() const { return restrict_flag; }


	//
public:
	void prolongCrseToFineBox(
		int ifine, DoubleBlockData &finedata,
		int icrse, const DoubleBlockData &crsedata,
		const IndexBox &finebox, 
		const Vector3i &fineoffset,
		int scomp, int ncomp) const;

	void restrictFineToCrseBox(
		int icrse, DoubleBlockData &crsedata,
		int ifine, const DoubleBlockData &finedata,
		const IndexBox &crsebox,
		const Vector3i &fineoffset,
		int scomp, int ncomp) const;

protected:
	/*
	 * cell-centered prolongation
	 */
	//
	void prolongCellBlockData_Injection(
		int ifine, DoubleBlockData &finedata,
		int icrse, const DoubleBlockData &crsedata,
		const IndexBox &finebox, const Vector3i &fineoffset,
		int comp) const;
	//
	void prolongCellBlockData_Bilinear(
		int ifine, DoubleBlockData &finedata,
		int icrse, const DoubleBlockData &crsedata,
		const IndexBox &finebox, const Vector3i &fineoffset,
		int comp) const;
	//
	void prolongCellBlockData_CenteredLinear(
		int ifine, DoubleBlockData &finedata,
		int icrse, const DoubleBlockData &crsedata,
		const IndexBox &finebox, const Vector3i &fineoffset,
		int comp) const;

	//
	void prolongCellBlockData_CenteredLinearLimited(
		int ifine, DoubleBlockData &finedata,
		int icrse, const DoubleBlockData &crsedata,
		const IndexBox &finebox, const Vector3i &fineoffset,
		int comp) const;

	/*
	 * cell-centered restriction
	 */
	void restrictCellBlockData_Average(
		int icrse, DoubleBlockData &crsedata,
		int ifine, const DoubleBlockData &finedata,
		const IndexBox &crsebox,
		const Vector3i &fineoffset,
		int comp, int do_average) const;

	/*
	 * face-centered
	 */
	//
	void prolongFaceBlockData_Injection(int dir,
		int ifine, DoubleBlockData &finedata,
		int icrse, const DoubleBlockData &crsedata,
		const IndexBox &finebox, const Vector3i &fineoffset,
		int comp) const;
	//
	void prolongFaceBlockData_CenteredLinear(int dir,
		int ifine, DoubleBlockData &finedata,
		int icrse, const DoubleBlockData &crsedata,
		const IndexBox &finebox, const Vector3i &fineoffset,
		int comp,
		int use_limiter) const;

	//
	void restrictFaceBlockData_Average(int dir,
		int icrse, DoubleBlockData &crsedata,
		int ifine, const DoubleBlockData &finedata,
		const IndexBox &crsebox, 
		const Vector3i &fineoffset,
		int comp, int do_average) const;


private:
	inline double centered_slope(double ul, double uc, double ur) const {
		double gradc = 0.5 * (ur - ul);
		return gradc;
	}
	inline double limited_slope(double ul, double uc, double ur) const {
		double gradc = 0.5 * (ur - ul);
		double gradl = uc - ul;
		double gradr = ur - uc;

		double signc = 0;
		if (gradc > 0) signc = 1;
		else if (gradc < 0) signc = -1;

		// MC centered limiter
		// min{abs(gc),2*s*gl,2*s*gr}
		double grad = abs(gradc);
		grad = std::min(grad, 2.0*signc*gradl);
		grad = std::min(grad, 2.0*signc*gradr);
		// max{0.0, _}
		grad = std::max(grad, 0.0);
		grad *= signc;

		return grad;
	}
}; // class_interppatch




} // namespace_sayaka

