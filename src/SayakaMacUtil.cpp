

#include "SayakaMacSolver.h"


namespace sayaka 
{

void MacSolver::mg_residual_level(int mg_level,
	TreeData &resid, TreeData &phi, const TreeData &rhsdata,
	int dstcomp, int srccomp, int rhscomp)
{
	// first apply MAC operator
	// save result in the residual data
	// i.e. resid = L(phi)
	mg_apply_level(mg_level, resid, phi, dstcomp, srccomp);

	// then subtract L(phi) from rhs
	// i.e. resid = rhs - resid = rhs - L(phi)

	const AmrTree &tree = getTree();
	const IndexBox &validbox = tree.validBlockCellBox();
	decl_box_range(validbox, i,j,k);

	// the MG level
	const MGLevelTower &mgLevelReg = getTower(mg_level);

	for (int igrid=0; igrid<mgLevelReg.numLevelBlock; igrid++) {
		const int iblock = mgLevelReg[igrid];
		assert(tree[iblock].getLevel() <= mg_level);

		DoubleBlockData &res = resid[iblock];
		const DoubleBlockData &rhs = rhsdata[iblock];
		const DoubleBlockData &vfrac = (*m_ib_volfrac)[iblock];

		for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
		for (int i=ilo; i<=ihi; i++) {
			double Lphi = res(i,j,k,dstcomp);
			if (vfrac(i,j,k,0) > 0) {
				res(i,j,k,dstcomp) = rhs(i,j,k,rhscomp) - Lphi;
			} else {
				res(i,j,k,dstcomp) = 0;
			}
		}
		}
		}
	}
}


// return Inf normal
double MacSolver::mg_norm_level(int mg_level,
	const TreeData &resid, int comp)
{
	const AmrTree &tree = getTree();
	const IndexBox &validbox = tree.validBlockCellBox();

	// the MG level
	const MGLevelTower &mgLevelReg = getTower(mg_level);
	assert(!mgLevelReg.isEmptyLevel());

	double norm_inf = 0;

	for (int igrid=0; igrid<mgLevelReg.numLevelBlock; igrid++) {
		const int iblock = mgLevelReg[igrid];
		assert(tree[iblock].getLevel() <= mg_level);

		const DoubleBlockData &r = resid[iblock];
		double rnorm = DoubleGridDataUtil.ReduceMaxAbs(r, validbox, comp);

		norm_inf = std::max(norm_inf, rnorm);
	}

	return norm_inf;
}

//
double MacSolver::mg_dotprod_level(int mg_level,
	const TreeData &adata, const TreeData &bdata, int acomp, int bcomp)
{
	const AmrTree &tree = getTree();
	const IndexBox &validbox = tree.validBlockCellBox();

	// the MG level
	const MGLevelTower &mgLevelReg = getTower(mg_level);
	assert(!mgLevelReg.isEmptyLevel());

	double dot = 0;

	for (int igrid=0; igrid<mgLevelReg.numLevelBlock; igrid++) {
		const int iblock = mgLevelReg[igrid];
		assert(tree[iblock].getLevel() <= mg_level);

		const DoubleBlockData &a = adata[iblock];
		const DoubleBlockData &b = bdata[iblock];

		double dot_block = DoubleGridDataUtil.DotProd(
			adata[iblock], bdata[iblock], acomp, bcomp,
			validbox);

		dot += dot_block;
	}

	return dot;
}

// sol += corr
void MacSolver::mg_correct_level(int mg_level, 
	TreeData &sol, const TreeData &corr,
	int solcomp, int corrcomp) 
{
	const AmrTree &tree = getTree();
	
	const IndexBox &validbox = tree.validBlockCellBox();
	//decl_box_range(validbox, i,j,k);

	const MGLevelTower &mgLevelReg = getTower(mg_level);
	assert(!mgLevelReg.isEmptyLevel());

	for (int igrid=0; igrid<mgLevelReg.numLevelBlock; igrid++) {
		const int iblock = mgLevelReg[igrid];

		DoubleBlockData &solbuf = sol[iblock];
		const DoubleBlockData &corrbuf = corr[iblock];

		// solution += correction
		DoubleGridDataUtil.AddEqual(
			solbuf, corrbuf, 
			solcomp, corrcomp, 1, 
			validbox);
	}
}

//
void MacSolver::mg_zero_level(int mg_level, TreeData &data, int scomp, int ncomp) {
	const AmrTree &tree = getTree();
	if (tree.getTreeLevel(mg_level).isEmptyLevel()) return;

	const MGLevelTower &mgLevelReg = getTower(mg_level);
	assert(!mgLevelReg.isEmptyLevel());

	for (int igrid=0; igrid<mgLevelReg.numLevelBlock; igrid++) {
		const int iblock = mgLevelReg[igrid];
		assert(iblock >= 0);
		assert(tree[iblock].getLevel() <= mg_level);
		
		data[iblock].setValue(scomp, ncomp, 0.0);
	}
}
//
void MacSolver::mg_zero_tree_level(int mg_level, int ilevel, 
	TreeData &data, int scomp, int ncomp) 
{
	const AmrTree &tree = getTree();
	if (tree.getTreeLevel(mg_level).isEmptyLevel()) return;

	const MGLevelTower &mgTower = getTower(mg_level);
	assert(!mgTower.isEmptyLevel());

	assert(1<=ilevel && ilevel<=mg_level);
	const int level_begin = mgTower.begin(ilevel);
	const int level_end = mgTower.end(ilevel);

	for (int igrid=level_begin; igrid<level_end; igrid++) {
		const int iblock = mgTower[igrid];
		assert(iblock >= 0);
		assert(tree[iblock].getLevel() == ilevel);
		
		data[iblock].setValue(scomp, ncomp, 0.0);
	}
}

//
void MacSolver::mg_copy_level(int mg_level, 
	TreeData &dst, TreeData &src, 
	int dcomp, int scomp, int ncomp, int ngrow) 
{
	//const AmrTree &tree = getTree();
	const MGLevelTower &tower = getTower(mg_level);
	if (tower.isEmptyLevel()) return;

	TreeData::Copy(dst, src, dcomp, scomp, ncomp, ngrow,
		tower.levelBlock, 0, tower.numLevelBlock);
}

void MacSolver::mg_saxby_level(int mg_level, TreeData &outdata,
	double a, const TreeData &xdata, double b, const TreeData &ydata,
	int outcomp, int xcomp, int ycomp)
{
	const AmrTree &tree = getTree();
	const MGLevelTower &tower = getTower(mg_level);
	if (tower.isEmptyLevel()) return;

	const IndexBox &validbox = tree.validBlockCellBox();

	for (int igrid=0; igrid<tower.numLevelBlock; igrid++) {
		const int iblock = tower[igrid];

		DoubleGridDataUtil.Calc_axpby(
			outdata[iblock],
			a, xdata[iblock],
			b, ydata[iblock],
			validbox, validbox, validbox,
			outcomp, xcomp, ycomp, 1);
	}
}


void MacSolver::mg_fillbndry_level(int mg_level, 
	TreeData &inout, int scomp, int ncomp, int nlayer) 
{
	if (nlayer <= 0) return;
	if (nlayer != 1) {
		LOGPRINTF("%s: nlayer=%d must be 1!\n", __FUNCTION__, nlayer);
		exit(1);
	}

	assert(inout.isCellData());
	assert(0<=scomp && scomp+ncomp<=inout.numComp());
	assert(nlayer <= inout.numGrow());

	const AmrTree &tree = getTree();
	if (tree.getTreeLevel(mg_level).isEmptyLevel()) return;

	const MGLevelTower &mgLevelReg = getTower(mg_level);
	assert(!mgLevelReg.isEmptyLevel());

	// 
	FillPatch fill(const_cast<AmrTree&>(tree), inout);
	// must use direct injection for fine-coarse boundary fill
	for (int comp=scomp; comp<scomp+ncomp; comp++) {
		fill.setProlongation(comp, InterpPatch::PROLONG_INJECTION);
	}

	// fill from coarse level
	for (int ilevel=1; ilevel<=mg_level; ilevel++) {
		
		if (ilevel<mg_level && mgLevelReg.size(ilevel+1)<=0) {
			assert(mgLevelReg.size(ilevel) == 0);
			// skip if next level is empty
			continue;
		}

		const int onlySameLevel = 0; // must fill from coarse level
		const int fillPhysBC = 0; // do not need BC
		const int fillCorner = 0; // do not need cross corner
		fill.fillLevelBoundary(ilevel, scomp, ncomp, nlayer,
			onlySameLevel, fillPhysBC, fillCorner);
	}
}




} // namespace sayaka



