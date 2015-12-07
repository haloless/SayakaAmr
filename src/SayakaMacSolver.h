#pragma once


#include "SayakaCommons.h"
#include "SayakaTree.h"
#include "SayakaTreeData.h"
#include "SayakaFillPatch.h"
#include "SayakaTower.h"

namespace sayaka
{




/**
 * Solve MAC projection-like operator
 * L(phi) = alpha.A.phi - beta.div(B.grad(phi)) = rhs
 */
class MacSolver
{
public:

	enum MacBCMode {
		BCMODE_HOMO,
		BCMODE_INHOMO,
	};

	enum SolverBufComp {
		COMP_BUF_RHS,
		COMP_BUF_SOL,
		COMP_BUF_RESID,
		COMP_BUF_CORR,
		COMP_BUF_TMP,

		// 
		COMP_CG_X0,
		COMP_CG_X,
		COMP_CG_R, COMP_CG_RHAT,
		COMP_CG_P, COMP_CG_PHAT,
		COMP_CG_S, COMP_CG_SHAT,
		COMP_CG_V,
		COMP_CG_T,

		NUM_COMP_BUF,

		NUM_GROW_BUF = 1,
	};

	MacSolver(const AmrTree &tree)
		: m_tree(&tree), 
		m_buf(NULL), m_fillbuf(NULL),
		m_ib_volfrac(NULL), m_ib_areafrac(),
		m_acoef(NULL), m_bcoef(),
		m_verbose(0)
	{
		// generate data structure
		define(tree);

		// define parameters
		// coefficients
		m_alpha = 0.0;
		m_beta = 1.0;

		// solver
		m_max_iter = 2000;
		m_tol_rel = 1.0e-6;
		m_tol_abs = 0.0;
		// MG
		m_mg_presmooth_num = 2;
		m_mg_postsmooth_num = 2;
		m_mg_bot_tol_abs = 1.0e-6;
		m_mg_bot_tol_rel = 1.0e-6;
		m_mg_bot_max_iter = 2000;
		m_mg_fixed_iter = 0;
		// CG
		m_cg_use_mgprecond = 0;

		// 
		m_use_bndryfill = 1;
	}

	~MacSolver() {
		if (m_buf) {
			delete m_buf; m_buf = NULL;
		}
		if (m_fillbuf) {
			delete m_fillbuf; m_fillbuf = NULL;
		}
		for (int dir=0; dir<NDIM; dir++) {
			if (m_crse_fine_save[dir]) {
				delete m_crse_fine_save[dir];
				m_crse_fine_save[dir] = NULL;
			}
		}

		if (m_ib_volfrac) {
			delete m_ib_volfrac; m_ib_volfrac = NULL;
		}
		TreeData::ReleaseDataPArray(m_ib_areafrac);

		if (m_acoef) {
			delete m_acoef; m_acoef = NULL;
		}
		TreeData::ReleaseDataPArray(m_bcoef);
	}

	void define(const AmrTree &tree);

	//
	void setRhs(const TreeData &rhsdata, int rhscomp) {
		TreeData &buf = getDataBuf();
		TreeData::Copy(buf, rhsdata, COMP_BUF_RHS, rhscomp, 1, 0);
	}
	void setInitSol(const TreeData &soldata, int solcomp) {
		TreeData &buf = getDataBuf();
		TreeData::Copy(buf, soldata, COMP_BUF_SOL, solcomp, 1, 0);
	}
	void setZeroSol() {
		TreeData &buf = getDataBuf();
		TreeData::SetValue(buf, 0.0, COMP_BUF_SOL, 1, 0);
	}
	void getSol(TreeData &soldata, int solcomp) {
		const TreeData &buf = getDataBuf();
		TreeData::Copy(soldata, buf, solcomp, COMP_BUF_SOL, 1, 0);
	}

	//
	void setUnitIBFrac();
	//
	void setUnitIBVolumeFrac(TreeData &volfrac, int volcomp=0) const;
	void setUnitIBAreaFrac(std::vector<TreeData*> &areafrac, int areacomp=0) const;
	void setUnitIBAreaFrac(int dir, TreeData &areafrac, int areacomp=0) const;
	//
	void setIBFracByLS(const TreeData &ls, int lscomp);
	void calcIBFracByLS(const TreeData &ls, 
		TreeData &volfrac, std::vector<TreeData*> &areafrac,
		int lscomp, int volcomp, int areacomp) const;
	void calcBlockIBVolumeFrac(int iblock,
		const TreeData &ls, TreeData &volfrac,
		int lscomp, int volcomp) const;
	void calcBlockIBAreaFrac(int iblock, int dir,
		const TreeData &ls, TreeData &areafrac,
		int lscomp, int areacomp) const;

	// alpha=0, beta=1, A=0, B[]=1
	void setDefaultCoefs();
	void setUniformCoefA(double a);
	void setUniformCoefB(double b);
	//
	void setCoefB(const std::vector<TreeData*> &bdata, int bcomp);
	void setCoefB(int dir, const TreeData &bdata, int bcomp);

	// must call each time before use!!!
	void prepareSolver();
	// must call each time after use!!!
	void clearSolver();

	void registerMGLevels();
	void registerMGLevel(int mg_level);

	//const MGLevelRegister& getMGLevelReg(int level) const { return m_level_reg[level]; }
	//MGLevelRegister& getMGLevelReg(int level) { return m_level_reg[level]; }
	// deprecated!!
	const MGLevelTower& getMGLevelReg(int level) const { return m_tower[level]; }
	MGLevelTower& getMGLevelReg(int level) { return m_tower[level]; }

	const MGLevelTower& getTower(int level) const { return m_tower[level]; }
	MGLevelTower& getTower(int level) { return m_tower[level]; }

	/*
	 *
	 */


	/*
	 * solver routines
	 */
	// MultiGrid on AMR structure
	int mg_solve();
	//
	int mg_solve_bottom(int bottom_level);

	// Relaxation on Leaf level
	int relax_solve();

	// BICGSTAB on Leaf level
	int cg_solve();



	int cg_solve_level(int mg_level,
		TreeData &sol, const TreeData &rhs, int solcomp, int rhscomp,
		double tol_rel, double tol_abs, int max_iter, 
		int precond,
		MacBCMode bc_mode);

	int mg_relax_vcycle(int finest_level,
		TreeData &sol, const TreeData &rhs, int solcomp, int rhscomp,
		double tol_rel, double tol_abs, int max_iter, int fixed_iter,
		MacBCMode bc_mode);

	int mg_precond_level(int mg_level,
		TreeData &sol, const TreeData &rhs, int solcomp, int rhscomp);


	/*
	 * Implementation routines
	 * Those 'use_fill' routines assumes ghost cells have been properly filled!!!
	 */
	// Smooth operator.
	// If boundary fill is used, the solver internal flux save data will be used
	// to hold the fine-crse contribution during a down-sweep from fine to coarse tree levels.
	void mg_smooth_level(int mg_level,
		TreeData &sol, const TreeData &rhs,
		int solcomp, int rhscomp,
		int redblack_order=0);
	// smooth, no bndry fill
	void mg_smooth_level_nofill(int mg_level, 
		TreeData &phi, const TreeData &rhs,
		int phicomp, int rhscomp);
	void mg_smooth_block_nofill(int mg_level, int iblock, 
		TreeData &phinew, const TreeData &rhs,
		int phicomp, int rhscomp);
	// smooth, use bndry fill
	void mg_smooth_level_usefill(int mg_level,
		TreeData &phi, const TreeData &rhs, int dstcomp, int rhscomp,
		std::vector<TreeData*> &fluxsave, int savecomp,
		int redblack_order=0);
	void mg_smooth_block_rbcycle_usefill(int mg_level, int iblock,
		TreeData &phi, const TreeData &rhs, int dstcomp, int rhscomp,
		std::vector<TreeData*> &fluxsave, int savecomp,
		int rb_phase);


	// Apply the MAC operator on mg_level.
	// out = L(phi)
	// When called on the finest level, can be used to estimate residual.
	// Currently B coefficient is not implemented.
	void mg_apply_level(int mg_level, 
		TreeData &out, TreeData &in,
		int outcomp, int incomp);
	// block-based apply, no bndry fill
	void mg_apply_block_nofill(int mg_level, int iblock,
		TreeData &outdata, const TreeData &indata,
		int dstcomp, int srccomp);
	// block-based apply, use bndry fill
	void mg_apply_block_usefill(int mg_level, int iblock,
		TreeData &outdata, const TreeData &indata, int dstcomp, int srccomp,
		std::vector<TreeData*> &savedata, int savecomp);

	// Derive the flux = -B.grad(phi) on the mg_level.
	// The result is of course face-centered flux.
	// When called on the finest level, can be used to correct divergence-free MAC velocity.
	// Currently B coefficient is not implemented.
	void mg_flux_level(int mg_level,
		int dir, TreeData &flux,
		TreeData &phi,
		int fluxcomp, int phicomp);
	// block-based flux, no bndry fill
	void mg_flux_block_nofill(int mg_level, int iblock,
		int dir, TreeData &flux, const TreeData &phi,
		int fluxcomp, int phicomp);
	void mg_flux_block_nofill2(int mg_level, int iblock,
		int dir, TreeData &flux, const TreeData &phi,
		int fluxcomp, int phicomp);
	// block-based flux, use bndry fill
	void mg_flux_block_usefill(int mg_level, int iblock,
		int dir, TreeData &flux, const TreeData &phi,
		int fluxcomp, int phicomp);

	/*
	 * 
	 */
	// Calculate residual = RHS - L(phi) on mg_level.
	void mg_residual_level(int mg_level,
		TreeData &resid, TreeData &phi, const TreeData &rhs,
		int dstcomp, int srccomp, int rhscomp);

	// Evaluate L(inf) norm of vector on mg_level.
	// Usually this is applied to the residual on finest level.
	double mg_norm_level(int mg_level,
		const TreeData &resid, int comp);
	// result = a . b
	double mg_dotprod_level(int mg_level,
		const TreeData &a, const TreeData &b,
		int acomp, int bcomp);
	// Add correction to solution on mg_level.
	// solution += correction
	void mg_correct_level(int mg_level, 
		TreeData &sol, const TreeData &corr,
		int solcomp, int corrcomp);

	// Utility funcition to zero blocks on mg_level.
	// Can be used to reset residual/correction/flux before calculation.
	void mg_zero_level(int mg_level, TreeData &data, int scomp, int ncomp);
	void mg_zero_tree_level(int mg_level, int ilevel, TreeData &data, int scomp, int ncomp);

	void mg_copy_level(int mg_level, TreeData &dst, TreeData &src, 
		int dcomp, int scomp, int ncomp, int ngrow=0);
	void mg_saxby_level(int mg_level, TreeData &out,
		double a, const TreeData &x, double b, const TreeData &y,
		int outcomp, int xcomp, int ycomp);


	// Fill ghost cells on mg_level.
	// This is done from coarset level.
	// Piecewise-Constant interpolation is used for fine/crse boundary.
	// When this function returns, fine/fine and fine/crse boundaries are ready;
	// however, crse/fine boundary is not properly treated.
	// It is expected that crse/fine boundary will be automatically treated 
	// as we process tree data in a sweep from the finest to the coarest levels.
	void mg_fillbndry_level(int mg_level, TreeData &inout, 
		int scomp, int ncomp, int nlayer=1);

	//
	void mg_restict_wgt(int crse_ilevel, TreeData &phi, int comp);
	void mg_prolong_wgt(int fine_ilevel, TreeData &phi, int comp);


public:
	//
	const double& getAlpha() const { return m_alpha; }
	const double& getBeta() const { return m_beta; }
	void setAlpha(double alpha) { m_alpha = alpha; }
	void setBeta(double beta) { m_beta = beta; }
	void setScalars(double alpha, double beta) { m_alpha = alpha; m_beta = beta; } 

	//
	void setMaxIteration(int max_iter) { m_max_iter = max_iter; }
	void setTolerance(double tol_rel, double tol_abs) { m_tol_rel = tol_rel; m_tol_abs = tol_abs; }
	
	void setPreSmoothNum(int n_smooth) { m_mg_presmooth_num = n_smooth; }
	void setPostSmoothNum(int n_smooth) { m_mg_postsmooth_num = n_smooth; }

	void setMGFixedIter(int fixed_iter) { m_mg_fixed_iter = fixed_iter; }
	bool mgUseFixedIter() const { return m_mg_fixed_iter > 0; }

	void setCGUseMGPrecond(int use_mgprecond) { m_cg_use_mgprecond = use_mgprecond; }

	//
	void setUseBndryFill(int use_bndryfill) { m_use_bndryfill = use_bndryfill; }
	int useBndryFill() const { return m_use_bndryfill; }

	//
	const AmrTree& getTree() const { return *m_tree; }
	TreeData& getDataBuf() { return *m_buf; }
	std::vector<TreeData*>& getFluxSaveDataBuf() { return m_crse_fine_save; }
	TreeData& getFluxSaveDataBuf(int dir) { return *m_crse_fine_save[dir]; }

	void setVerbose(int verbose) { m_verbose = verbose; }
	int isVerbose() const { return m_verbose; }
	int Verbose() const { return m_verbose; }
	
protected:
	const AmrTree *m_tree;

	// solver local buffer
	TreeData *m_buf;
	FillPatch *m_fillbuf;
	std::vector<TreeData*> m_crse_fine_save;

	// volume/area IB fraction
	TreeData *m_ib_volfrac;
	std::vector<TreeData*> m_ib_areafrac;

	MGLevelTowerArray m_tower;

	// coefficients
	double m_alpha, m_beta;
	TreeData *m_acoef;
	std::vector<TreeData*> m_bcoef;

	// solver parameters
	int m_max_iter;
	double m_tol_abs, m_tol_rel;
	//
	// NOTE the current MG algorithm does not have pre-smooth
	int m_mg_presmooth_num;
	int m_mg_postsmooth_num;
	double m_mg_bot_tol_abs;
	double m_mg_bot_tol_rel;
	int m_mg_bot_max_iter;
	int m_mg_fixed_iter; // useful for MG at preconditioner

	//
	int m_cg_use_mgprecond;

	// if true, use ghost cell data to compute apply/smooth/flux
	// if false, use direct index of neighborhood cells
	// default set to 1
	int m_use_bndryfill;

	//
	int m_verbose;

private:
	// disable copy
	MacSolver(const MacSolver&);
	MacSolver& operator=(const MacSolver&);
}; // class_sayakamacsolver



} // namespace_sayaka
