
#define USE_BNDRY_FILL



#include "SayakaMacSolver.h"


namespace sayaka 
{

void MacSolver::prepareSolver() {
	registerMGLevels();

	// default MAC coefficients
	setDefaultCoefs();

	// non-IB volume/area fraction
	setUnitIBVolumeFrac(*m_ib_volfrac);
	setUnitIBAreaFrac(m_ib_areafrac);

	if (!useBndryFill()) {
		LOGPRINTF("MacSolver: set use_bndry_fill=1\n");
		exit(1);
	}
}
void MacSolver::clearSolver() {
	// do something...
}

//void MacSolver::setHomoBC(int phys_bc, int math_bc) {
//	BCRegister &bcreg = getBCReg();
//	for (int ksurr=-ZDIM; ksurr<=ZDIM; ksurr++) {
//	for (int jsurr=-YDIM; jsurr<=YDIM; jsurr++) {
//	for (int isurr=-XDIM; isurr<=XDIM; isurr++) {
//		bcreg.setBCMap(isurr,jsurr,ksurr, phys_bc, math_bc);
//		bcreg.setBCVal(isurr,jsurr,ksurr, phys_bc, 0.0);
//	}
//	}
//	}
//}

int MacSolver::mg_solve() {
	const AmrTree &tree = getTree();

	const int finest_ilevel = tree.currentFinestLevel();
	const int coarset_ilevel = 1;
	assert(tree.minAmrLevel<=finest_ilevel && finest_ilevel<=tree.maxAmrLevel);

	// prepare MG levels
	// NOT needed, called by the user in prepareSolver()
	// registerMGLevels();

	//
	TreeData &data = getDataBuf();

	// derive initial residual 
	mg_residual_level(finest_ilevel, data, data, data, 
		COMP_BUF_RESID, COMP_BUF_SOL, COMP_BUF_RHS);
	// compute initial residual norm
	const double resid_norm0 = mg_norm_level(finest_ilevel, data, COMP_BUF_RESID);
	if (isVerbose()) {
		LOGPRINTF("MG: |resid0| = %e\n", resid_norm0);
	}

	int max_iter = m_max_iter;
	if (mgUseFixedIter()) max_iter = m_mg_fixed_iter;
	const double tol_rel = m_tol_rel;
	const double tol_abs = m_tol_abs;
	const int nsmooth_post = m_mg_postsmooth_num;

	bool mg_conv = false;

	for (int iter=1; iter<=max_iter; iter++) {
		if (isVerbose()) {
			LOGPRINTF("MG: begin iter=%d/%d\n", iter, max_iter);
		}

		// clear correction = 0
		TreeData::SetValue(data, 0.0, COMP_BUF_CORR, 1, 0);
		
		for (int ilevel=finest_ilevel; ilevel>=coarset_ilevel+1; ilevel--) { 
			// pre-smooth
			// NOTE that the current MG does not have pre-smooth
			// residual is directly averaged down to coarse level

			// restrict residual to coarser levels
			//m_fillbuf->restrictLevel(ilevel-1, COMP_BUF_RESID, 1);
			mg_restict_wgt(ilevel-1, data, COMP_BUF_RESID);

			if (isVerbose()) {
				LOGPRINTF("Restrict RESID: level%d -> level%d\n", ilevel, ilevel-1);
			}
		}

		// solve at bottom level
		{
			int ret_bot = mg_solve_bottom(coarset_ilevel);
			if (ret_bot != 0) {
				LOGPRINTF("MG: Bottom break down\n");
				exit(1);
			}
		}

		for (int ilevel=coarset_ilevel+1; ilevel<=finest_ilevel; ilevel++) {
			// prolong correction to finer levels
			//m_fillbuf->prolongLevel(ilevel, COMP_BUF_CORR, 1);
			mg_prolong_wgt(ilevel, data, COMP_BUF_CORR);

			if (isVerbose()) {
				LOGPRINTF("Prolong CORR: level%d -> level%d\n", ilevel-1, ilevel);
			}

			// post-smooth
			for (int ismooth=1; ismooth<=nsmooth_post; ismooth++) {
				mg_smooth_level(ilevel, data, data, COMP_BUF_CORR, COMP_BUF_RESID);

				if (isVerbose()) {
					LOGPRINTF("PostSmooth: smooth=%d/%d\n", ismooth, nsmooth_post);
				}
			}
		}

		// correct on finest level
		mg_correct_level(finest_ilevel, data, data, COMP_BUF_SOL, COMP_BUF_CORR);

		if (isVerbose()) {
			LOGPRINTF("Correct: finest_level=%d\n", finest_ilevel);
		}

		// calculate current residual
		mg_residual_level(finest_ilevel, data, data, data,
			COMP_BUF_RESID, COMP_BUF_SOL, COMP_BUF_RHS);
		double resid_norm = mg_norm_level(finest_ilevel, data, COMP_BUF_RESID);
		if (isVerbose()) {
			LOGPRINTF("Residual: |resid|=%e\n", resid_norm);
		}

		if (resid_norm < resid_norm0*tol_rel + tol_abs) {
			mg_conv = true;
		} 
		if (isVerbose()) {
			LOGPRINTF("MG: end iter=%d/%d\n", iter, max_iter);
		}

		if (!mgUseFixedIter()) {
			if (mg_conv) break;
		}
	} // end MG iteration

	if (!mgUseFixedIter()) {
		if (mg_conv) {
			if (isVerbose()) {
				LOGPRINTF("MG: converge\n");
			}
		} else {
			LOGPRINTF("MG: failure\n");
		}
	}

	return 0;
}

int MacSolver::mg_solve_bottom(int bottom_level) {
	assert(bottom_level >= 1);

	int ret = 0;

	const int max_iter = m_mg_bot_max_iter;
	const double tol_abs = m_mg_bot_tol_abs;
	const double tol_rel = m_mg_bot_tol_rel;

	TreeData &data = getDataBuf();
	const int solcomp = COMP_BUF_CORR;
	const int rhscomp = COMP_BUF_RESID;
	const int tmpcomp = COMP_BUF_TMP;

	// calculate residual norm
	// NOTE residual must be ready during the error restriction
	double resid_norm0 = mg_norm_level(bottom_level, data, rhscomp);
	if (isVerbose()) {
		LOGPRINTF("Bottom: bottom_level=%d, |resid0|=%e\n", 
			bottom_level, resid_norm0);
	}

	// set initial correction to zero
	mg_zero_level(bottom_level, data, solcomp, 1);

	bool conv = false;

	if (0) { // use smoother
		int iter = 0;
		double resid_norm = 0;
		for (iter=1; iter<=max_iter; iter++) {

			mg_smooth_level(bottom_level, 
				data, data, solcomp, rhscomp);

			// calculate residual 
			mg_residual_level(bottom_level, data, data, data,
				tmpcomp, solcomp, rhscomp);
			resid_norm = mg_norm_level(bottom_level, data, tmpcomp);

			if (iter%50==0 && isVerbose()) {
				LOGPRINTF("Bottom: bot_iter=%d/%d, |resid|=%e\n", 
					iter, max_iter, resid_norm);
			}

			if (resid_norm < resid_norm0*tol_rel + tol_abs) {
				conv = true;
			}

			if (conv) break;
		}

		if (isVerbose()) {
			LOGPRINTF("Bottom: final iter=%d, |resid|=%e\n", 
				iter, resid_norm);
		}
	} else { // use CG
		int precond = 0; // do not use preconditioner

		ret = cg_solve_level(bottom_level, 
			data, data, solcomp, rhscomp, 
			tol_rel, tol_abs, max_iter,
			precond, 
			BCMODE_HOMO);

		if (ret == 0) { conv = true; }
	}

	if (conv) {
		if (isVerbose()) {
			LOGPRINTF("Bottom: converged\n");
		}
	} else {
		LOGPRINTF("Bottom: failed\n");
		ret = 1;
	}

	return ret;
}


int MacSolver::mg_relax_vcycle(int finest_level,
	TreeData &sol, const TreeData &rhs, int solcomp, int rhscomp,
	double tol_rel, double tol_abs, int max_iter, int fixed_iter,
	MacBCMode bc_mode)
{
	const AmrTree &tree = getTree();

	const int coarset_ilevel = 1;
	assert(tree.minAmrLevel<=finest_level && finest_level<=tree.maxAmrLevel);
	assert(finest_level <= tree.currentFinestLevel());

	//
	TreeData &data = getDataBuf();

	// copy initial solution to buffer
	mg_copy_level(finest_level, data, sol, COMP_BUF_SOL, solcomp, 1, 0);

	// derive initial residual 
	mg_residual_level(finest_level, data, data, rhs, 
		COMP_BUF_RESID, COMP_BUF_SOL, rhscomp);

	// compute initial residual norm
	const double rnorm0 = mg_norm_level(finest_level, data, COMP_BUF_RESID);
	if (isVerbose()) {
		LOGPRINTF("MG: |resid0| = %e\n", rnorm0);
	}

	// initialize correction = 0
	//mg_zero_level(finest_level, data, COMP_BUF_CORR, 1);

	// 
	int ret = 0;
	bool mg_conv = false;

	for (int iter=1; iter<=max_iter; iter++) {
		if (isVerbose()) {
			LOGPRINTF("MG: begin iter=%d/%d\n", iter, max_iter);
		}

		// clear correction = 0
		mg_zero_level(finest_level, data, COMP_BUF_CORR, 1);
	
		// V-cycle down
		for (int ilevel=finest_level; ilevel>=coarset_ilevel+1; ilevel--) { 
			// pre-smooth
			// NOTE that the current MG does not have pre-smooth
			// residual is directly averaged down to coarse level

			// restrict residual to coarser levels
			//m_fillbuf->restrictLevel(ilevel-1, COMP_BUF_RESID, 1);
			mg_restict_wgt(ilevel-1, data, COMP_BUF_RESID);

			if (isVerbose()) {
				LOGPRINTF("Restrict RESID: level%d -> level%d\n", ilevel, ilevel-1);
			}
		}

		// solve at bottom level
		int ret_bot = mg_solve_bottom(coarset_ilevel);
		if (ret_bot != 0) {
			LOGPRINTF("MG: Bottom break down\n");
			ret = 1; break;
		}

		// V-cycle up
		for (int ilevel=coarset_ilevel+1; ilevel<=finest_level; ilevel++) {
			// prolong correction to finer levels
			//m_fillbuf->prolongLevel(ilevel, COMP_BUF_CORR, 1);
			mg_prolong_wgt(ilevel, data, COMP_BUF_CORR);

			if (isVerbose()) {
				LOGPRINTF("Prolong CORR: level%d -> level%d\n", ilevel-1, ilevel);
			}

			// post-smooth
			const int nsmooth_post = m_mg_postsmooth_num;
			for (int ismooth=1; ismooth<=nsmooth_post; ismooth++) {
				mg_smooth_level(ilevel, data, data, COMP_BUF_CORR, COMP_BUF_RESID);

				if (isVerbose()) {
					LOGPRINTF("PostSmooth: smooth=%d/%d\n", ismooth, nsmooth_post);
				}
			}
		}

		// correct on finest level
		mg_correct_level(finest_level, data, data, COMP_BUF_SOL, COMP_BUF_CORR);

		if (isVerbose()) {
			LOGPRINTF("Correct: finest_level=%d\n", finest_level);
		}

		// calculate current residual
		mg_residual_level(finest_level, data, data, data,
			COMP_BUF_RESID, COMP_BUF_SOL, COMP_BUF_RHS);
		double rnorm = mg_norm_level(finest_level, data, COMP_BUF_RESID);
		if (isVerbose()) {
			LOGPRINTF("Residual: |resid|=%e\n", rnorm);
		}

		if (rnorm < rnorm*tol_rel + tol_abs) {
			mg_conv = true;
		} 
		if (isVerbose()) {
			LOGPRINTF("MG: end iter=%d/%d\n", iter, max_iter);
		}

		if (fixed_iter <= 0) {
			if (mg_conv) break;
		}
	} // end MG iteration

	if (fixed_iter <= 0) {
		if (mg_conv) {
			if (isVerbose()) {
				LOGPRINTF("MG: converge\n");
			}
		} else {
			ret = 1;
			LOGPRINTF("MG: failure\n");
		}
	}

	// update solution
	if (mg_conv || fixed_iter>0) {
		mg_copy_level(finest_level, sol, data, solcomp, COMP_BUF_SOL, 1, 0);
	}

	return ret;
}


int MacSolver::relax_solve()
{
	const AmrTree &tree = getTree();

	// perform relaxation on the finest level
	const int finest_ilevel = tree.currentFinestLevel();
	assert(tree.minAmrLevel<=finest_ilevel && finest_ilevel<=tree.maxAmrLevel);
	if (isVerbose()) {
		LOGPRINTF("%s: finest_level=%d\n", __FUNCTION__, finest_ilevel);
	}

	//
	TreeData &data = getDataBuf();


	// initialize correction = 0
	TreeData::SetValue(data, 0.0, COMP_BUF_CORR, 1, 0);

	// derive initial residual 
	mg_residual_level(finest_ilevel, data, data, data, COMP_BUF_RESID, COMP_BUF_SOL, COMP_BUF_RHS);
	// compute initial residual norm
	const double resid_norm0 = mg_norm_level(finest_ilevel, data, COMP_BUF_RESID);
	if (isVerbose()) {
		LOGPRINTF("%s: |resid0| = %e\n", __FUNCTION__, resid_norm0);
	}

	const int max_iter = m_max_iter;
	const double tol_rel = m_tol_rel;
	const double tol_abs = m_tol_abs;

	bool conv = false;

	for (int iter=1; iter<=max_iter; iter++) {
		if (isVerbose()) {
			LOGPRINTF("%s: iter=%d/%d\n", __FUNCTION__, iter, max_iter);
		}

		// relax solution
		if (useBndryFill()) { // use ghost
			std::vector<TreeData*> &cfsave = getFluxSaveDataBuf();
			const int savecomp = 0;

			mg_smooth_level_usefill(finest_ilevel,
				data, data, COMP_BUF_SOL, COMP_BUF_RHS,
				cfsave, savecomp);
		} else { // do not use ghost cell
			mg_smooth_level_nofill(finest_ilevel, data, data,
				COMP_BUF_SOL, COMP_BUF_RHS);
		}

		// calculate current residual
		mg_residual_level(finest_ilevel, data, data, data,
			COMP_BUF_RESID, COMP_BUF_SOL, COMP_BUF_RHS);
		// evaluate residual norm
		double resid_norm = mg_norm_level(finest_ilevel, data, COMP_BUF_RESID);
		if (isVerbose()) {
			LOGPRINTF("%s: |resid|=%e\n", __FUNCTION__, resid_norm);
		}

		if (resid_norm < resid_norm0*tol_rel + tol_abs) {
			conv = true;
		}

		if (conv) break;
	} // end iteration

	if (conv) {
		if (isVerbose()) {
			LOGPRINTF("%s: converge\n", __FUNCTION__);
		}
	} else {
		LOGPRINTF("%s: failure\n", __FUNCTION__);
	}

	return 0;
}


int MacSolver::cg_solve()
{
	const AmrTree &tree = getTree();

	// perform CG solver on the finest level
	const int finest_ilevel = tree.currentFinestLevel();
	assert(tree.minAmrLevel<=finest_ilevel && finest_ilevel<=tree.maxAmrLevel);
	if (isVerbose()) {
		LOGPRINTF("%s: finest_level=%d\n", __FUNCTION__, finest_ilevel);
	}

	//
	const MGLevelTower tower = getTower(finest_ilevel);
	assert(!tower.isEmptyLevel());

	//
	TreeData &data = getDataBuf();
	//
	const int comp_sol = COMP_BUF_SOL;
	const int comp_rhs = COMP_BUF_RHS;
	const int comp_resid = COMP_BUF_RESID;
	const int comp_corr = COMP_BUF_CORR;
	//
	const int comp_x0 = COMP_CG_X0;
	const int comp_x = COMP_CG_X;
	const int comp_r = COMP_CG_R;
	const int comp_rhat = COMP_CG_RHAT;
	const int comp_p = COMP_CG_P;
	const int comp_phat = COMP_CG_PHAT;
	const int comp_s = COMP_CG_S;
	const int comp_shat = COMP_CG_SHAT;
	const int comp_v = COMP_CG_V;
	const int comp_t = COMP_CG_T;

	//
	const int max_iter = m_max_iter;
	const double tol_rel = m_tol_rel;
	const double tol_abs = m_tol_abs;
	const int precond = m_cg_use_mgprecond;

#if (0)
	// let x0 holds the initial guess
	TreeData::Copy(data, data, comp_x0, comp_sol, 1, 0);

	// r = initial residual 
	mg_residual_level(finest_ilevel, data, data, data, 
		comp_r, comp_x0, comp_rhs);

	// rhat = r
	TreeData::Copy(data, data, comp_rhat, comp_r, 1, 0);

	// x = 0
	TreeData::SetValue(data, 0.0, comp_x, 1, 0);
	
	// initial residual norm
	const double resid_norm0 = mg_norm_level(finest_ilevel, data, comp_r);
	if (isVerbose()) {
		LOGPRINTF("%s: |resid0| = %e\n", __FUNCTION__, resid_norm0);
	}


	// BICGSTAB variables
	double rho1 = 0;
	double alpha = 0;
	double omega = 0;

	//
	bool conv = false;
	int stat = 0;
	int iter = 0;
	double rnorm = 0;

	if (precond) {
		setMGFixedIter(1);
	}

	for (iter=1; iter<=max_iter; iter++) {
		if (isVerbose()) {
			LOGPRINTF("%s: iter=%d/%d\n", __FUNCTION__, iter, max_iter);
		}

		// rho = rhat . r
		const double rho = mg_dotprod_level(finest_ilevel, data, data, comp_rhat, comp_r);
		if (rho == 0) { // this corresponds to a solver breakdown
			stat = 1; break;
		}

		if (iter == 1) {
			// p = r
			TreeData::Copy(data, data, comp_p, comp_r, 1, 0);
		} else {
			const double beta = (rho/rho1) * (alpha/omega);
			// p = p - omega*v
			mg_saxby_level(finest_ilevel, 
				data, 1.0, data, -omega, data,
				comp_p, comp_p, comp_v);
			// p = r + beta*p
			mg_saxby_level(finest_ilevel,
				data, 1.0, data, beta, data,
				comp_p, comp_r, comp_p);
		}

		// precondition M.phat = p
		if (precond == 0) {
			// phat = p
			TreeData::Copy(data, data, comp_phat, comp_p, 1, 0);
		} else { // precondition
			// phat = inv(M) . p
			setZeroSol();
			setRhs(data, comp_p);
			int verbose_save = m_verbose;
			setVerbose(0);
			mg_solve();
			getSol(data, comp_phat);
			setVerbose(verbose_save);
		}

		// v = A . phat
		mg_apply_level(finest_ilevel, data, data, comp_v, comp_phat);

		// rhat . v
		alpha = mg_dotprod_level(finest_ilevel, data, data, comp_rhat, comp_v);
		if (alpha == 0) { // another breakdown?
			stat = 1; break;
		}
		// alpha = rho / (rhat.v)
		alpha = rho / alpha;

		// s = r - alpha*v
		mg_saxby_level(finest_ilevel, 
			data, 1.0, data, -alpha, data,
			comp_s, comp_r, comp_v);
		// x = x + alpha*phat
		mg_saxby_level(finest_ilevel,
			data, 1.0, data, alpha, data,
			comp_x, comp_x, comp_phat);

		// check |s|
		rnorm = mg_norm_level(finest_ilevel, data, comp_s);
		if(isVerbose()) {
			LOGPRINTF("%s: half |resid|=%e\n", __FUNCTION__, rnorm);
		}
		if (rnorm<resid_norm0*tol_rel || rnorm<tol_abs) {
			// BICGSTAB special convergence
			conv = true; break;
		}

		// precondition M.shat = s
		if (precond == 0) {
			// shat = s
			TreeData::Copy(data, data, comp_shat, comp_s, 1, 0);
		} else {
			// shat = inv(M).s
			setZeroSol();
			setRhs(data, comp_s);
			int verbose_save = m_verbose;
			setVerbose(0);
			mg_solve();
			getSol(data, comp_shat);
			setVerbose(verbose_save);
		}

		// t = A . shat
		mg_apply_level(finest_ilevel, data, data, comp_t, comp_shat);

		// t.s & t.t
		const double t_dot_s = mg_dotprod_level(finest_ilevel, data, data, comp_t, comp_s);
		const double t_dot_t = mg_dotprod_level(finest_ilevel, data, data, comp_t, comp_t);
		if (t_dot_t == 0) { // another breakdown?
			stat = 1; break;
		} else {
			omega = t_dot_s / t_dot_t;
		}

		// x = x + omega*shat
		mg_saxby_level(finest_ilevel, 
			data, 1.0, data, omega, data,
			comp_x, comp_x, comp_shat);
		// r = s - omega*t
		mg_saxby_level(finest_ilevel,
			data, 1.0, data, -omega, data,
			comp_r, comp_s, comp_t);

		// current residual norm
		rnorm = mg_norm_level(finest_ilevel, data, comp_r);
		if (isVerbose()) {
			LOGPRINTF("%s: |resid|=%e\n", __FUNCTION__, rnorm);
		}

		// check |r|
		if (rnorm < resid_norm0*tol_rel + tol_abs) {
			// BICGSTAB converged
			conv = true; break;
		}

		// for continuation omega!=0
		if (omega == 0) {
			stat = 1; break;
		}

		// save for next iteration
		rho1 = rho;
	} // end iteration

	if (conv) {
		// add x (i.e. the correction) to initial solution
		mg_correct_level(finest_ilevel, data, data, comp_x0, comp_x);
		// copy x0 to solution buffer
		mg_copy_level(finest_ilevel, data, data, comp_sol, comp_x0, 1, 0);

		if (isVerbose()) {
			LOGPRINTF("%s: converge\n", __FUNCTION__);
		}
	} else {
		LOGPRINTF("%s: failure\n", __FUNCTION__);
	}

	return 0;
#else
	return cg_solve_level(finest_ilevel, 
		data, data, comp_sol, comp_rhs,
		tol_rel, tol_abs, max_iter,
		precond, BCMODE_INHOMO);
#endif

}


int MacSolver::cg_solve_level(int mg_level,
	TreeData &sol, const TreeData &rhs, int solcomp, int rhscomp,
	double tol_rel, double tol_abs, int max_iter, 
	int precond, MacBCMode bc_mode)
{
	int stat = 0;

	const AmrTree &tree = getTree();

	// perform CG solver on the current level
	assert(1<=mg_level && mg_level<=tree.maxAmrLevel);
	if (isVerbose()) {
		LOGPRINTF("%s: mg_level=%d\n", __FUNCTION__, mg_level);
	}

	//
	const MGLevelTower tower = getTower(mg_level);
	assert(!tower.isEmptyLevel());
	if (tower.isEmptyLevel()) {
		LOGPRINTF("%s: mg_level=%d is empty\n", __FUNCTION__, mg_level);
		exit(1);
	}

	//
	TreeData &data = getDataBuf();
	//
	const int comp_x0 = COMP_CG_X0;
	const int comp_x = COMP_CG_X;
	const int comp_r = COMP_CG_R;
	const int comp_rhat = COMP_CG_RHAT;
	const int comp_p = COMP_CG_P;
	const int comp_phat = COMP_CG_PHAT;
	const int comp_s = COMP_CG_S;
	const int comp_shat = COMP_CG_SHAT;
	const int comp_v = COMP_CG_V;
	const int comp_t = COMP_CG_T;

	// let x0 holds the initial guess
	// we need x0 because the input solution data may not have ghost cells
	mg_copy_level(mg_level, data, sol, comp_x0, solcomp, 1, 0);

	// r = initial residual 
	mg_residual_level(mg_level, data, data, rhs, 
		comp_r, comp_x0, rhscomp);

	// rhat = r
	mg_copy_level(mg_level, data, data, comp_rhat, comp_r, 1, 0);

	// x = 0
	mg_zero_level(mg_level, data, comp_x, 1);
	
	// initial residual norm
	const double rnorm0 = mg_norm_level(mg_level, data, comp_r);
	if (isVerbose()) {
		LOGPRINTF("%s: |resid0| = %e\n", __FUNCTION__, rnorm0);
	}
	if (rnorm0 == 0) {
		if (isVerbose()) {
			LOGPRINTF("%s: initial converge\n", __FUNCTION__);
		}
		return stat;
	}

	// BICGSTAB variables
	double rho1 = 0;
	double alpha = 0;
	double omega = 0;

	//
	bool conv = false;
	int iter = 0;
	double rnorm = 0;

	if (precond) { 
		// TODO
		if (mg_level == 1) { // at bottom level
			LOGPRINTF("%s: cannot use MG precondition at bottom level\n", __FUNCTION__);
			precond = 0;
		} else {
			setMGFixedIter(1);
		}
	}

	for (iter=1; iter<=max_iter; iter++) {
		if (isVerbose()) {
			LOGPRINTF("%s: iter=%d/%d\n", __FUNCTION__, iter, max_iter);
		}

		// rho = rhat . r
		const double rho = mg_dotprod_level(mg_level, data, data, comp_rhat, comp_r);
		if (rho == 0) { // this corresponds to a solver breakdown
			if (isVerbose()) LOGPRINTF("%s: rho=0\n", __FUNCTION__);
			stat = 1; break;
		}

		if (iter == 1) {
			// p = r
			mg_copy_level(mg_level, data, data, comp_p, comp_r, 1, 0);
		} else {
			const double beta = (rho/rho1) * (alpha/omega);
			// p = p - omega*v
			mg_saxby_level(mg_level, 
				data, 1.0, data, -omega, data,
				comp_p, comp_p, comp_v);
			// p = r + beta*p
			mg_saxby_level(mg_level,
				data, 1.0, data, beta, data,
				comp_p, comp_r, comp_p);
		}

		// precondition M.phat = p
		if (precond == 0) {
			// phat = p
			mg_copy_level(mg_level, data, data, comp_phat, comp_p, 1, 0);
		} else { // precondition
			mg_precond_level(mg_level, data, data, comp_phat, comp_p);
			//// phat = inv(M) . p
			//// TODO
			//setZeroSol();
			////setRhs(data, comp_p);
			//mg_copy_level(mg_level, data, data, COMP_BUF_RHS, comp_p, 1, 0);
			//int verbose_save = m_verbose;
			////setVerbose(0);
			//mg_solve();
			////getSol(data, comp_phat);
			//mg_copy_level(mg_level, data, data, comp_phat, COMP_BUF_SOL, 1, 0);
			//setVerbose(verbose_save);
		}

		// v = A . phat
		mg_apply_level(mg_level, data, data, comp_v, comp_phat);

		// rhat . v
		alpha = mg_dotprod_level(mg_level, data, data, comp_rhat, comp_v);
		if (alpha == 0) { // another breakdown?
			if (isVerbose()) LOGPRINTF("%s: alpha=0\n", __FUNCTION__);
			stat = 1; break;
		}
		// alpha = rho / (rhat.v)
		alpha = rho / alpha;

		// s = r - alpha*v
		mg_saxby_level(mg_level, 
			data, 1.0, data, -alpha, data,
			comp_s, comp_r, comp_v);
		// x = x + alpha*phat
		mg_saxby_level(mg_level,
			data, 1.0, data, alpha, data,
			comp_x, comp_x, comp_phat);

		// check |s|
		rnorm = mg_norm_level(mg_level, data, comp_s);
		if(isVerbose() && Verbose()>=2) {
			LOGPRINTF("%s: half |resid|=%e\n", __FUNCTION__, rnorm);
		}
		if (rnorm<rnorm0*tol_rel || rnorm<tol_abs) {
			// BICGSTAB special convergence
			conv = true; break;
		}

		// precondition M.shat = s
		if (precond == 0) {
			// shat = s
			mg_copy_level(mg_level, data, data, comp_shat, comp_s, 1, 0);
		} else {
			mg_precond_level(mg_level, data, data, comp_shat, comp_s);
			//int verbose_save = m_verbose;
			////setVerbose(0);
			//// shat = inv(M).s
			//// TODO
			//setZeroSol();
			////setRhs(data, comp_s);
			//mg_copy_level(mg_level, data, data, COMP_BUF_RHS, comp_s, 1, 0);
			//
			//mg_solve();
			////getSol(data, comp_shat);
			//mg_copy_level(mg_level, data, data, comp_shat, COMP_BUF_SOL, 1, 0);
			//setVerbose(verbose_save);
		}

		// t = A . shat
		mg_apply_level(mg_level, data, data, comp_t, comp_shat);

		// t.s & t.t
		const double t_dot_s = mg_dotprod_level(mg_level, data, data, comp_t, comp_s);
		const double t_dot_t = mg_dotprod_level(mg_level, data, data, comp_t, comp_t);
		if (t_dot_t == 0) { // another breakdown?
			stat = 1; break;
		} else {
			omega = t_dot_s / t_dot_t;
		}

		// x = x + omega*shat
		mg_saxby_level(mg_level, 
			data, 1.0, data, omega, data,
			comp_x, comp_x, comp_shat);
		// r = s - omega*t
		mg_saxby_level(mg_level,
			data, 1.0, data, -omega, data,
			comp_r, comp_s, comp_t);

		// current residual norm
		rnorm = mg_norm_level(mg_level, data, comp_r);
		if (isVerbose()) {
			LOGPRINTF("%s: |resid|=%e\n", __FUNCTION__, rnorm);
		}

		// check |r|
		if (rnorm < rnorm0*tol_rel + tol_abs) {
			// BICGSTAB converged
			conv = true; break;
		}

		// for continuation omega!=0
		if (omega == 0) {
			if (isVerbose()) LOGPRINTF("%s: omega=0\n", __FUNCTION__);
			stat = 1; break;
		}

		// save for next iteration
		rho1 = rho;
	} // end iteration

	if (conv) {
		// add x (i.e. the correction) to initial solution
		mg_correct_level(mg_level, data, data, comp_x0, comp_x);
		// copy x0 to solution buffer
		mg_copy_level(mg_level, sol, data, solcomp, comp_x0, 1, 0);

		if (isVerbose()) {
			LOGPRINTF("%s: converge\n", __FUNCTION__);
		}
	} else {
		LOGPRINTF("%s: failure\n", __FUNCTION__);
	}

	return stat;
}


// sol = inv(M) . rhs
int MacSolver::mg_precond_level(int mg_level,
	TreeData &sol, const TreeData &rhs, int solcomp, int rhscomp)
{
	int verbose_save = m_verbose;
	setVerbose(0);

	// use zero guess
	mg_zero_level(mg_level, sol, solcomp, 1);
	
	// dummy 
	double tol_rel = m_tol_rel;
	double tol_abs = m_tol_abs;
	// only one pass
	int max_iter = 1;
	int fixed_iter = 1;

	//
	int ret = mg_relax_vcycle(mg_level, 
		sol, rhs, solcomp, rhscomp,
		tol_rel, tol_abs, max_iter, fixed_iter,
		BCMODE_HOMO);

	setVerbose(verbose_save);

	return ret;
}


} // namespace_sayaka



