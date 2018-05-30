#pragma once


#include <cmath>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>

#include "log.h"
#include "vector3d.h"

#include "SayakaCommons.h"
#include "Sayaka.h"
#include "SayakaTree.h"
#include "SayakaTreeData.h"
#include "SayakaBoundaryPatch.h"
#include "SayakaFillPatch.h"
#include "SayakaVarLocInterpolator.h"
#include "SayakaMacSolver.h"
#include "SayakaWriter.h"

using namespace sayaka;


class Test7 
{
protected:

	//
	enum {
		Exterior = -100,
	};

	int m_visc_integ; // 

	Vector3d m_grav;
	
	double m_dt;
	double m_nu;

	//
	std::unique_ptr<AmrTree> m_tree;

	int m_min_level;
	int m_max_level;
	int m_finest_level;
	
	int m_ng;

	IndexBox m_valid_box;
	IndexBox m_grown_box;

	// variables
	TreeData m_pres; // cell pressure
	TreeData m_vel[NDIM]; // face velocity
	TreeData m_vel_tmp[NDIM];
	TreeData m_gradp[NDIM];

	TreeData m_rhs; // 


	// BC
	std::unique_ptr<BCPatch> m_bc_pres;

	BCPatch m_bc_vel[NDIM];

	// fill
	std::unique_ptr<FillPatch> m_fill_pres;
	std::vector<std::unique_ptr<FillPatch> > m_fill_vel;

	// solver
	std::unique_ptr<MacSolver> m_solver;

	MGLevelTower *m_tower;

public:

	void init(int argc, char *argv[]) {

		m_visc_integ = 1;

		m_grav = { 1, 0, 0 };
		m_dt = 1.0e-3;
		m_nu = 0.05;

		//
		m_min_level = 2;
		m_max_level = 3;

		const int nbx = 4;
		const int nby = 4;
		const int nbz = 4;

		//
		for (int i = 1; i < argc; i += 2) {
			if (strcmp(argv[i], "-viscinteg") == 0) {
				m_visc_integ = atoi(argv[i + 1]);
			}
			else if (strcmp(argv[i], "-dt") == 0) {
				m_dt = atof(argv[i + 1]);
			}
			else if (strcmp(argv[i], "-nu") == 0) {
				m_nu = atof(argv[i + 1]);
			}
			else if (strcmp(argv[i], "-minlvl") == 0) {
				m_min_level = atoi(argv[i + 1]);
			}
			else if (strcmp(argv[i], "-maxlvl") == 0) {
				m_max_level = atoi(argv[i + 1]);
			}
			else {
				std::cerr << "Invalid argument" << std::endl;
				std::abort();
			}
		}


		//
		m_ng = 2;
		m_valid_box = IndexBox({ 0, 0, 0 }, { nbx - 1, nby - 1, nbz - 1 });
		m_grown_box = IndexBox(m_valid_box).extend(m_ng);

		//
		this->init_tree();
		this->refine_tree();

		this->init_data();

		this->init_solver();

		//
		{
			//m_vel[0].setValue(0.1);
			fill_vel_all();
		}

		for (int step = 1; step <= 100; step++) {

			if (m_visc_integ == 1) {
				this->add_vel_diffuse();
			}
			else if (m_visc_integ == 2) {
				this->add_vel_diffuse2(1, m_dt);
			}

			this->add_vel_source();

			this->calc_pres_rhs();

			this->solve_pres();

			this->corr_vel();

			if (step % 1 == 0) {
				TreeData savedata = TreeData::MakeCellData(*m_tree, 1 + NDIM, 0);

				TreeData::Copy(savedata, m_pres, 0, 0, 1, 0);

				VarLocInterpolator interp(*m_tree);
				for (int dir = 0; dir < NDIM; dir++) {
					interp.interp_face_to_cell(dir,
						savedata, m_vel[dir],
						1 + dir, 0, 1, 0);
				}

				WriterMeta wm;
				wm.guessName(savedata);

				std::ostringstream oss;
				oss << "test7/hoge"
					<< std::setw(4) << std::setfill('0') << step
					<< ".h5";
				WriteDataHdf5(savedata, oss.str().c_str(), wm);
			}

		}

		// 
		{
			std::cout << "saving" << std::endl;

			TreeData savedata = TreeData::MakeCellData(*m_tree, 1 + NDIM, 0);

			TreeData::Copy(savedata, m_pres, 0, 0, 1, 0);
			
			VarLocInterpolator interp(*m_tree);
			for (int dir = 0; dir < NDIM; dir++) {
				interp.interp_face_to_cell(dir,
					savedata, m_vel[dir], 
					1 + dir, 0, 1, 0);
			}

			WriterMeta wm;
			wm.guessName(savedata);
			WriteDataHdf5(savedata, "test7/hoge.h5", wm);
		}
		if (0) {
			fill_vel_all();

			TreeData node_vel = TreeData::MakeNodeData(*m_tree, NDIM, 0);

			VarLocInterpolator interp(*m_tree);
			for (int dir = 0; dir < NDIM; dir++) {
				interp.interp_face_to_node(dir,
					node_vel, m_vel[dir], 
					dir, 0, 1, 0);
			}

			WriterMeta wm;
			wm.guessName(node_vel);
			WriteDataHdf5(node_vel, "test7/vel.h5", wm);

		}
	}

protected:

	void init_tree() {

		m_tree = std::make_unique<AmrTree>();

		m_tree->minAmrLevel = m_min_level;
		m_tree->maxAmrLevel = m_max_level;

		m_tree->blockValidIndex = m_valid_box;
		m_tree->blockGrownIndex = m_grown_box;
		m_tree->blockNumGrow = m_ng;

		// allocate 
		m_tree->init();

		// setup root level
		{
			const int root_level = 1;

#define TEST7_ROOT_CONFIG 1
#if (TEST7_ROOT_CONFIG == 1)
			const int nroot = 7;
			const double xlo[][3] = {
				0, 0, 0,
				1, 0, 0,
				1, 1, 0,
				2, 1, 0,
				3, 1, 0,
				3, 0, 0,
				4, 0, 0,
			};
			const double xhi[][3] = {
				1, 1, 1,
				2, 1, 1,
				2, 2, 1,
				3, 2, 1,
				4, 2, 1,
				4, 1, 1,
				5, 1, 1,
			};
			const int conn[][6] = {
				//-100, 1, -100, -100, -100, -100,
				6, 1, -100, -100, -100, -100,
				0, -100, -100, 2, -100, -100,
				-100, 3, 1, -100, -100, -100,
				2, 4, -100, -100, -100, -100,
				3, -100, 5, -100, -100, -100,
				-100, 6, -100, 4, -100, -100,
				//5, -100, -100, -100, -100, -100,
				5, 0, -100, -100, -100, -100,
			};
#elif (TEST7_ROOT_CONFIG == 2)
			const int nroot = 5;
			const double xlo[][3] = {
				0, 0, 0,
				1, 0, 0,
				2, 0, 0,
				3, 0, 0,
				4, 0, 0,
			};
			const double xhi[][3] = {
				1, 1, 1,
				2, 1, 1,
				3, 1, 1,
				4, 1, 1,
				5, 1, 1,
			};
			const int conn[][6] = {
				4, 1, -100, -100, -100, -100,
				0, 2, -100, -100, -100, -100,
				1, 3, -100, -100, -100, -100,
				2, 4, -100, -100, -100, -100,
				3, 0, -100, -100, -100, -100,
			};
#endif
#undef TEST7_ROOT_CONFIG

			m_tree->numBlocks = nroot;

			for (int iroot = 0; iroot < nroot; iroot++) {
				auto &root = m_tree->blocks[iroot];
				
				root.setBound(xlo[iroot], xhi[iroot]);

				root.nodeType = TreeNodeType_IsLeaf;
				root.setLevel(root_level);

				for (FaceIndex iface = 0; iface < FaceIndex::NumFace; ++iface) {
					root.neighbor[iface] = conn[iroot][iface];
				}
			}
		}


	}

	void refine_tree() {

		// static refine
		m_tree->initRefineToMinLevel();

		// dynamic
		for (int ilevel = m_min_level; ilevel <= m_max_level; ilevel++) {
			for (int iblock = 0; iblock < m_tree->numBlocks; iblock++) {
				if ((*m_tree)[iblock].getLevel() == ilevel) {
					test_block_refine(iblock);
				}
			}

			m_tree->regrid();
		}

		// current finest level
		m_finest_level = m_tree->currentFinestLevel();

		std::cout << __FUNCTION__
			<< ": minlvl=" << m_min_level
			<< "; maxlvl=" << m_max_level
			<< "; finest=" << m_finest_level
			<< std::endl;
	}

	void init_data() {
		
		const int ng = m_ng;

		// buffers
		m_pres = TreeData::MakeCellData(*m_tree, 1, ng);
		m_pres.setValue(0);

		for (int dir = 0; dir < NDIM; dir++) {
			m_vel[dir] = TreeData::MakeFaceData(*m_tree, dir, 1, ng);
			m_vel[dir].setValue(0);
		}

		for (int dir = 0; dir < NDIM; dir++) {
			m_vel_tmp[dir] = TreeData::MakeFaceData(*m_tree, dir, 1, 0);
			m_vel_tmp[dir].setValue(0);
		}

		for (int dir = 0; dir < NDIM; dir++) {
			m_gradp[dir] = TreeData::MakeFaceData(*m_tree, dir, 1, ng);
			m_gradp[dir].setValue(0);
		}

		m_rhs = TreeData::MakeCellData(*m_tree, 1, 0);
		m_rhs.setValue(0);

		// BC
		m_bc_pres.reset(new BCPatch(*m_tree, m_pres, m_pres.indexType()));
		{
			BCRegister &bcreg = m_bc_pres->boundaryRegister(0);
			for (SurroundingIndex isurr = 0; isurr < SurroundingIndex::NumSurr; ++isurr) {
				bcreg.setBCMap(isurr.ipos(), isurr.jpos(), isurr.kpos(), 
					Exterior, BCType::Neumann);
			}
		}

		for (int dir = 0; dir < NDIM; dir++) {
			m_bc_vel[dir].define(*m_tree, m_vel[dir].indexType(),
				m_vel[dir].numComp(), m_vel[dir].numGrow());

			auto &bcreg = m_bc_vel[dir].boundaryRegister(0);
			for (SurroundingIndex isurr = 0; isurr < SurroundingIndex::NumSurr; ++isurr) {
				bcreg.setBCMap(isurr, Exterior, BCType::Dirichlet);
				bcreg.setBCVal(isurr, Exterior, 0.0);
			}
		}

		// fill
		for (int dir = 0; dir < NDIM; dir++) {
			m_fill_vel.emplace_back(new FillPatch(*m_tree, m_vel[dir]));
			m_fill_vel[dir]->setProlongation(0, InterpPatch::PROLONG_CENTERED_LINEAR);
			m_fill_vel[dir]->setRestriction(0, InterpPatch::RESTRICT_AVERAGE);
			//m_fill_vel[dir]->setBC(m_bc_vel[dir]);
		}

		m_fill_pres.reset(new FillPatch(*m_tree, m_pres));
		m_fill_pres->setProlongation(0, InterpPatch::PROLONG_CENTERED_LINEAR);
		m_fill_pres->setRestriction(0, InterpPatch::RESTRICT_AVERAGE);
		m_fill_pres->setBC(*m_bc_pres);
	}

	void init_solver() {
		// 
		m_solver.reset(new MacSolver(*m_tree));
		m_solver->setVerbose(1);

		//
		m_tower = &m_solver->getTower(m_finest_level);
	}

private:

	void test_block_refine(int iblock) {
		auto &block = (*m_tree)[iblock];
		const int level = block.getLevel();

		const IndexBox &validbox = m_tree->validBlockCellBox();

		if (block.isLeaf() || block.isParentOfLeafChild()) {
			const Vector3d blo = block.boundBox.lo();
			const Vector3d dh = m_tree->getBlockCellSize(iblock);
			const double dx = dh.x;

			int need_refine = 0;
			int need_coarsen = 0;

			for (int k = validbox.klo(); k <= validbox.khi(); k++) {
				for (int j = validbox.jlo(); j <= validbox.jhi(); j++) {
					for (int i = validbox.ilo(); i <= validbox.ihi(); i++) {
						double z = blo.z + ((double)k + 0.5) * dh.z;
						double y = blo.y + ((double)j + 0.5) * dh.y;
						double x = blo.x + ((double)i + 0.5) * dh.x;

						{
							double dst = hypot(x - 1, y - 1);
							if (dst < dx * 2) need_refine = 1;
						}
						{
							double dst = hypot(x - 2, y - 1);
							if (dst < dx * 2) need_refine = 1;
						}
						{
							double dst = hypot(x - 3, y - 1);
							if (dst < dx * 2) need_refine = 1;
						}
						{
							double dst = hypot(x - 4, y - 1);
							if (dst < dx * 2) need_refine = 1;
						}
					}
				}
			}

			if (level < m_tree->maxRefineLevel()) {
				if (need_refine) m_tree->refineFlag[iblock] = 1;
			}
			if (level > m_tree->minAmrLevel && !m_tree->refineFlag[iblock]) {
				if (need_coarsen) m_tree->coarsenFlag[iblock] = 1;
			}
		}
	}

	void add_vel_diffuse() {
		
		fill_vel_all();

		for (int iblock = 0; iblock < m_tree->numBlocks; iblock++) {
			if (!(*m_tree)[iblock].isLeaf()) continue;

			add_vel_diffuse_block(iblock, m_dt);
		}

		//for (int dir = 0; dir < NDIM; dir++) {
		//	m_fill_vel[dir]->restrictAll(0, 1);
		//}

		fill_vel_all();
	}

	void add_vel_diffuse2(int ilevel, double dt) {
		assert(1 <= ilevel && ilevel <= m_finest_level);

		if (ilevel == 1) {
			fill_vel_all();
		}

		// fill boundary
		fill_vel_ghost_level(ilevel, false);

		for (int iblock : m_tree->getTreeLevel(ilevel).blockIds()) {
			if (!(*m_tree)[iblock].isLeaf()) continue;
			add_vel_diffuse_block(iblock, dt);
		}

		//fill_vel_ghost_level(ilevel, false);

		if (ilevel < m_finest_level) {
			const int ncycle = 4; // subcycle for finer levels
			for (int icycle = 0; icycle < ncycle; icycle++) {
				add_vel_diffuse2(ilevel + 1, dt / ncycle);
			}
		}

		if (ilevel == 1) { 
			fill_vel_all();
		}
	}


	void add_vel_diffuse_block(int iblock, double dt) {
		const double nu = m_nu;
		//const double dt = m_dt;
		const Vector3d dh = m_tree->getBlockCellSize(iblock);

		auto &block = (*m_tree)[iblock];

		for (int dir = 0; dir < NDIM; dir++) {
			auto &u = m_vel[dir][iblock];
			auto &udiff = m_vel_tmp[dir][iblock];
			udiff.setValue(0);

			IndexBox b = get_interior_box(m_vel[dir].validBox(), block);

			BEGIN_FOR_BOX_RANGE(b, i, j, k);
			{
				double up = u(i, j, k, 0);
				udiff(i, j, k, 0) =
					(u(i + 1, j, k, 0) + u(i - 1, j, k, 0) - 2 * up) / (dh.x*dh.x) +
					(u(i, j + 1, k, 0) + u(i, j - 1, k, 0) - 2 * up) / (dh.y*dh.y) +
					(u(i, j, k + 1, 0) + u(i, j, k - 1, 0) - 2 * up) / (dh.z*dh.z);
			}
			END_FOR_BOX_RANGE;

			DoubleGridDataUtil.AddEqual(
				u, udiff, dt*nu,
				0, 0, 1, b);
		}
	}

	void add_vel_source() {

		for (int iblock = 0; iblock < m_tree->numBlocks; iblock++) {
			if (!(*m_tree)[iblock].isLeaf()) continue;

			add_vel_source_block(iblock);
		}

		//for (int dir = 0; dir < NDIM; dir++) {
		//	m_fill_vel[dir]->restrictAll(0, 1);
		//}

		fill_vel_all();
	}

	void add_vel_source_block(int iblock) {
		auto &block = (*m_tree)[iblock];

		for (int dir = 0; dir < NDIM; dir++) {
			// staggered velocity
			IndexBox b = get_interior_box(m_vel[dir].validBox(), block);
			auto &d = m_vel[dir][iblock];

			BEGIN_FOR_BOX_RANGE(b, i, j, k);
			{
				d(i, j, k, 0) += m_grav(dir) * m_dt;
			}
			END_FOR_BOX_RANGE;
		}
	}

	void calc_pres_rhs() {
		for (int iblock = 0; iblock < m_tree->numBlocks; iblock++) {
			//if (!(*m_tree)[iblock].isLeaf()) continue;

			calc_pres_rhs_block(iblock);
		}
	}

	void calc_pres_rhs_block(int iblock) {
		auto &block = (*m_tree)[iblock];
		Vector3d dh = m_tree->getBlockCellSize(iblock);

		const IndexBox &b = m_rhs.validBox();
		auto &d = m_rhs[iblock];
		d.setValue(0);

		const double dt_inv = 1.0 / m_dt;

		for (int dir = 0; dir < NDIM; dir++) {
			int ii, jj, kk;
			std::tie(ii, jj, kk) = sayaka::SelectStaggerIncr(dir);

			BEGIN_FOR_BOX_RANGE(b, i, j, k);
			{
				double ul = m_vel[dir][iblock](i, j, k, 0);
				double ur = m_vel[dir][iblock](i + ii, j + jj, k + kk, 0);
				d(i, j, k, 0) += -dt_inv / dh(dir) * (ur - ul);
			}
			END_FOR_BOX_RANGE;
		}
	}

	void solve_pres() {
		//m_solver.reset(new MacSolver(*m_tree));
		//m_solver->setVerbose(1);

		MacSolver *mac = m_solver.get();

		mac->prepareSolver();

		mac->setBC(m_bc_pres->boundaryRegister(0));

		mac->setRhs(m_rhs, 0);
		mac->setZeroSol();
		mac->setTolerance(1.0e-7, 1.0e-8);

		mac->setCGUseMGPrecond(1);

		mac->cg_solve();

		mac->getSol(m_pres, 0);

		mac->clearSolver();

		m_fill_pres->fillBoundary(0, 1, 1);

		std::cout << "solved" << std::endl;
	}

	void corr_vel() {
		const int finest_level = m_finest_level;
		for (int dir = 0; dir < NDIM; dir++) {
			m_solver->mg_flux_level(finest_level,
				dir, m_gradp[dir], m_pres, 0, 0);

			auto &tower = m_solver->getTower(finest_level);

			// u += -dt * grad(p)
			TreeData::AddEqual(m_vel[dir], m_gradp[dir], m_dt,
				0, 0, 1, 0, tower.blockIds());
		}

		fill_vel_all();
	}

	// synchronize coarse-fine face velocity, all level
	void sync_vel_all() {
		const int finest_level = m_finest_level;
		
		for (int ilevel = finest_level; ilevel >= 2; ilevel--) {
			sync_vel_level(ilevel, ilevel - 1);
		}
	}

	// synchronize velocity on fine-coarse faces between two levels
	void sync_vel_level(int fine_level, int crse_level) {
		assert(1 < fine_level && fine_level <= m_finest_level);
		assert(crse_level == fine_level - 1);

		for (int dir = 0; dir < NDIM; dir++) {
			m_tower->syncSubLevelCrseFineFlux(
				fine_level, crse_level,
				dir, m_vel[dir], 0,
				MGLevelTower::FLUX_SYNC_AVG);
		}
	}

	void avgdown_vel_all() {
		for (int ilevel = m_finest_level - 1; ilevel >= 1; ilevel--) {
			avgdown_vel_level(ilevel + 1, ilevel);
		}
	}

	void avgdown_vel_level(int fine_level, int crse_level) {
		assert(1 < fine_level && fine_level <= m_finest_level);
		assert(crse_level == fine_level - 1);

		for (int dir = 0; dir < NDIM; dir++) {
			m_fill_vel[dir]->restrictLevel(crse_level, 0, 1);
		}
	}


	void fill_vel_all() {
		const int finest_level = m_finest_level;

		// 
		for (int ilevel = finest_level - 1; ilevel >= 1; ilevel--) {
			avgdown_vel_level(ilevel + 1, ilevel);
			sync_vel_level(ilevel + 1, ilevel);
		}

		// fill ghost
		fill_vel_ghost_all();
	}

	// fill ghost points
	void fill_vel_ghost_all() {
		for (int ilevel = 1; ilevel <= m_finest_level; ilevel++) {
			fill_vel_ghost_level(ilevel);
		}
	}

	// fill ghost points
	void fill_vel_ghost_level(int ilevel, 
		bool fill_fine_crse = true)
	{
		assert(1 <= ilevel && ilevel <= m_finest_level);

		for (int dir = 0; dir < NDIM; dir++) {
			// fine-fine and fine-coarse boundary
			m_fill_vel[dir]->fillLevelBoundary(ilevel, 0, 1, 1,
				!fill_fine_crse, false, true);
			
			// physical boundary
			for (int iblock : m_tree->getTreeLevel(ilevel).blockIds()) {
				fill_vel_bc_block(iblock, dir);
			}
		}
	}

	void fill_vel_bc_block(int iblock, int dir) {
		auto &block = (*m_tree)[iblock];
		auto &vel = m_vel[dir][iblock];

		const IndexBox &validbox = m_vel[dir].validBox();
		const int ilo = validbox.ilo();
		const int ihi = validbox.ihi();
		const int jlo = validbox.jlo();
		const int jhi = validbox.jhi();
		const int klo = validbox.klo();
		const int khi = validbox.khi();

		if (dir == 0) { // u
			if (block.neighbor[0] == Exterior) {
				for (int k = klo-1; k <= khi+1; k++) {
					for (int j = jlo-1; j <= jhi+1; j++) {
						vel(ilo, j, k, 0) = 0;
					}
				}
			}
			if (block.neighbor[1] == Exterior) {
				for (int k = klo-1; k <= khi+1; k++) {
					for (int j = jlo-1; j <= jhi+1; j++) {
						vel(ihi, j, k, 0) = 0;
					}
				}
			}
			if (block.neighbor[2] == Exterior) {
				for (int k = klo-1; k <= khi+1; k++) {
					for (int i = ilo-1; i <= ihi+1; i++) {
						vel(i, jlo - 1, k, 0) = -vel(i, jlo, k, 0);
					}
				}
			}
			if (block.neighbor[3] == Exterior) {
				for (int k = klo-1; k <= khi+1; k++) {
					for (int i = ilo-1; i <= ihi+1; i++) {
						vel(i, jhi + 1, k, 0) = -vel(i, jhi, k, 0);
					}
				}
			}
			if (block.neighbor[4] == Exterior) {
				for (int j = jlo-1; j <= jhi+1; j++) {
					for (int i = ilo-1; i <= ihi+1; i++) {
						vel(i, j, klo - 1, 0) = -vel(i, j, klo, 0);
					}
				}
			}
			if (block.neighbor[5] == Exterior) {
				for (int j = jlo-1; j <= jhi+1; j++) {
					for (int i = ilo-1; i <= ihi+1; i++) {
						vel(i, j, khi + 1, 0) = -vel(i, j, khi, 0);
					}
				}
			}
		}
		else if (dir == 1) { // v
			if (block.neighbor[2] == Exterior) {
				for (int k = klo - 1; k <= khi + 1; k++) {
					for (int i = ilo - 1; i <= ihi + 1; i++) {
						vel(i, jlo, k, 0) = 0;
					}
				}
			}
			if (block.neighbor[3] == Exterior) {
				for (int k = klo - 1; k <= khi + 1; k++) {
					for (int i = ilo - 1; i <= ihi + 1; i++) {
						vel(i, jhi, k, 0) = 0;
					}
				}
			}
			if (block.neighbor[0] == Exterior) {
				for (int k = klo - 1; k <= khi + 1; k++) {
					for (int j = jlo - 1; j <= jhi + 1; j++) {
						vel(ilo - 1, j, k, 0) = -vel(ilo, j, k, 0);
					}
				}
			}
			if (block.neighbor[1] == Exterior) {
				for (int k = klo - 1; k <= khi + 1; k++) {
					for (int j = jlo - 1; j <= jhi + 1; j++) {
						vel(ihi + 1, j, k, 0) = -vel(ihi, j, k, 0);
					}
				}
			}
			if (block.neighbor[4] == Exterior) {
				for (int j = jlo - 1; j <= jhi + 1; j++) {
					for (int i = ilo - 1; i <= ihi + 1; i++) {
						vel(i, j, klo - 1, 0) = -vel(i, j, klo, 0);
					}
				}
			}
			if (block.neighbor[5] == Exterior) {
				for (int j = jlo - 1; j <= jhi + 1; j++) {
					for (int i = ilo - 1; i <= ihi + 1; i++) {
						vel(i, j, khi + 1, 0) = -vel(i, j, khi, 0);
					}
				}
			}
		}
		else if (dir == 2) { // w
			if (block.neighbor[4] == Exterior) {
				for (int j = jlo - 1; j <= jhi + 1; j++) {
					for (int i = ilo - 1; i <= ihi + 1; i++) {
						vel(i, j, klo, 0) = 0;
					}
				}
			}
			if (block.neighbor[5] == Exterior) {
				for (int j = jlo - 1; j <= jhi + 1; j++) {
					for (int i = ilo - 1; i <= ihi + 1; i++) {
						vel(i, j, khi, 0) = 0;
					}
				}
			}
			if (block.neighbor[0] == Exterior) {
				for (int k = klo - 1; k <= khi + 1; k++) {
					for (int j = jlo - 1; j <= jhi + 1; j++) {
						vel(ilo - 1, j, k, 0) = -vel(ilo, j, k, 0);
					}
				}
			}
			if (block.neighbor[1] == Exterior) {
				for (int k = klo - 1; k <= khi + 1; k++) {
					for (int j = jlo - 1; j <= jhi + 1; j++) {
						vel(ihi + 1, j, k, 0) = -vel(ihi, j, k, 0);
					}
				}
			}
			if (block.neighbor[2] == Exterior) {
				for (int k = klo - 1; k <= khi + 1; k++) {
					for (int i = ilo - 1; i <= ihi + 1; i++) {
						vel(i, jlo - 1, k, 0) = -vel(i, jlo, k, 0);
					}
				}
			}
			if (block.neighbor[3] == Exterior) {
				for (int k = klo - 1; k <= khi + 1; k++) {
					for (int i = ilo - 1; i <= ihi + 1; i++) {
						vel(i, jhi + 1, k, 0) = -vel(i, jhi, k, 0);
					}
				}
			}
		}

	}

	IndexBox get_interior_box(
		const IndexBox &validbox, 
		const AmrTreeNode &block) 
	{
		IndexBox b = validbox;

		if (b.isFaceBox()) {
			int dir = b.type().getFaceVarDir();
			if (dir == 0) {
				if (block.neighbor[0] == Exterior) b.shrinkLo(dir, 1);
				if (block.neighbor[1] == Exterior) b.shrinkHi(dir, 1);
			}
			if (dir == 1) {
				if (block.neighbor[2] == Exterior) b.shrinkLo(dir, 1);
				if (block.neighbor[3] == Exterior) b.shrinkHi(dir, 1);
			}
			if (dir == 2) {
				if (block.neighbor[4] == Exterior) b.shrinkLo(dir, 1);
				if (block.neighbor[5] == Exterior) b.shrinkHi(dir, 1);
			}
		}

		return b;
	}
};



