#pragma once

#include <cmath>

#include <iostream>
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

#define TESTDIR(path) ("test6" path)

/**
 * Droplet (or bubble) 
 * with variable density and IB plane
 */
class Test6
{
	enum {
		PhysBC_Periodic = -100,
		PhysBC_Fixed,
		PhysBC_Symmetry,
	};

	enum CellComp {
		Comp_Phi = 0,
		Comp_Divu,
		Comp_Resid,
		Comp_LS,
		Comp_IB,
		Comp_Exact,
		NumComp,
	};
	enum FaceComp {
		FaceCompVel = 0,
		FaceCompAdv,
		FaceCompRho,
		NumFaceComp,
	};
	enum NodeComp {
		NodeCompP,
		NodeCompLS,
		NodeCompU,
		NumNodeComp = NodeCompU + NDIM,
	};

	enum ProbType {
		PROB_CIRCLE,
		//PROB_STAR,
		//PROB_ELLIPSE,
	};

	AmrTree *m_tree;
	TreeData *m_data;
	BCPatch *m_bcpatch;
	FillPatch *m_fillpatch;
	TreeStateData *m_state;

	TreeData *m_face_data[NDIM];
	FillPatch *m_face_fillpatch[NDIM];

	// for visualization
	TreeData *m_node_data;

	TreeData *m_flux[NDIM];

	VarLocInterpolator *m_varloc_interp;

	int m_step;
	double m_time;

	Vector3d probmin, probmax;
	RealBox probbox;
	int probbc[NDIM*2];
	bool probperiodic[NDIM * 2];

	int probtype;
	bool prob_use_shape;
	double prob_circ_rad;

	//
	bool prob_use_density;
	double dens_liq, dens_gas;

	int n0x, n0y, n0z; // tile at root level

	int nb;
	int nbx, nby, nbz;
	int ng;
	int ngx, ngy, ngz;
	IndexBox validBox, grownBox;

	int minlevel, maxlevel;

public:

	void init(int argc, char *argv[]) {
		m_step = 0;
		m_time = 0;

		probmin = { -0.5, -0.5, 0.0 };
		probmax = { 0.5, 0.5, 0.0 };
		if (NDIM == 3) {
			probmin.z = -0.5;
			probmax.z = 0.5;
		}
		probbox = RealBox(probmin, probmax);

		// problem
		probtype = PROB_CIRCLE;
		//probtype = PROB_STAR;
		//probtype = PROB_ELLIPSE;
		prob_circ_rad = 0.2;

		//prob_use_shape = false;
		prob_use_shape = true;

		prob_use_density = true;
		//prob_use_density = false;
		dens_liq = 1.0e3;
		dens_gas = 1.0;

		for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
			//probbc[face] = PhysBC_Periodic;
			//probbc[face] = PhysBC_Symmetry;
			probbc[face] = PhysBC_Fixed;
			probperiodic[face] = false;
		}

		//
		n0x = 1;
		n0y = 1;
		n0z = 1;

		const int ndymref = 3;
		//minlevel = 1;
		minlevel = 3;
		//maxlevel = 3;
		//maxlevel = 7;
		//minlevel = 8;
		maxlevel = minlevel + ndymref;

		nb = 4;
		assert(nb%2 == 0);
		//
		nbx = nb;
		nby = nb;
		nbz = NDIM<=2 ? 1 : nb;
		
		//
		ng = 2;
		ngx = ng;
		ngy = ng;
		ngz = NDIM<2 ? 0 : ng;

		validBox = IndexBox(
			{ 0, 0, 0 },
			{ nbx - 1, nby - 1, nbz - 1 });
		grownBox = IndexBox(
			IndexBox(validBox).extend(ngx,ngy,ngz));

		//
		this->init_tree();

		//
		this->init_data();
		//
		this->init_patch();
		// 
		this->init_state();

		//this->move_umac_to_unode();
		if (NDIM == 2) {
		// save data
		m_tree->writeTreeLeafBlockGrid(TESTDIR("/amr_leaf_grid_init.tec.dat"), m_step+1, m_time);
		m_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_data_init.tec.dat"), m_step+1, m_time);
		m_node_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_nodedata_init.tec.dat"), m_step+1, m_time);
		}

		// dynamic refinement
		for (int level=minlevel; level<maxlevel+5; level++) {
			for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
				if ((*m_tree)[iblock].getLevel() == level) {
					this->test_block_refine(iblock);
				}
			}

			m_tree->regrid();

			//
			m_tree->initNewBlockDataAfterRegrid();

			for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
				calc_block_ls_sign(iblock);
			}

			//m_fillpatch->restrictAll(Comp_Phi, NumComp);
			for (int dir=0; dir<NDIM; dir++) {
				m_face_fillpatch[dir]->restrictAll(0, NumFaceComp);
			}
		}

		if (1) {
			// fill boundary
			m_fillpatch->fillBoundary(0, NumComp, 2);
			// calculate face data
			for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
				for (int dir=0; dir<NDIM; dir++) {
					calc_block_mac_value(dir, iblock);
				}
			}
		}


		// calculate div(u) cell data
		for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
			calc_block_rhs_value(iblock);
		}
		// calculate exact solution
		for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
			calc_block_exact(iblock);
		}
		
		// fill boundary
		m_fillpatch->fillBoundary(0, NumComp, 2);

		//this->move_umac_to_unode();
		// save data
		if (NDIM == 2) {
		m_tree->writeTreeLeafBlockGrid(TESTDIR("/amr_leaf_grid.tec.dat"), m_step+1, m_time);
		m_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_data.tec.dat"), m_step+1, m_time);
		m_node_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_nodedata.tec.dat"), m_step+1, m_time);

		for (int level=1; level<=maxlevel; level++) {
			char filename[64];
			sprintf(filename, TESTDIR("/amr_level%02d_nodedata.tec.dat"),level);
			m_node_data->writeTreeLevelBlockData(level, filename, m_step+1, m_time);
		}
		}

		// solving...
		if (1) { 
			std::unique_ptr<MacSolver> mac_op(new MacSolver(*m_tree));
			mac_op->setVerbose(1);

			mac_op->prepareSolver();

			if (prob_use_shape) {
				mac_op->setIBFracByLS(*m_data, Comp_IB);
			}
			if (prob_use_density) {
				std::vector<TreeData*> rho_inv;
				for (int dir=0; dir<NDIM; dir++) {
					rho_inv.push_back(m_face_data[dir]);
				}
				mac_op->setCoefB(rho_inv, FaceCompRho);
			}

			// set BC
			mac_op->setBC(m_bcpatch->boundaryRegister(Comp_Phi));
			// correct RHS with BC
			TreeData::SetValue(*m_data, 0.0, Comp_Phi, 1, 0);
			m_fillpatch->fillBoundary(Comp_Phi, 1, 1);
			mac_op->fixBC(*m_data, *m_data, Comp_Phi, Comp_Divu);

			mac_op->setRhs(*m_data, Comp_Divu);
			mac_op->setZeroSol();
			mac_op->setTolerance(0.0, 1.0e-8);
			
			mac_op->setCGUseMGPrecond(1);

			const int solver = 1;
			if (solver == 0) {
				mac_op->relax_solve();
			} else if (solver == 1) {
				mac_op->cg_solve();
			}
			else {
				mac_op->mg_solve();
			}

			// get solution
			mac_op->getSol(*m_data, Comp_Phi);

			// 
			mac_op->clearSolver();

			if (0) { // Neumann BC, shift the solution
				double sumsol = 0;
				double sumana = 0;
				double volume = 0;

				const IndexBox &vbox = m_tree->validBlockCellBox();

				for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
					if ((*m_tree)[iblock].isLeaf()) {
						const DoubleBlockData &db = (*m_data)[iblock];

						double dv = m_tree->getBlockCellVolume(iblock);

						for (int k=vbox.klo(); k<=vbox.khi(); k++) {
						for (int j=vbox.jlo(); j<=vbox.jhi(); j++) {
						for (int i=vbox.ilo(); i<=vbox.ihi(); i++) {
							if (!prob_use_shape || db(i,j,k,Comp_IB)>0) {
								// only consider valid solution
								sumsol += db(i,j,k,Comp_Phi) * dv;
								sumana += db(i,j,k,Comp_Exact) * dv;
								volume += dv;
							}
						}
						}
						}
					}
				}

				double avgsol = sumsol / volume;
				double avgana = sumana / volume;

				for (int i=0; i<m_tree->numBlocks; i++) {
					if ((*m_tree)[i].isLeaf()) {
						DoubleGridDataUtil.AddEqual(
							(*m_data)[i], Comp_Phi, 
							m_tree->validBlockCellBox(), 
							avgana-avgsol);
					}
				}


				// calculate the error
				double err_Linf = 0;
				double err_L1 = 0;
				double err_L2 = 0;
				for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
					if ((*m_tree)[iblock].isLeaf()) {
						const DoubleBlockData &db = (*m_data)[iblock];

						double dv = m_tree->getBlockCellVolume(iblock);

						for (int k=vbox.klo(); k<=vbox.khi(); k++) {
						for (int j=vbox.jlo(); j<=vbox.jhi(); j++) {
						for (int i=vbox.ilo(); i<=vbox.ihi(); i++) {
							if (!prob_use_shape || db(i,j,k,Comp_IB)>0) {
								double sol = db(i,j,k,Comp_Phi);
								double ana = db(i,j,k,Comp_Exact);
								double err = abs(sol - ana);

								err_Linf = std::max(err_Linf, err);
								err_L1 += err * dv;
								err_L2 += err*err * dv;
							}
						}
						}
						}
					}
				}
				
				//
				err_L1 /= volume;
				err_L2 = sqrt(err_L2 / volume);

				LOGPRINTF("level_min=%d, level_max=%d\n", minlevel, maxlevel);
				LOGPRINTF("Linf=%e, L1=%e, L2=%e\n", err_Linf, err_L1, err_L2);
			}

			// now the data only lives on leaf nodes
			// have to average down to parents
			m_fillpatch->restrictAll(Comp_Phi, 1);
			move_pcell_to_pnode();

			if (0) { // flux 
				const int finest_level = m_tree->currentFinestLevel();
				for (int dir=0; dir<NDIM; dir++) {
					mac_op->mg_flux_level(finest_level, 
						dir, *(m_flux[dir]), *m_data,
						0, Comp_Phi);
				}

				for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
					if ((*m_tree)[iblock].isLeaf()) {
						calc_block_flux_resid(iblock);
					}
				}
			}

			// save data
			if (NDIM == 2) {
			//m_tree->writeTreeLeafBlockGrid(TESTDIR("/amr_leaf_grid.tec.dat"), m_step+1, m_time);
			m_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_data_sol.tec.dat"), m_step+1, m_time);
			m_node_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_nodedata_sol.tec.dat"), m_step+1, m_time);
			}

			WriterMeta cellmeta;
			cellmeta.guessName(*m_data);
			cellmeta.name[Comp_Phi] = "pres";
			cellmeta.name[Comp_Divu] = "rhs";
			cellmeta.name[Comp_Resid] = "resid";
			cellmeta.name[Comp_LS] = "ls";
			cellmeta.name[Comp_IB] = "sdf";
			cellmeta.name[Comp_Exact] = "pexact";
			//WriteLeafDataVtk(*m_data, TESTDIR("/cell_leaf_sol.vtu"), cellmeta);
			WriteDataHdf5(*m_data, TESTDIR("/cell_leaf_sol.h5"), cellmeta);

			WriterMeta nodemeta;
			nodemeta.guessName(*m_node_data);
			nodemeta.name[NodeCompP] = "pres";
			//WriteLeafDataVtk(*m_node_data, TESTDIR("/node_leaf_sol.vtu"), nodemeta);
			WriteDataHdf5(*m_node_data, TESTDIR("/node_leaf_sol.h5"), nodemeta);
		}
	}


protected:

	void init_tree() {
		m_tree = new AmrTree;
		//
		m_tree->minAmrLevel = minlevel;
		m_tree->maxAmrLevel = maxlevel;
		//
		m_tree->blockValidIndex = validBox;
		m_tree->blockGrownIndex = grownBox;
		m_tree->blockNumGrow = ng;
		
		// allocate memory
		m_tree->init();

		if (0) { // setup root level
			const int root_level = 1;
			const int nroot = n0x * n0y * n0z;

			m_tree->numBlocks = nroot;

			for (int k=0; k<n0z; k++) {
			for (int j=0; j<n0y; j++) {
			for (int i=0; i<n0x; i++) {
				int iroot = i + j*n0x + k*n0x*n0y;

				AmrTreeNode &root = m_tree->blocks[iroot];

				Vector3d block_lo = probbox.lo();
				Vector3d block_hi = probbox.lo();
				block_lo.x += probbox.length(0) / n0x * i;
				block_hi.x += probbox.length(0) / n0x * (i+1);
				block_lo.y += probbox.length(1) / n0y * j;
				block_hi.y += probbox.length(1) / n0y * (j+1);
				if (NDIM == 3) {
					block_lo.z += probbox.length(2) / n0z * k;
					block_hi.z += probbox.length(2) / n0z * (k+1);
				}

				RealBox rootbox(block_lo, block_hi);
				root.boundBox = rootbox;
				root.blockCenter = rootbox.center();
				root.blockLength = rootbox.length();

				root.nodeType = TreeNodeType_IsLeaf;
				root.setLevel(root_level);

				// set BC to all-periodic
				for (FaceIndex iface=0; iface<FaceIndex::NumFace; iface++) {
					int ineigh = -1;
					if (iface == 0) {
						if (i == 0) {
							if (probbc[iface] == PhysBC_Periodic) {
								ineigh = (n0x-1) + j*n0x + k*n0x*n0y;
							} else {
								ineigh = probbc[iface];
							}
						} else {
							ineigh = (i-1) + j*n0x + k*n0x*n0y;
						}
					} else if (iface == 1) {
						if (i == n0x-1) {
							if (probbc[iface] == PhysBC_Periodic)
								ineigh = 0 + j*n0x + k*n0x*n0y;
							else
								ineigh = probbc[iface];
						} else {
							ineigh = (i+1) + j*n0x + k*n0x*n0y;
						}
					} else if (iface == 2) {
						if (j == 0) {
							if (probbc[iface] == PhysBC_Periodic) 
								ineigh = i + (n0y-1)*n0x + k*n0x*n0y;
							else
								ineigh = probbc[iface];
						} else {
							ineigh = i + (j-1)*n0x + k*n0x*n0y;
						}
					} else if (iface == 3) {
						if (j == n0y-1) {
							if (probbc[iface] == PhysBC_Periodic)
								ineigh = i + (0)*n0x + k*n0x*n0y;
							else
								ineigh = probbc[iface];
						} else {
							ineigh = i + (j+1)*n0x + k*n0x*n0y;
						}
					} else if (iface == 4) {
						if (k == 0) {
							if (probbc[iface] == PhysBC_Periodic)
								ineigh = i + j*n0x + (n0z-1)*n0x*n0y;
							else
								ineigh = probbc[iface];
						} else {
							ineigh = i + j*n0x + (k-1)*n0x*n0y;
						}
					} else if (iface == 5) {
						if (k == n0z-1) {
							if (probbc[iface] == PhysBC_Periodic)
								ineigh = i + j*n0x + (0)*n0x*n0y;
							else
								ineigh = probbc[iface];
						} else {
							ineigh = i + j*n0x + (k+1)*n0x*n0y;
						}
					}
					//assert(0<=ineigh && ineigh<n0x*n0y*n0z);
					root.neighbor[iface] = ineigh;
				} // end loop face
			}
			}
			}
		}
		else {
			m_tree->initUniformRootLevel({ n0x,n0y,n0z }, probbox, probbc, probperiodic);
		}

		// static refine up to min level
		m_tree->initRefineToMinLevel();
	}

	// must be called after tree initialization
	void init_data() {
		// cell data
		m_data = new TreeData(*m_tree, validBox, grownBox, ng, NumComp);
		m_data->setValue(0);
		
		// set LS sign
		for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
			calc_block_ls_sign(iblock);
		}

		// face-based data
		for (int dir=0; dir<NDIM; dir++) {
			IndexBox faceBox = IndexBox::Stagger(validBox, dir);
			const int face_ngrow = ng;
			IndexBox faceGrownBox = IndexBox::Extend(faceBox, face_ngrow);

			m_face_data[dir] = new TreeData(
				*m_tree, faceBox, faceGrownBox, face_ngrow, NumFaceComp);
			m_face_data[dir]->setValue(0);
		}

		// set value for MAC data
		for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
			for (int dir=0; dir<NDIM; dir++) {
				calc_block_mac_value(dir, iblock);
			}
		}

		// calculate div(u) cell data
		for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
			calc_block_rhs_value(iblock);
		}

		// node-centered data
		{
			IndexBox nodeBox = validBox;
			nodeBox.staggerAll();
			const int node_ngrow = 0;
			m_node_data = new TreeData(
				*m_tree, nodeBox, nodeBox, node_ngrow, NumNodeComp);
			m_node_data->setValue(0.0);
		}

		// flux 
		for (int dir=0; dir<NDIM; dir++) {
			IndexBox facebox = IndexBox::Stagger(validBox, dir);
			m_flux[dir] = new TreeData(
				*m_tree, facebox, facebox, 0, 1);
			m_flux[dir]->setValue(0.0);
		}
	}

	void init_patch() {
		m_bcpatch = new BCPatch(*m_tree, *m_data, VARLOC_CELL);
		for (int comp=0; comp<NumComp; comp++) {
			BCRegister& bcreg = m_bcpatch->boundaryRegister(comp);
			for (SurroundingIndex isurr=0; isurr<SurroundingIndex::NumSurr; isurr++) {
				int ipos = isurr.ipos();
				int jpos = isurr.jpos();
				int kpos = isurr.kpos();

				if (comp == Comp_Phi) {
					bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Fixed, Dirichlet);
					bcreg.setBCVal(ipos, jpos, kpos, PhysBC_Fixed, 10.0);
					bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Symmetry, Neumann);
				} else {
					bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Symmetry, Neumann);
					bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Fixed, Neumann);
				}
			}
		}

		m_fillpatch = new FillPatch(*m_tree, *m_data/*, VARLOC_CELL*/);
		for (int comp=0; comp<NumComp; comp++) {
			int prolong = InterpPatch::PROLONG_CENTERED_LINEAR;
			int restrict = InterpPatch::RESTRICT_AVERAGE;

			if (comp == Comp_Divu) {
				prolong = InterpPatch::PROLONG_IGNORE;
			} else if (comp == Comp_LS) {
				//prolong = InterpPatch::PROLONG_INJECTION;
				prolong = InterpPatch::PROLONG_CENTERED_LINEAR;
			}

			m_fillpatch->setProlongation(comp, prolong);
			m_fillpatch->setRestriction(comp, restrict);
		}
		// bind BC
		m_fillpatch->bc_patch = m_bcpatch;

		// face fill-patch
		for (int dir=0; dir<NDIM; dir++) {
			m_face_fillpatch[dir] = new FillPatch(*m_tree, *m_face_data[dir]);
			for (int comp=0; comp<NumFaceComp; comp++) {
				m_face_fillpatch[dir]->setProlongation(comp, InterpPatch::PROLONG_CENTERED_LINEAR);
				m_face_fillpatch[dir]->setRestriction(comp, InterpPatch::RESTRICT_AVERAGE);
			}
			// we do not need BC for this problem
			m_face_fillpatch[dir]->bc_patch = NULL;
		}

		//
		m_varloc_interp = new VarLocInterpolator(*m_tree);
	}

	void init_state() {
		// cell state
		m_state = new TreeStateData;
		// set variables
		m_state->new_data = m_data;
		m_state->bc_patch = m_bcpatch;
		m_state->fill_new_patch = m_fillpatch;
		// add to tree
		m_tree->treeStateData.push_back(*m_state);

		// MAC state
		if (0) {
		for (int dir=0; dir<NDIM; dir++) {
			TreeStateData face_state;
			face_state.new_data = m_face_data[dir];
			face_state.bc_patch = NULL;
			face_state.fill_new_patch = m_face_fillpatch[dir];

			m_tree->treeStateData.push_back(face_state);
		}
		}
	}

	void test_block_refine(int iblock) {
		AmrTreeNode &block = (*m_tree)[iblock];
		int level = block.getLevel();

		const IndexBox &validbox = m_tree->validBlockCellBox();

		const DoubleBlockData &data = (*m_data)[iblock];
		const int lscomp = Comp_LS;
		const int ibcomp = Comp_IB;

		if (block.isLeaf() || block.isParentOfLeafChild()) {
			const RealBox &bbox = block.boundBox;
			const Vector3d &blo = bbox.lo();
			const Vector3d dh = m_tree->getBlockCellSize(iblock);
			const double dx = dh.x;

			int needRefine = 0;
			int needCoarsen = 0;

			for (int k=validbox.klo(); k<=validbox.khi(); k++) {
			for (int j=validbox.jlo(); j<=validbox.jhi(); j++) {
			for (int i=validbox.ilo(); i<=validbox.ihi(); i++) {
				double z = blo.z + ((double) ZDIM) * ((double) k + 0.5) * dh.z;
				double y = blo.y + ((double) j + 0.5) * dh.y;
				double x = blo.x + ((double) i + 0.5) * dh.x;

				double ls = data(i,j,k,lscomp);
				double ib = data(i,j,k,ibcomp);
				const double lscutoff = dx * nb;
				const double ibcutoff = dx * nb;
				if (abs(ls)<lscutoff) {
					needRefine = 1;
				}
				if (prob_use_shape && abs(ib)<ibcutoff) {
					needRefine = 1;
				}
			}
			}
			}


			if (level < m_tree->maxRefineLevel()) {
				if (needRefine) {
					m_tree->refineFlag[iblock] = 1;
				}
			}
			if (level>m_tree->minAmrLevel && !m_tree->refineFlag[iblock]) {
				if (needCoarsen) {
					m_tree->coarsenFlag[iblock] = 1;
				}
			}
		}
	}

	// assign LS value (only sign is valid)
	void calc_block_ls_sign(int iblock) {
		AmrTreeNode &block = (*m_tree)[iblock];
		int level = block.getLevel();

		const IndexBox &validbox = m_tree->validBlockCellBox();

		DoubleBlockData &ls = (*m_data)[iblock];
		const int lscomp = Comp_LS;
		const int ibcomp = Comp_IB;

		{
			const RealBox &bbox = block.boundBox;
			const Vector3d &blo = bbox.lo();
			const Vector3d dh = m_tree->getBlockCellSize(iblock);

			for (int k=validbox.klo(); k<=validbox.khi(); k++) {
			for (int j=validbox.jlo(); j<=validbox.jhi(); j++) {
			for (int i=validbox.ilo(); i<=validbox.ihi(); i++) {
				double z = blo.z + ((double) ZDIM) * ((double) k + 0.5) * dh.z;
				double y = blo.y + ((double) j + 0.5) * dh.y;
				double x = blo.x + ((double) i + 0.5) * dh.x;

				Vector3d pt { x, y, z };
				
				// interface LS
				double val = test_ls_sign(pt);
				ls(i,j,k,lscomp) = val;

				// IB LS
				val = test_ib_sign(pt);
				ls(i,j,k,ibcomp) = val;
			}
			}
			}
		}
	}

	void calc_block_exact(int iblock) {
		AmrTreeNode &block = (*m_tree)[iblock];
		int level = block.getLevel();

		const IndexBox &validbox = m_tree->validBlockCellBox();

		DoubleBlockData &exact = (*m_data)[iblock];
		const int lscomp = Comp_Exact;

		{
			const RealBox &bbox = block.boundBox;
			const Vector3d &blo = bbox.lo();
			const Vector3d dh = m_tree->getBlockCellSize(iblock);

			for (int k=validbox.klo(); k<=validbox.khi(); k++) {
			for (int j=validbox.jlo(); j<=validbox.jhi(); j++) {
			for (int i=validbox.ilo(); i<=validbox.ihi(); i++) {
				double z = blo.z + ((double) ZDIM) * ((double) k + 0.5) * dh.z;
				double y = blo.y + ((double) j + 0.5) * dh.y;
				double x = blo.x + ((double) i + 0.5) * dh.x;

				//Vector3d pt = Vector3d::VecMake(x, y, z);
				
				double val = 0.25 * (cos(4.0*M_PI*x) + cos(4.0*M_PI*y));

				exact(i,j,k,lscomp) = val;
			}
			}
			}
		}
	}

	// calculate the div(uadv)
	void calc_block_rhs_value(int iblock) {
		const IndexBox &vbox = m_tree->validBlockCellBox();
		const Vector3d dh = m_tree->getBlockCellSize(iblock);

		const int comp = Comp_Divu;

		DoubleBlockData &divu = (*m_data)[iblock];
		divu.setValue(comp, 1, 0.0);

		for (int dir=0; dir<NDIM; dir++) {
			int ii, jj, kk;
			sayaka::SelectStaggerIncr(dir, ii,jj,kk);

			const DoubleBlockData &umac = (*(m_face_data[dir]))[iblock];
			const int maccomp = FaceCompAdv;
			const int rhocomp = FaceCompRho; // in fact 1/rho

			for (int k=vbox.klo(); k<=vbox.khi(); k++) {
				for (int j=vbox.jlo(); j<=vbox.jhi(); j++) {
					for (int i=vbox.ilo(); i<=vbox.ihi(); i++) {
						double ul = umac(i,j,k,maccomp);
						double ur = umac(i+ii,j+jj,k+kk,maccomp);
						if (prob_use_density) {
							ul *= umac(i,j,k,rhocomp);
							ur *= umac(i+ii,j+jj,k+kk,rhocomp);
						}

						divu(i,j,k,comp) += 1.0/dh(dir) * (ur-ul);
					}
				}
			}
		}
	}

	// residual after flux correction
	void calc_block_flux_resid(int iblock) {
		const int residcomp = Comp_Resid;
		const int maccomp = FaceCompAdv;
		const int fluxcomp = 0;

		const IndexBox &vbox = m_tree->validBlockCellBox();
		const Vector3d dh = m_tree->getBlockCellSize(iblock);

		DoubleBlockData &buf = (*m_data)[iblock];
		buf.setValue(residcomp, 1, 0.0);

		for (int dir=0; dir<NDIM; dir++) {
			int ii, jj, kk;
			sayaka::SelectStaggerIncr(dir, ii,jj,kk);

			const DoubleBlockData &flux = (*(m_flux[dir]))[iblock];
			const DoubleBlockData &umac = (*(m_face_data[dir]))[iblock];

			for (int k=vbox.klo(); k<=vbox.khi(); k++) {
				for (int j=vbox.jlo(); j<=vbox.jhi(); j++) {
					for (int i=vbox.ilo(); i<=vbox.ihi(); i++) {
						double ul = umac(i,j,k,maccomp) - flux(i,j,k,fluxcomp);
						double ur = umac(i+ii,j+jj,k+kk,maccomp) - flux(i+ii,j+jj,k+kk,fluxcomp);
						//double ul = umac(i,j,k,maccomp);
						//double ur = umac(i+ii,j+jj,k+kk,maccomp);
						//double ul = flux(i,j,k,fluxcomp);
						//double ur = flux(i+ii,j+jj,k+kk,fluxcomp);
						buf(i,j,k,residcomp) += 1.0/dh(dir) * (ur - ul);
					}
				}
			}
		}
	}

	void calc_block_mac_value(int dir, int iblock) {
		const Vector3i nb = m_tree->validBlockCellBox().size();
		const RealBox &bbox = (*m_tree)[iblock].boundBox;
		const Vector3d dh = m_tree->getBlockCellSize(iblock);
		const double dx = bbox.length(0) / nb.x;
		const double dy = bbox.length(1) / nb.y;
		const double dz = bbox.length(2) / nb.z;
		const double &xlo = bbox.lo().x;
		const double &ylo = bbox.lo().y;
		const double &zlo = bbox.lo().z;

		{
			int ii, jj, kk;
			sayaka::SelectStaggerIncr(dir, ii,jj,kk);

			TreeData &face_data = *(m_face_data[dir]);

			const IndexBox &face_box = face_data.validBox();

			DoubleBlockData &block_data = face_data[iblock];

			const DoubleBlockData &ls = (*m_data)[iblock];
			const int lscomp = Comp_LS;
			const int ibcomp = Comp_IB;

			for (int k=face_box.klo(); k<=face_box.khi(); k++) {
			for (int j=face_box.jlo(); j<=face_box.jhi(); j++) {
			for (int i=face_box.ilo(); i<=face_box.ihi(); i++) {
				double z = zlo + ((double) k + 0.5*(1-kk)) * dz;
				double y = ylo + ((double) j + 0.5*(1-jj)) * dy;
				double x = xlo + ((double) i + 0.5*(1-ii)) * dx;

				for (int comp=0; comp<NumFaceComp; comp++) {
					//double val = calc_mac_value(dir, comp, x,y,z);
					double val = 0;

					if (comp == FaceCompAdv) {
						if (probtype==PROB_CIRCLE) {
							const double rad = prob_circ_rad;
							const double kappa = 1.0/rad * (NDIM-1);
							double hl = ls(i-ii,j-jj,k-kk,lscomp)>=0 ? 1 : 0;
							double hr = ls(i,j,k,lscomp)>=0 ? 1 : 0;
							val = -kappa * (hr-hl) / dh(dir);
						}

						if (prob_use_shape) {
							if (ls(i-ii,j-jj,k-kk,ibcomp)<=0 || ls(i,j,k,ibcomp)<=0) {
								val = 0;
							}
						}
					} else if (comp == FaceCompRho) {
						double phil = ls(i-ii,j-jj,k-kk,lscomp);
						double phir = ls(i,j,k,lscomp);
						double theta = height_frac(phil, phir);
						double rho = dens_liq*theta + dens_gas*(1.0-theta);
						val = 1.0 / rho;
					}

					
					block_data(i,j,k,comp) = val;
				}
			}
			}
			}
		}
	}

	//double calc_mac_value(int dir, int comp, double x, double y, double z) const {
	//	const double k = 2.0 * M_PI;
	//	
	//	const double skx = sin(k*x);
	//	const double ckx = cos(k*x);
	//	const double sky = sin(k*y);
	//	const double cky = cos(k*y);

	//	const double u = skx * cky;
	//	const double v = -ckx * sky;
	//	const double w = 0;

	//	double val = 0;
	//	if (comp == FaceCompVel) {
	//		if (dir == 0) {
	//			val = u;
	//		} else if (dir == 1) { 
	//			val = v;
	//		} else {
	//			val = w;
	//		}
	//	} else if (comp == FaceCompAdv) {
	//		if (dir == 0) {
	//			const double dudx = k * ckx * cky;
	//			const double dudy = -k * skx * sky;
	//			val = u*dudx + v*dudy;
	//		} else if (dir == 1) {
	//			const double dvdx = k * skx * sky;
	//			const double dvdy = -k * ckx * cky;
	//			val = u*dvdx + v*dvdy;
	//		} else {
	//			val = 0;
	//		}
	//	}
	//	return val;
	//}

	void move_umac_to_unode() {
		for (int dir=0; dir<NDIM; dir++) {
			m_face_fillpatch[dir]->fillBoundary(FaceCompVel, 1, 1);
			
			for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
				//
				m_varloc_interp->interp_block_face_to_node(dir,
					(*m_node_data)[iblock], (*m_face_data[dir])[iblock],
					m_node_data->validBox(), 
					NodeCompU+dir, FaceCompVel);
			}
		}
	}
	void move_pcell_to_pnode() {
		m_fillpatch->fillBoundary(Comp_Phi, 1, 2);
		for (int i=0; i<m_tree->numBlocks; i++) {
			m_varloc_interp->interp_block_cell_to_node(
				(*m_node_data)[i], (*m_data)[i], 
				m_node_data->validBox(),
				NodeCompP, Comp_Phi);
		}

		m_fillpatch->fillBoundary(Comp_LS, 1, 2);
		for (int i=0; i<m_tree->numBlocks; i++) {
			m_varloc_interp->interp_block_cell_to_node(
				(*m_node_data)[i], (*m_data)[i], 
				m_node_data->validBox(),
				NodeCompLS, Comp_LS);
		}
	}

	double test_ls_sign(const Vector3d &ptin) const {
		Vector3d pos = ptin;

		double ls = 0;
		if (probtype == PROB_CIRCLE) {
			const Vector3d cen = Vector3d::VecMake(0.0, 0.0, 0.0);
			const double rad = prob_circ_rad;
			double r = Vec3Length(pos, cen);
			ls = rad - r;
		} 
		//else if (probtype == PROB_STAR) {
		//	pos.z = 0;
		//	// r(theta) = ra + rb*cos(omega*theta)
		//	const double ra = 0.237;
		//	const double rb = 0.079;
		//	const double omega = 6.0;
		//	//static const Vector3d cen = Vector3d::VecMake(0.0, 0.0, 0.0);

		//	double theta = atan2(pos.y, pos.x);
		//	double rstar = ra + rb*cos(omega*theta);
		//	double rpt = Vec3Length(pos);
		//	ls = rpt - rstar;
		//} else if (probtype == PROB_ELLIPSE) {
		//	pos.z = 0;

		//	const double a = 0.75 / 2;
		//	const double b = 0.625 / 2;

		//	double theta = atan2(pos.y, pos.x);
		//	double ct = cos(theta);
		//	double st = sin(theta);
		//	double rellipse = sqrt(1.0 / (ct*ct/(a*a) + st*st/(b*b)));
		//	double rpt = Vec3Length(pos);
		//	ls = rellipse - rpt;
		//}
		else {
			LOGPRINTF("Unknown probtype=%d\n", probtype);
			exit(1);
		}
		return ls;
	}

	double test_ib_sign(const Vector3d &ptin) const {
		Vector3d pos = ptin;

		double ls = 0;
		if (probtype == PROB_CIRCLE) {
			//const Vector3d cen = Vector3d::VecMake(0.0, 0.0, 0.0);
			const double rad = prob_circ_rad;
			const double yplane = -rad / 2;
			ls = pos.y - yplane;
		} 
		//else if (probtype == PROB_STAR) {
		//	pos.z = 0;
		//	// r(theta) = ra + rb*cos(omega*theta)
		//	const double ra = 0.237;
		//	const double rb = 0.079;
		//	const double omega = 6.0;
		//	//static const Vector3d cen = Vector3d::VecMake(0.0, 0.0, 0.0);

		//	double theta = atan2(pos.y, pos.x);
		//	double rstar = ra + rb*cos(omega*theta);
		//	double rpt = Vec3Length(pos);
		//	ls = rpt - rstar;
		//} else if (probtype == PROB_ELLIPSE) {
		//	pos.z = 0;

		//	const double a = 0.75 / 2;
		//	const double b = 0.625 / 2;

		//	double theta = atan2(pos.y, pos.x);
		//	double ct = cos(theta);
		//	double st = sin(theta);
		//	double rellipse = sqrt(1.0 / (ct*ct/(a*a) + st*st/(b*b)));
		//	double rpt = Vec3Length(pos);
		//	ls = rellipse - rpt;
		//}
		else {
			LOGPRINTF("Unknown probtype=%d\n", probtype);
			exit(1);
		}
		return ls;
	}

	double height_frac(double phil, double phir) const {
		double theta = 0;
		if (phil>=0 && phir>=0) {
			theta = 1.0;
		} else if (phil<0 && phir<0) {
			theta = 0.0;
		} else {
			double lp = std::max(phil, 0.0);
			double rp = std::max(phir, 0.0);
			double la = abs(phil);
			double ra = abs(phir);
			theta = (lp+rp) / (la+ra);
		}
		return theta;
	}
};
#ifdef TESTDIR
#undef TESTDIR
#endif


