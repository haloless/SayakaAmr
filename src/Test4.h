#pragma once


#include <cmath>

#include <iostream>

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


using namespace sayaka;



#define TESTDIR(path) ("test4" path)
class Test4
{
	enum {
		PhysBC_Fixed = -100,
		PhysBC_Symmetry,
	};

	enum CellComp {
		Comp_Phi = 0,
		Comp_Divu,
		Comp_Resid,
		NumComp,
	};
	enum FaceComp {
		FaceCompVel = 0,
		FaceCompAdv,
		NumFaceComp,
	};
	enum NodeComp {
		NodeCompP,
		NodeCompU,
		NumNodeComp = NodeCompU + NDIM,
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

		probmin = Vector3d::VecMake(0.0, 0.0, 0.0);
		probmax = Vector3d::VecMake(1.0, 1.0, 0.0);
		probbox = RealBox(probmin, probmax);

		// problem

		//
		//n0x = 6;
		//n0y = 6;
		//n0x = 4;
		//n0y = 4;
		n0x = 1;
		n0y = 1;
		n0z = 1;

		//minlevel = 1;
		//maxlevel = 3;
		minlevel = 3;
		maxlevel = 5;
		//maxlevel = 7;

		nb = 4;
		nbx = nb;
		nby = nb;
		nbz = 1;
		ng = 2;
		ngx = ng;
		ngy = ng;
		ngz = 0;

		validBox = IndexBox(
			Vector3i::VecMake(0,0,0),
			Vector3i::VecMake(nbx-1,nby-1,nbz-1));
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

		this->move_umac_to_unode();
		// save data
		m_tree->writeTreeLeafBlockGrid(TESTDIR("/amr_leaf_grid_init.tec.dat"), m_step+1, m_time);
		m_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_data_init.tec.dat"), m_step+1, m_time);
		m_node_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_nodedata_init.tec.dat"), m_step+1, m_time);
		

		// dynamic refinement
		for (int level=minlevel; level<maxlevel+10; level++) {
			for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
				if ((*m_tree)[iblock].getLevel() == level) {
				this->test_block_refine(iblock);
				}
			}

			m_tree->regrid();

			//
			m_tree->initNewBlockDataAfterRegrid();

			//m_fillpatch->restrictAll(Comp_Phi, NumComp);
			for (int dir=0; dir<NDIM; dir++) {
				m_face_fillpatch[dir]->restrictAll(0, NumFaceComp);
			}
		}

		if (1) {
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
		
		this->move_umac_to_unode();
		// save data
		m_tree->writeTreeLeafBlockGrid(TESTDIR("/amr_leaf_grid.tec.dat"), m_step+1, m_time);
		m_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_data.tec.dat"), m_step+1, m_time);
		m_node_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_nodedata.tec.dat"), m_step+1, m_time);

		for (int level=1; level<=maxlevel; level++) {
			char filename[64];
			sprintf(filename, TESTDIR("/amr_level%02d_nodedata.tec.dat"),level);
			m_node_data->writeTreeLevelBlockData(level, filename, m_step+1, m_time);
		}

		if (1) {
			std::auto_ptr<MacSolver> mac_op(new MacSolver(*m_tree));
			mac_op->setVerbose(1);

			mac_op->prepareSolver();

			//
			//mac_op->setScalars(0.0, 0.25);
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

			if (1) { // Neumann BC, shift the solution
				double sumsol = 0;
				double volume = 0;
				for (int i=0; i<m_tree->numBlocks; i++) {
					if ((*m_tree)[i].isLeaf()) {
						double sum = DoubleGridDataUtil.ReduceSum(
							(*m_data)[i], m_tree->validBlockCellBox(), Comp_Phi);

						double dv = m_tree->getBlockCellVolume(i);
						
						sumsol += sum * dv;
						volume += dv * m_tree->validBlockCellBox().stride();
					}
				}

				double avgsol = sumsol / volume;
				for (int i=0; i<m_tree->numBlocks; i++) {
					if ((*m_tree)[i].isLeaf()) {
						DoubleGridDataUtil.AddEqual(
							(*m_data)[i], Comp_Phi, 
							m_tree->validBlockCellBox(), 
							-avgsol);
					}
				}
			}

			// now the data only lives on leaf nodes
			// have to average down to parents
			m_fillpatch->restrictAll(Comp_Phi, 1);
			move_pcell_to_pnode();

			if (1) { // flux 
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
			//m_tree->writeTreeLeafBlockGrid(TESTDIR("/amr_leaf_grid.tec.dat"), m_step+1, m_time);
			m_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_data_sol.tec.dat"), m_step+1, m_time);
			m_node_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_nodedata_sol.tec.dat"), m_step+1, m_time);
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

		{ // setup root level
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
						block_lo.y += probbox.length(1) / n0y * j;
						block_hi.x += probbox.length(0) / n0x * (i+1);
						block_hi.y += probbox.length(1) / n0y * (j+1);

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
									ineigh = (n0x-1) + j*n0x + k*n0x*n0y;
								} else {
									ineigh = (i-1) + j*n0x + k*n0x*n0y;
								}
							} else if (iface == 1) {
								if (i == n0x-1) {
									ineigh = 0 + j*n0x + k*n0x*n0y;
								} else {
									ineigh = (i+1) + j*n0x + k*n0x*n0y;
								}
							} else if (iface == 2) {
								if (j == 0) {
									ineigh = i + (n0y-1)*n0x + k*n0x*n0y;
								} else {
									ineigh = i + (j-1)*n0x + k*n0x*n0y;
								}
							} else if (iface == 3) {
								if (j == n0y-1) {
									ineigh = i + (0)*n0x + k*n0x*n0y;
								} else {
									ineigh = i + (j+1)*n0x + k*n0x*n0y;
								}
							}
							assert(0<=ineigh && ineigh<n0x*n0y*n0z);
							root.neighbor[iface] = ineigh;
						} // end loop face
					}
				}
			}
		}

		// static refine up to min level
		m_tree->initRefineToMinLevel();
	}

	// must be called after tree initialization
	void init_data() {
		// cell data
		m_data = new TreeData(*m_tree, validBox, grownBox, ng, NumComp);
		m_data->setValue(0);
		
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

				if (1) {
					bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Fixed, Dirichlet);
				} else {
					bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Fixed, SimpleFill);
				}

				bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Symmetry, Neumann);
			}
		}

		m_fillpatch = new FillPatch(*m_tree, *m_data/*, VARLOC_CELL*/);
		for (int comp=0; comp<NumComp; comp++) {
			int prolong = InterpPatch::PROLONG_CENTERED_LINEAR;
			int restrict = InterpPatch::RESTRICT_AVERAGE;

			if (comp == Comp_Divu) {
				prolong = InterpPatch::PROLONG_IGNORE;
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
		for (int dir=0; dir<NDIM; dir++) {
			TreeStateData face_state;
			face_state.new_data = m_face_data[dir];
			face_state.bc_patch = NULL;
			face_state.fill_new_patch = m_face_fillpatch[dir];

			m_tree->treeStateData.push_back(face_state);
		}

	}

	void test_block_refine(int iblock) {
		AmrTreeNode &block = (*m_tree)[iblock];
		int level = block.getLevel();

		if (block.isLeaf() || block.isParentOfLeafChild()) {
			const Vector3d &blockCenter = (*m_tree)[iblock].blockCenter;
			const Vector3d &blockLength = (*m_tree)[iblock].blockLength;

			const Vector3d spot0 = Vector3d::VecMake(0.25, 0.25, 0.0);
			const Vector3d spot1 = Vector3d::VecMake(0.75, 0.75, 0.0);
			//const Vector3d spot2 = Vector3d::VecMake(0.5, 0.5, 0.0);


			int needRefine = 0;
			int needCoarsen = 0;

			if (Vec3Length(blockCenter, spot0) < blockLength.x) {
				needRefine = 1;
			} else if (Vec3Length(blockCenter, spot1) < blockLength.x) {
				needRefine = 1;
			}
			//else if (Vec3Length(blockCenter, spot2) < blockLength.x) {
			//	needRefine = 1;
			//}

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

			for (int k=vbox.klo(); k<=vbox.khi(); k++) {
				for (int j=vbox.jlo(); j<=vbox.jhi(); j++) {
					for (int i=vbox.ilo(); i<=vbox.ihi(); i++) {
						divu(i,j,k,comp) += 1.0/dh(dir) * 
							(umac(i+ii,j+jj,k+kk,maccomp)-umac(i,j,k,maccomp));
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

			for (int k=face_box.klo(); k<=face_box.khi(); k++) {
				double z = zlo + ((double) k + 0.5*(1-kk)) * dz;
				for (int j=face_box.jlo(); j<=face_box.jhi(); j++) {
					double y = ylo + ((double) j + 0.5*(1-jj)) * dy;
					for (int i=face_box.ilo(); i<=face_box.ihi(); i++) {
						double x = xlo + ((double) i + 0.5*(1-ii)) * dx;

						for (int comp=0; comp<NumFaceComp; comp++) {
							double val = calc_mac_value(dir, comp, x,y,z);
							block_data(i,j,k,comp) = val;
						}
					}
				}
			}
		}
	}

	double calc_mac_value(int dir, int comp, double x, double y, double z) const {
		const double k = 2.0 * M_PI;
		
		const double skx = sin(k*x);
		const double ckx = cos(k*x);
		const double sky = sin(k*y);
		const double cky = cos(k*y);

		const double u = skx * cky;
		const double v = -ckx * sky;

		double val = 0;
		if (comp == FaceCompVel) {
			if (dir == 0) {
				val = u;
			} else { 
				val = v;
			}
		} else if (comp == FaceCompAdv) {
			if (dir == 0) {
				const double dudx = k * ckx * cky;
				const double dudy = -k * skx * sky;
				val = u*dudx + v*dudy;
			} else {
				const double dvdx = k * skx * sky;
				const double dvdy = -k * ckx * cky;
				val = u*dvdx + v*dvdy;
			}
		}
		return val;
	}

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
	}
};
#ifdef TESTDIR
#undef TESTDIR
#endif

