#pragma once

/*
Refined distance function of a circle.
Show how to search a point in tree blocks.
*/

#include <cmath>

#include <array>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>

#include "log.h"
#include "vector3d.h"

#include "SayakaCommons.h"
#include "Sayaka.h"
#include "SayakaTree.h"
#include "SayakaTreeData.h"
#include "SayakaBoundaryPatch.h"
#include "SayakaFillPatch.h"
#include "SayakaVarLocInterpolator.h"
//#include "SayakaMacSolver.h"
#include "SayakaWriter.h"

using namespace sayaka;

class Test8 {
protected:

	double m_dist_range;

	//----------------------------------------
	Vector3i m_base_num;
	RealBox m_base_dom;
	Vector3d m_base_dh;
	
	IndexBox m_valid_box;

	int m_finest_level;



	//----------------------------------------
	std::unique_ptr<AmrTree> m_tree;

	TreeData m_sdf;

public:

	void init(int argc, char *argv[]) {

		{
			m_dist_range = 0.1;
		}

		{
			m_base_dom.vlo = { -1, -1, -0.25 };
			m_base_dom.vhi = { 1, 1, 0.25 };

			m_base_num = { 4, 4, 1 };

			for (int dir = 0; dir < NDIM; dir++) {
				m_base_dh(dir) = m_base_dom.length(dir) / m_base_num(dir);
			}

			m_valid_box.vlo = { 0, 0, 0 };
			m_valid_box.vhi = { 3, 3, 3 };
		}

		{
			m_tree = std::make_unique<AmrTree>();

			m_tree->minAmrLevel = 1;
			m_tree->maxAmrLevel = 4;

			m_tree->blockValidIndex = m_valid_box;
			m_tree->blockGrownIndex = m_valid_box;
			m_tree->blockNumGrow = 0;

			m_tree->init();

			std::array<int, NUM_FACE> bc_type;
			bc_type.fill(-100);
			std::array<bool, NUM_FACE> is_per;
			is_per.fill(0);

			m_tree->initUniformRootLevel(m_base_num, m_base_dom,
				bc_type.data(), is_per.data());
		}
		{
			m_tree->initRefineToMinLevel();

			for (int ilevel = m_tree->minAmrLevel; ilevel <= m_tree->maxAmrLevel; ilevel++) {
				for (int iblock = 0; iblock < m_tree->numBlocks; iblock++) {
					if ((*m_tree)[iblock].getLevel() == ilevel) {
						test_refine_block(iblock);
					}
				}

				m_tree->regrid();
			}

			m_finest_level = m_tree->currentFinestLevel();
		}

		{
			m_sdf = TreeData::MakeNodeData(*m_tree, 1, 0);
			m_sdf.fillValue(0);

			for (int iblock = 0; iblock < m_tree->numBlocks; iblock++) {
				if (!(*m_tree)[iblock].isLeaf()) continue;
				calc_dist_block(iblock);
			}
		}


		{
			WriterMeta wm;
			wm.guessName(m_sdf);
			WriteDataHdf5(m_sdf, "hoge.h5", wm);
		}

		{
			std::ofstream ofs("hoge.csv");
			ofs << "x,y,z,d" << std::endl;
			const int nx = 32;
			const int ny = 32;
			for (int j = 0; j <= ny; j++) {
				for (int i = 0; i <= nx; i++) {
					double x = -0.8 + 1.6 / nx * i;
					double y = -0.8 + 1.6 / ny * j;
					double z = 0;
					double d = eval_dist({ x,y,z });
					ofs << x << ',' << y << ',' << z << ',' << d << std::endl;
				}
			}
		}


	}


	double eval_dist(const Vector3d &vpos) const {
		if (!m_base_dom.contains(vpos)) {
			throw std::out_of_range("Point out of base range");
		}

		// find base level box
		int ibase = (int)((vpos.x - m_base_dom.vlo.x) / m_base_dh.x);
		int jbase = (int)((vpos.y - m_base_dom.vlo.y) / m_base_dh.y);
		int kbase = (int)((vpos.z - m_base_dom.vlo.z) / m_base_dh.z);
		
		int iblock = ibase + jbase * m_base_num.x + kbase * m_base_num.x*m_base_num.y;
		assert((*m_tree)[iblock].boundBox.contains(vpos));

		return eval_dist_block(vpos, iblock);
	}

protected:

	double eval_dist_block(const Vector3d &vpos, int iblock) const {
		auto &block = (*m_tree)[iblock];
		if (block.isLeaf()) {
			return interp_dist_block(vpos, iblock);
		}
		else {
			for (int ichild = 0; ichild < NUM_CHILD; ichild++) {
				int jblock = block.getChild(ichild);
				if ((*m_tree)[jblock].boundBox.contains(vpos)) {
					return eval_dist_block(vpos, jblock);
				}
			}
			return 99999; // should never reach here
		}
	}

	double interp_dist_block(const Vector3d &vpos, int iblock) const {
		auto &block = (*m_tree)[iblock];
		auto &dist = m_sdf[iblock];

		const Vector3d xlo = block.boundBox.lo();
		const Vector3d dh = m_tree->getBlockCellSize(iblock);

		Vector3d vdif = vpos - xlo;
		int ii = (int)(vdif.x / dh.x);
		int jj = (int)(vdif.y / dh.y);
		int kk = (int)(vdif.z / dh.z);

		double rx = vdif.x / dh.x - ii; assert(0 <= rx && rx <= 1);
		double ry = vdif.y / dh.y - jj; assert(0 <= ry && ry <= 1);
		double rz = vdif.z / dh.z - kk; assert(0 <= rz && rz <= 1);
		double rx1 = 1.0 - rx;
		double ry1 = 1.0 - ry;
		double rz1 = 1.0 - rz;
		//return dist(ii, jj, kk, 0);
		return rx1 * ry1 * rz1 * dist(ii, jj, kk, 0)
			+ rx * ry1 * rz1 *dist(ii + 1, jj, kk, 0)
			+ rx1 * ry * rz1 * dist(ii, jj + 1, kk, 0)
			+ rx * ry * rz1 * dist(ii + 1, jj + 1, kk, 0)
			+ rx1 * ry1 * rz * dist(ii, jj, kk+1, 0)
			+ rx * ry1 * rz *dist(ii + 1, jj, kk+1, 0)
			+ rx1 * ry * rz * dist(ii, jj + 1, kk+1, 0)
			+ rx * ry * rz * dist(ii + 1, jj + 1, kk+1, 0)
			;
	}

	double calc_dist(const Vector3d &v) const {
		return hypot(v.x, v.y) - 0.5;
	}

	void calc_dist_block(int iblock) {
		auto &block = (*m_tree)[iblock];
		auto &dist = m_sdf[iblock];

		const Vector3d xlo = block.boundBox.lo();
		const Vector3d dh = m_tree->getBlockCellSize(iblock);

		IndexBox b = m_sdf.validBox();
		BEGIN_FOR_BOX_RANGE(b, i, j, k);
		{
			double xx = xlo.x + dh.x * i;
			double yy = xlo.y + dh.y * j;
			double zz = xlo.z + dh.z * k;
			dist(i, j, k, 0) = calc_dist({ xx,yy,zz });
		}
		END_FOR_BOX_RANGE;

	}

	void test_refine_block(int iblock) {
		auto &block = (*m_tree)[iblock];

		const IndexBox &validbox = m_tree->validBlockCellBox();

		if (block.isLeaf() || block.isParentOfLeafChild()) {
			const Vector3d blo = block.boundBox.lo();
			const Vector3d bhi = block.boundBox.hi();
			//const Vector3d dh = m_tree->getBlockCellSize(iblock);
			//const double dx = dh.x;

			int need_refine = 0;
			int need_coarsen = 0;

			//
			double dist[8];
			for (int kk = 0; kk <= 1; kk++) {
				for (int jj = 0; jj <= 1; jj++) {
					for (int ii = 0; ii <= 1; ii++) {
						double xx = blo.x*(1 - ii) + bhi.x*ii;
						double yy = blo.y*(1 - jj) + bhi.y*jj;
						double zz = blo.z*(1 - kk) + bhi.z*kk;
						double dd = calc_dist({ xx, yy, zz });
						dist[ii + jj * 2 + kk * 4] = dd;
					}
				}
			}

			bool has_pos = std::any_of(dist, dist + 8, [](double d) {return d >= 0; });
			bool has_neg = std::any_of(dist, dist + 8, [](double d) {return d < 0; });
			if (has_pos && has_neg) {
				need_refine = 1;
			}

			const double range = m_dist_range;
			bool is_close = std::any_of(dist, dist + 8,
				[=](double d) {return abs(d) < range; });
			if (is_close) {
				need_refine = 1;
			}

			int ilevel = block.getLevel();
			if (ilevel < m_tree->maxRefineLevel()) {
				if (need_refine) m_tree->refineFlag[iblock] = 1;
			}
			if (ilevel > m_tree->minAmrLevel && !m_tree->refineFlag[iblock]) {
				if (need_coarsen) m_tree->coarsenFlag[iblock] = 1;
			}
		}
	}
};







