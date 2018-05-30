#pragma once

#include <cstdlib>

#include "SayakaTree.h"
#include "SayakaTreeData.h"

SAYAKA_NS_BEGIN;


/**
 * Interpolate between cell/face/node data.
 * This only operates block-by-block on the same tree level.
 * It is required that boundary cells have been properly populated before this.
 */
struct VarLocInterpolator 
{
	const AmrTree *m_tree;

	VarLocInterpolator() : m_tree(nullptr) {}

	VarLocInterpolator(const AmrTree &tree)
		: m_tree(&tree)
	{}

	virtual ~VarLocInterpolator() {}

public:

	/*
	 * Low-level, block-based
	 */

	//
	void interp_block_cell_to_face(int dir,
		DoubleBlockData &face_data,
		const DoubleBlockData &cell_data, 
		const IndexBox &face_box,
		int dstcomp, int srccomp) const;
	//
	void interp_block_face_to_cell(int dir,
		DoubleBlockData &cell_data,
		const DoubleBlockData &face_data,
		const IndexBox &cell_box,
		int dstcomp, int srccomp) const;
	//
	void interp_block_cell_to_node(
		DoubleBlockData &node_data,
		const DoubleBlockData &cell_data,
		const IndexBox &node_box,
		int dstcomp, int srccomp) const;
	//
	void interp_block_node_to_cell(
		DoubleBlockData &cell_data,
		const DoubleBlockData &node_data,
		const IndexBox &cell_box,
		int dstcomp, int srccomp) const;
	//
	void interp_block_face_to_node(int dir,
		DoubleBlockData &node_data,
		const DoubleBlockData &face_data,
		const IndexBox &node_box,
		int dstcomp, int srccomp) const;

	/*
	 * Work on all tree data
	 */

	//
	void interp_cell_to_face(int dir,
		TreeData &face_data, const TreeData &cell_data,
		int dstcomp, int srccomp, int ncomp, int ngrow) const;
	//
	void interp_face_to_cell(int dir,
		TreeData &cell_data, const TreeData &face_data,
		int dstcomp, int srccomp, int ncomp, int ngrow) const;
	//
	void interp_cell_to_node(
		TreeData &node_data, const TreeData &cell_data,
		int dstcomp, int srccomp, int ncomp, int ngrow) const;

	//
	void interp_face_to_node(int dir,
		TreeData &node_data, const TreeData &face_data,
		int dstcomp, int srccomp, int ncomp, int ngrow) const;

}; // struct_varlocinterpolator






SAYAKA_NS_END;


