/*
 * structure_alignment.h
 *
 *  Created on: Jan 25, 2022
 *      Author: hildilab
 */

#ifndef INCLUDE_STRUCTURE_ALIGNMENT_H_
#define INCLUDE_STRUCTURE_ALIGNMENT_H_

#include "matrix.t.h"
#include "dynamic_programing_matrix_element.h"
#include "triplet.t.h"

//#include "rmsd.h"


class 3DAlignScore
: public Function
{
private:
	double                                              m_PrevScore;
	int                                                 m_NrPrev;
	ShPtr< Matrix< DynamicProgramingMatrixElement> >    m_DPM;
	vec<V3>                                             m_First;
	vec<V3>                                             m_Second;
public:
	3DAlignScore(){};
};


void BuildStructureAlignmentMatrix( const std::vector<std::string> &FILES, Matrix< DynamicProgrammingMatrixElement> &MATRIX)
{
	// read pdbs
	// - sequence from ATOM section
	// - pos into scorefct? (GAA, DPME)
	// - filter: CA atoms! 1 per residue  || user defined

	ReadPDB( FILES[0], first_seq, first_pos);
	ReadPDB( FILES[1], second_seq, second_pos);

	3DAlignScore score;

	NeedlemanWunsch ali( ...penalty, first_seq, second_seq, score, MATRIX);

	ali.CalculateMatrix();
	auto alignment = ali.TraceBack();

	// score function:
	// - store: vec<V3> first, second
	// - vec<double> || double,int prev score and nr of aligned atoms, length equals diagonals(aligned) bits, gap set to zero
	// - ptr to matrix to go back
	// - operator()(GAA1, GAA2)    (normalize!)
	//    - find out which path best scores comes from
	// - calc rmsd
	// needleman wunsch: calc matrix()
	// - construct: Sequences, gaps, scorefct, ptr to matrix

}


#endif /* INCLUDE_STRUCTURE_ALIGNMENT_H_ */
