/*
 * structure_alignment.h
 *
 *  Created on: Jan 25, 2022
 *      Author: rene staritzbichler
 */

#ifndef INCLUDE_STRUCTURE_ALIGNMENT_H_
#define INCLUDE_STRUCTURE_ALIGNMENT_H_

#include <vector>
#include <map>

#include "function.t.h"
#include "amino_acid.h"

#include "matrix.t.h"
#include "dynamic_programing_matrix_element.h"
#include "triplet.t.h"

#include "sequence.h"
#include "alignment_variables.h"

#include <Eigen/Dense>

// How to get domain alignment? potentially via filter when reading

struct Rotation
{
	Eigen::Vector3d center1;
	Eigen::Matrix3d rotation;
	Eigen::Vector3d center2;

	Rotation( const Eigen::Vector3d &C1, const Eigen::Matrix3d &ROT, const Eigen::Vector3d &C2);

	Eigen::Vector3d
	Position( const Eigen::Vector3d &POS) const;

	std::vector< Eigen::Vector3d >
	Position( const std::vector< Eigen::Vector3d >  &POS) const;

	// Eigen::MatrixXd ??
};

class Align3DScore
: public Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double>
{
private:
	Matrix< DynamicProgrammingMatrixElement>                                &m_Matrix;
	Matrix<double>                                                           m_RMSD;
	std::vector<Eigen::Vector3d>                                             m_First;
	std::vector<Eigen::Vector3d>                                             m_Second;
public:
	Align3DScore(  Matrix< DynamicProgrammingMatrixElement> & M, const std::vector<Eigen::Vector3d> &FIRST, const std::vector<Eigen::Vector3d> &SECOND);

	virtual double operator() ( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA);
};

std::string&
Remove( std::string &STR, char *TO_REMOVE);

std::map< std::string, char>
CreateAAMap();

char
AA3to1( const std::string AA, const std::map< std::string, char> &MAP);

void
ReadPDB( const std::string &NAME , std::vector< Eigen::Vector3d> &POS, Sequence &SEQ );

void
WriteRotatedPDB( const std::string &NAME , /*const std::vector< Eigen::Vector3d> &POS,*/ const Rotation &ROT, const std::string &OUT );

Rotation
Superimpose( std::vector< Eigen::Vector3d> &A,  std::vector< Eigen::Vector3d> &B);

double
RMSD( const std::vector< Eigen::Vector3d> &A, const std::vector< Eigen::Vector3d> &B);

void BuildStructureAlignmentMatrix( const AlignmentVariables &VARS, Matrix< DynamicProgrammingMatrixElement> &MATRIX);


#endif /* INCLUDE_STRUCTURE_ALIGNMENT_H_ */
