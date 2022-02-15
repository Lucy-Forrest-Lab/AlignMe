/*
 * test.cpp
 *
 *  Created on: 29.01.2022
 *      Author: rene
 */


#include "structure_alignment.h"
#include "sequence.h"

int main( int argc, char * argv[])
{
	std::vector< Eigen::Vector3d> pos1, pos2;
	Sequence seq1, seq2;

	ReadPDB( std::string( argv[1]), pos1, seq1);
	ReadPDB( std::string( argv[2]), pos2, seq2);

	std::cout << pos1[0].transpose() << std::endl;
	std::cout << pos1[1].transpose() << std::endl;
	std::cout << pos1[2].transpose() << std::endl << std::endl;
	std::cout << pos2[0].transpose() << std::endl;
	std::cout << pos2[1].transpose() << std::endl;
	std::cout << pos2[2].transpose() << std::endl << std::endl;

	Rotation
		rot = Superimpose( pos1, pos2);

	std::cout << pos1[0].transpose() << std::endl;
	std::cout << pos1[1].transpose() << std::endl;
	std::cout << pos1[2].transpose() << std::endl << std::endl;
	std::cout << pos2[0].transpose() << std::endl;
	std::cout << pos2[1].transpose() << std::endl;
	std::cout << pos2[2].transpose() << std::endl << std::endl;

	WritePDB( argv[1], pos1, rot, argv[3]);
//	std::cout << "main rot: " << rot.rotation << std::endl;

	std::vector< Eigen::Vector3d> pos;
	Eigen::Vector3d v;
	for( std::vector< Eigen::Vector3d>::const_iterator itr = pos1.begin(); itr != pos1.end(); ++itr)
	{
		v = *itr - rot.center1;
		v = rot.rotation * v;
		v += rot.center2;
		pos.push_back( v);
	}

	std::cout << "rmsd: " << RMSD( pos, pos2) << std::endl;
}

