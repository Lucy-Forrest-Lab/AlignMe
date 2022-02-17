/*
 * structure_alignment.h
 *
 *  Created on: Jan 25, 2022
 *      Author: rene staritzbichler
 */


#include "structure_alignment.h"

#include "needleman_wunsch.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <boost/format.hpp>

#include "alignment_write.h"

/*
import numpy as np

def PCA(X , num_components):

    #Step-1
    X_meaned = X - np.mean(X , axis = 0)

    #Step-2
    cov_mat = np.cov(X_meaned , rowvar = False)

    #Step-3
    eigen_values , eigen_vectors = np.linalg.eigh(cov_mat)

    #Step-4
    sorted_index = np.argsort(eigen_values)[::-1]
    sorted_eigenvalue = eigen_values[sorted_index]
    sorted_eigenvectors = eigen_vectors[:,sorted_index]

    #Step-5
    eigenvector_subset = sorted_eigenvectors[:,0:num_components]

    #Step-6
    X_reduced = np.dot(eigenvector_subset.transpose() , X_meaned.transpose() ).transpose()

    return X_reduced
*/



Rotation::Rotation( const Eigen::Vector3d &C1, const Eigen::Matrix3d &ROT, const Eigen::Vector3d &C2)
: center1( C1), rotation(ROT), center2( C2)
{}



Eigen::Vector3d Rotation::Position( const Eigen::Vector3d &POS) const
{
	Eigen::Vector3d
		v = POS - center1;
	v = rotation * v;
	v += center2;
	return v;
}



std::vector< Eigen::Vector3d >
Rotation::Position( const std::vector< Eigen::Vector3d >  &POS) const
{
	std::vector< Eigen::Vector3d >
		v( POS.size());
	std::vector< Eigen::Vector3d >::const_iterator ptr = POS.begin();
	for( std::vector< Eigen::Vector3d >::iterator itr = v.begin(); itr != v.end(); ++itr, ++ptr)
	{
		*itr = Position( *ptr);
	}
	return v;
}



Align3DScore::Align3DScore(  Matrix< DynamicProgrammingMatrixElement> & M, const std::vector<Eigen::Vector3d> &FIRST, const std::vector<Eigen::Vector3d> &SECOND)
: m_Matrix( M ),
  m_RMSD( M.GetNumberOfRows(), M.GetNumberOfColumns()),
  m_First(FIRST),
  m_Second(SECOND)
  {}



double Align3DScore::operator() ( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA)
{
	///  ***********************************************************************
	double base = 10;            ///    MAKE MEMBER !!! adjust to gap penalties
	double default_return = 0.0; ///    MAKE MEMBER !!! use GAP PENALTY ????
	double default_rmsd = 0;    ///    MAKE MEMBER !!! use GAP PENALTY ????
	///  ***********************************************************************

	std::cout << __FUNCTION__ << std::endl;
	int
		i = AA.first.GetSeqID(),
		j = AA.second.GetSeqID(),
		ip = i - 1,
		jp = j - 1;
	std::vector< Eigen::Vector3d>
		pos1,
		pos2;
	if( i < 3 || j < 3 )
//		if( i < 2 || j < 2 )
	{
		std::cout << "index: " << i << " " << j << std::endl;
		m_RMSD(i,j) = default_rmsd;
		return default_return;
	}

	std::cout << "push back 1st: " << i << " " << ip << " " << m_First.size() << std::endl;
	pos1.push_back( m_First[i]);
	pos1.push_back( m_First[ip--]); // ??
	pos1.push_back( m_First[ip]);

	std::cout << "push back 2nd: " << j << " " << jp << " " << m_Second.size() << std::endl;
	pos2.push_back( m_Second[j]);
	pos2.push_back( m_Second[jp--]);
	pos2.push_back( m_Second[jp]);

	std::cout << "now: " << ip << " " << jp << std::endl;
	int counter = 0;
	std::cout << "loop" << std::endl;
	do{
		std::pair< size_t, size_t>
			prev = m_Matrix(ip,jp).GetIndicesOfPreviousElement();
		std::cout << "prev: " << prev.first << " " << prev.second << std::endl;
		if( ip - prev.first == 1 && jp - prev.second == 1){                            ///  **********************************
			std::cout << "add" << std::endl;
			pos1.push_back( m_First[prev.first]);
			pos2.push_back( m_Second[prev.second]);
		}
		else{
			if( counter++ > 1){
				break;}           	             ///  kills backward moving once 2 gaps are introduced  ??
		}
		ip = prev.first;
		jp = prev.second;
	} while( ip > 0 && jp > 0);


	if( pos1.size() != pos2.size() || pos1.size() < 3)
	{
		std::cout << "WARNING: size issue in " << __FUNCTION__ << " " << pos1.size() << " " << pos2.size() << std::endl;
		m_RMSD(i,j) = default_rmsd;
		return default_return;
	}
	Rotation rot = Superimpose( pos1, pos2);
	pos1 = rot.Position( pos1);
	double rmsd = RMSD( pos1, pos2);
	std::cout << __FUNCTION__ << " len: " << pos1.size() << " rmsd: " << rmsd << " ids: " << i << " " << j << std::endl;
	m_RMSD(i,j) = rmsd;                    // store actual RMSD in matrix
	std::cout << "assigned, max: " << m_RMSD.GetNumberOfRows() << " " << m_RMSD.GetNumberOfColumns() << std::endl;
//		return m_RMSD(i-1,j-1) - rmsd;         // return diff with prev element as score  **********************************
	return base * exp(-0.69315 * rmsd);         // lambda = ln2 / r1/2, ln2 = 0.69315, r1/2 = 1   **********************************
}



std::map< std::string, char>
CreateAAMap()
{
	std::map< std::string, char> aa;
	aa.insert( std::make_pair( "CYS", 'C' ) );
	aa.insert( std::make_pair( "ILE", 'I' ) );
	aa.insert( std::make_pair( "SER", 'S' ) );
	aa.insert( std::make_pair( "GLN", 'Q' ) );
	aa.insert( std::make_pair( "MET", 'M' ) );
	aa.insert( std::make_pair( "ASN", 'N' ) );
	aa.insert( std::make_pair( "PRO", 'P' ) );
	aa.insert( std::make_pair( "LYS", 'K' ) );
	aa.insert( std::make_pair( "ASP", 'D' ) );
	aa.insert( std::make_pair( "THR", 'T' ) );
	aa.insert( std::make_pair( "PHE", 'F' ) );
	aa.insert( std::make_pair( "ALA", 'A' ) );
	aa.insert( std::make_pair( "GLY", 'G' ) );
	aa.insert( std::make_pair( "HIS", 'H' ) );
	aa.insert( std::make_pair( "LEU", 'L' ) );
	aa.insert( std::make_pair( "ARG", 'R' ) );
	aa.insert( std::make_pair( "TRP", 'W' ) );
	aa.insert( std::make_pair( "VAL", 'V' ) );
	aa.insert( std::make_pair( "GLU", 'E' ) );
	aa.insert( std::make_pair( "TYR", 'Y' ) );
	return aa;
}



char
SaveGet( const std::string AA, const std::map< std::string, char> &MAP)
{
	std::map< std::string, char>::const_iterator atr = MAP.find( AA);
	if( atr == MAP.end())
	{
		return 'x';
	}
	return atr->second;
}



// local function ??
std::string
Erase( std::string STR, const std::string &TO_REMOVE)
{
	 for ( unsigned int i = 0; i < TO_REMOVE.size(); ++i ) {
	      STR.erase( remove(STR.begin(), STR.end(), TO_REMOVE[i]), STR.end() );
	   }
	return STR;
}



void
ReadPDB( const std::string &NAME , std::vector< Eigen::Vector3d> &POS, Sequence &SEQ )
{
	std::cout << __FUNCTION__ <<  " " << NAME << std::endl;
	std::map< std::string, char> aa = CreateAAMap();
	std::ifstream read( NAME);
	if( !read){ std::cerr << "ERROR: file not found: " << NAME << std::endl; exit(1);}
	std::string
		line;
	int
		prev_id = -999999, id;
	char
		prev_chain = '%', chain;

	int count = 1;
	while( std::getline( read, line))
	{
		if( line.substr( 0, 4) == "ATOM" && Erase( line.substr( 12, 4 ), " ") == "CA")
		{
			chain = line[21];
			id = std::stoi( line.substr(22,4));
			if( id != prev_id || prev_chain != chain)
			{
				SEQ.push_back( GeneralizedAminoAcid( SaveGet( line.substr( 17,3), aa ), 0, count++ ));
				POS.push_back( Eigen::Vector3d( std::stod( line.substr(30,8)), std::stod( line.substr( 38, 8)), std::stod( line.substr( 46, 8))));
				prev_chain = chain;
				prev_id = id;
			}
		}
	}
	std::cout << __FUNCTION__ << " #lines: " << POS.size() << std::endl;
	read.close();
	read.clear();
}



void
WriteRotatedPDB( const std::string &NAME /*, const std::vector< Eigen::Vector3d> &POS*/, const Rotation &ROT, const std::string &OUT )
{
	std::cout << __FUNCTION__ << std::endl;
	std::ifstream in( NAME);
	std::ofstream out( OUT);

	if( !in){ std::cerr << "ERROR: file not found: " << NAME << std::endl; exit(1);}
	if( !out){ std::cerr << "ERROR: file not found: " << OUT << std::endl; exit(1);}
	std::string
		line;

	int
		prev_id = -999999, id;
	char
		prev_chain = '%', chain;
	Eigen::Vector3d
		pos;
	std::cout << "rot: " << ROT.rotation << std::endl<< std::endl;

	int count = 0;
	while( std::getline( in, line))
	{
		if( line.substr( 0, 4) == "ATOM" || line.substr(0,6) == "HETATM")
		{
				pos = Eigen::Vector3d( std::stod( line.substr(30,8)), std::stod( line.substr( 38, 8)), std::stod( line.substr( 46, 8)));
				pos = ROT.Position( pos);
//				if( count++ < 3){
//					std::cout << pos.transpose() << std::endl;}
				out << line.substr(0, 30);
				out << boost::format("%1$8.3f%2$8.3f%3$8.3f") % pos[0] % pos[1] % pos[2];
				out << line.substr( 54, line.size() - 54 ) << std::endl;
		}
		else
		{
			out << line;
		}
	}
//	std::cout << __FUNCTION__ << " #lines: " << POS.size() << std::endl;
	in.close();
	in.clear();
	out.close();
	out.clear();
}



Eigen::Vector3d Center( const std::vector< Eigen::Vector3d> &P)
{
	std::cout << __FUNCTION__ << std::endl;

	Eigen::Vector3d sum( 0.0, 0.0, 0.0);
	for( std::vector< Eigen::Vector3d>::const_iterator itr = P.begin(); itr != P.end(); ++itr)
	{
		sum += *itr;
	}
	return sum / (double) P.size();
}



Eigen::Vector3d Center( const Eigen::MatrixXd &P)
{
	std::cout << __FUNCTION__ << " sizes: " << P.rows() << " " << P.cols() << std::endl;

	Eigen::Vector3d sum( 0.0, 0.0, 0.0);
	for( int i = 0; i < P.rows(); ++i)
	{
		sum += P.row(i);
	}
	return sum / (double) P.rows();
}



void OrthogonalizeRows( Eigen::Matrix3d &M, int ROW = 2)
{
	std::cout << __FUNCTION__ << " row: " << ROW << std::endl;

	M.row( ROW % 3) = M.row( (ROW+1)%3 ).cross( M.row( (ROW+2)%3 ) );
}



void SwapRows( Eigen::Matrix3d &M, int R1, int R2)
{
	std::cout << __FUNCTION__ << std::endl;
	Eigen::Vector3d
		m1 = M.row(R1),
		m2 = M.row(R2);
	M.row(R1) = m2;
	M.row(R2) = m1;
}



Rotation
Superimpose(  std::vector< Eigen::Vector3d> &A,  std::vector< Eigen::Vector3d> &B)
{
	std::cout << __FUNCTION__ << std::endl;

	Eigen::Vector3d
		c1 = Center( A),
		c2 = Center( B);

//	Eigen::MatrixXd
//		m1 = Eigen::Map<Eigen::MatrixXd>(A[0].data(), 3, A.size()),
//		m2 = Eigen::Map<Eigen::MatrixXd>(B[0].data(), 3, B.size()),
//		mol1 = m1.transpose(),
//		mol2 = m2.transpose();
	Eigen::MatrixXd
		mol1 = Eigen::Map<Eigen::MatrixXd>(A[0].data(), A.size(), 3),
		mol2 = Eigen::Map<Eigen::MatrixXd>(B[0].data(), B.size(), 3);

	for( int i = 0; i < A.size(); ++i)
	{
		mol1.row(i) = A[i];
		mol2.row(i) = B[i];
	}

	for( int i = 0; i < mol1.rows(); ++i)
	{
		mol1.row(i) -= c1.transpose();
		mol2.row(i) -= c2.transpose();
	}

	Eigen::MatrixXd
		moment = mol1.transpose() * mol2;

	// symmetric
	Eigen::MatrixXd
		rotate = moment * moment.transpose();
//	std::cout << "rotate sizes: " << rotate.rows() << " " << rotate.cols() << std::endl;
//	std::cout << "rotate: \n" << rotate << std::endl;

//	std::cout << "solve "<< std::endl;
	// solve eigensystem
	Eigen::EigenSolver<Eigen::MatrixXd> eigensolver( rotate );

//	std::cout << "eigenvalues: " << eigensolver.eigenvalues().transpose() << std::endl;
//	std::cout << "eigenvectors: " << eigensolver.eigenvectors().real().transpose() << std::endl;

	Eigen::VectorXd
		eigenvalues = eigensolver.eigenvalues().real();
	Eigen::Matrix3d
		eigenvectors = eigensolver.eigenvectors().real().transpose();
//	std::cout << "determinant eigenvectors: " << eigenvectors.determinant() << std::endl;

    if( eigenvalues(0) < eigenvalues(1)){
//        std::cout << "prior swap 1: " << eigenvalues.transpose() << std::endl;
//        std::cout << "prior swap 1: " << eigenvectors << std::endl;
        SwapRows( eigenvectors, 0,1);
    	double e0 = eigenvalues(0);
    	eigenvalues(0) = eigenvalues(1);
    	eigenvalues(1) = e0;
//        std::cout << "post swap 1: " << eigenvalues.transpose() << std::endl;
//        std::cout << "post swap 1: " << eigenvectors << std::endl;
    }
    if( eigenvalues(1) < eigenvalues(2)){
//        std::cout << "prior swap 2: " << eigenvalues.transpose() << std::endl;
//        std::cout << "prior swap 2: " << eigenvectors << std::endl;
        SwapRows( eigenvectors, 1,2);
    	double e1 = eigenvalues(1);
    	eigenvalues(1) = eigenvalues(2);
    	eigenvalues(2) = e1;
//        std::cout << "post swap 2: " << eigenvalues.transpose() << std::endl;
//        std::cout << "post swap 2: " << eigenvectors << std::endl;
    }
    if( eigenvalues(0) < eigenvalues(1)){
//        std::cout << "prior swap 3: " << eigenvalues.transpose() << std::endl;
//        std::cout << "prior swap 3: " << eigenvectors << std::endl;
        SwapRows( eigenvectors, 0,1);
    	double e0 = eigenvalues(0);
    	eigenvalues(0) = eigenvalues(1);
    	eigenvalues(1) = e0;
//        std::cout << "post swap 3: " << eigenvalues.transpose() << std::endl;
//        std::cout << "post swap 3: " << eigenvectors << std::endl;
    }

    if( eigenvalues(0) <= 0.0 || eigenvalues(1) <= 0.0)
    {
    	std::cerr << "ERROR in eigenvalues" << std::endl;
    	exit(1);
    }

    OrthogonalizeRows( eigenvectors);
    Eigen::Matrix3d
		rot = eigenvectors * moment;

    for( int i = 0; i < 2; ++i){
    	double
			t = 1.0 / sqrt( eigenvalues(i));
//    	std::cout << "t: " << t << std::endl;
        for( int j = 0; j < 3; ++j){
            rot(i,j) *= t;  // todo: more efficient ways?
        }
    }
//    std::cout << "rot rep: " << std::endl;
//    std::cout << rot << std::endl;
    OrthogonalizeRows( rot );
//    std::cout << "build rotation matrix" << std::endl;
    rot = rot.transpose() * eigenvectors;
//    std::cout << "rotation: " << rot << std::endl;
//    std::cout << "determinant: " << rot.determinant() << std::endl;
//    std::cout << "rotate coors " << std::endl;
    return Rotation( c1, rot, c2);
}



double
RMSD( const std::vector< Eigen::Vector3d> &A, const std::vector< Eigen::Vector3d> &B)
{
	std::cout << __FUNCTION__ << std::endl;
    // check size of vectors
    if( A.size() != B.size())
    {
    	std::cout << "WARNING: vectors do match in " << __FUNCTION__ << std::endl;
    	return std::numeric_limits<double>::max();
    }
    // compute RMSD
    size_t number = A.size();
    double sumsquarermsd( 0);
    for( size_t i = 0; i < number; i++)
      sumsquarermsd += ( A[i] - B[i] ).squaredNorm();
    return sqrt( sumsquarermsd / number);
}



void BuildStructureAlignmentMatrix( const AlignmentVariables &VARS, Matrix< DynamicProgrammingMatrixElement> &MATRIX)
{
	std::cout << __FUNCTION__ << std::endl;
	// read pdbs
	// - sequence from ATOM section
	// - pos into scorefct? (GAA, DPME)
	// - filter: CA atoms! 1 per residue  || user defined

	std::vector< Eigen::Vector3d>
		first_pos ,
		second_pos;
	Sequence
		first_seq,
		second_seq;
	ReadPDB( VARS.first_pdb_file, first_pos, first_seq);
	ReadPDB( VARS.second_pdb_file, second_pos, second_seq);
	std::cout << "first 0: " << first_pos[0].transpose() << std::endl;
	std::cout << "second 0: " << second_pos[0].transpose() << std::endl;
	std::cout << "lengths: " << first_pos.size() << " " << first_seq.size() << " , " << second_pos.size() << " " << second_seq.size() << std::endl;

	std::cout << "score" << std::endl;
	ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
		score( new Align3DScore( MATRIX, first_pos, second_pos));
	std::cout << "needle" << std::endl;
	NeedlemanWunsch ali(
			VARS.gap_opening_penalty,
			VARS.gap_extension_penalty,
			VARS.termini_gap_opening_penalty,
			VARS.termini_gap_extension_penalty,
			first_seq,
			second_seq,
			score,
			MATRIX);

	std::cout << "calc matrix" << std::endl;
	ali.CalculateMatrix();
	std::cout << "trace back" << std::endl;
	std::pair< double, std::vector< std::pair< int, int> > >
		alignment = ali.TraceBack();

	std::cout << "dyn prog matrix:" << std::endl;
	MATRIX.Write( std::cout );
	std::cout << std::endl;


	std::cout << "final score: " << alignment.first << std::endl;

	/// TODO MOVE THIS INTO BACK OF PROGRAMM TO INCLUDE SEQ ALIG !!!!  ****************************************************
	std::cout << "alignment:" << std::endl;
	std::vector< Eigen::Vector3d>
		mol1, mol2;
	for( std::vector< std::pair< int, int> >::const_iterator itr = alignment.second.begin(); itr != alignment.second.end(); ++itr)
	{
		std::cout << itr->first << " :: " << itr->second << std::endl;
		if( itr->first < first_pos.size() && itr->second < second_pos.size())
		{
			mol1.push_back( first_pos[itr->first]);
			mol2.push_back( second_pos[itr->second]);
		}
	}
	Rotation rot = Superimpose( mol1, mol2);
	WriteRotatedPDB( VARS.first_pdb_file, rot,  "aligned.pdb");
	mol1 = rot.Position( mol1);
	std::cout << "rmsd: " << RMSD( mol1, mol2) << std::endl;
	WriteAlignedSequencesInClustalwFormat( alignment, first_seq, second_seq, ">first", ">second", "aligned3D.fa", 60, VARS.gap_extension_penalty, VARS.anchors);

	///// **************************************88888

	// score function:
	// - store: std::vector<Eigen::Vector3d> first, second
	// - std::vector<double> || double,int prev score and nr of aligned atoms, length equals diagonals(aligned) bits, gap set to zero
	// - ptr to matrix to go back
	// - operator()(GAA1, GAA2)    (normalize!)
	//    - find out which path best scores comes from
	// - calc rmsd
	// needleman wunsch: calc matrix()
	// - construct: Sequences, gaps, scorefct, ptr to matrix

}

