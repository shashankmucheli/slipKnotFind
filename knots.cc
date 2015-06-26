// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <interactive/util/filter/knots.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray3P.hh>
#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/ObjexxFCL.hh>
#include <ObjexxFCL/format.hh>

// Numeric Headers
#include <numeric/constants.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyzVector.hh>

// Utility Headers
#include <utility/basic_sys_util.hh>

// C++ Headers
#include <algorithm>   
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>


using namespace ObjexxFCL;
using namespace ObjexxFCL::format;



namespace interactive {
namespace util {
namespace filter {
			
///     Knotfind algorithm:
///
///     knot_filter
///     for each chain in position array
///        run knotfind
///        print if chain is knotted
///     return true if any chain is knotted
///
///
///     knotfind:
///       calculate distances for all c-alphas from res1 to res1+2, sort by
///        shortest distance to simplify smallest triangles first
///     while the chain can be simplified
///       simplify chain
///     if chain has been simplified to a line between the N and C termini
///        return false
///     else
///        try to simplify chain one last time using the area of the triangles
///         to determine which triangles are smallest
///        return true
///
///
///     simplify chain:
///           res1 = first residue of smallest triangle in the chain
///           while res1 != final residue in chain
///              create a plane using the coordinates of res1 and the
///              two residues that follow it in the chain
///           if no pair of adjacent residues in the chain (other than these three)
///            forms a line segment that crosses through this plane
///              form a straight line between res1 and res1+2  (removing res1+1
///               from the chain and all triangles involving it)
///              res1 = first residue of smallest unsimplified triangle in chain
///            else
///              leave these three residues alone
///              res1 = first residue of next smallest triangle in the chain
///

struct Pair {
	int res1;
	int res2;
	float dist;

	//constructor
		explicit
		Pair( int res1_a, int res2_a, float dist_a ):
			res1( res1_a ),
			res2( res2_a ),
			dist( dist_a )
		{}
};

typedef std::vector< Pair > Pairs;
typedef std::vector< Pair >::iterator iterator;

static
bool
compare_Pair(
			 Pair first,
			 Pair second
)
{
	return first.dist < second.dist;
}


	
//forward declarations

static void
insert_into_correct_position(
	Pairs & planeDists,
	float const dist1,
	float const dist2,
	bool done1,
	bool done2,
	int const res1,
	int const res3,
	int const old1st,
	int const old3rd
);

static void
add_to_end_in_order(
  float const dist1,
  float const dist2,
  int const res1,
  int const res3,
  int const old1st,
  int const old3rd,
  Pairs & planeDists
);

static void
calculate_distance(
	FArray2A_float chain, // coordinates of CA-trace
	int const nres,        // logical size of chain
	int const first_res,
	int const second_res,  //only used if by_area is true
	int const third_res,
	float & dist,          // will be changed on output
	bool const by_area     // is this sorting planes by area or by distance
);

static void
find_new_ends_of_plane(
  FArray2A_float chain, // coordinates of CA-trace
	int const nres,        // logical size of chain
	int & res1,
	int & res3,
	bool & res1start,
	bool & res3end,
	int const simplify_start,
	int const simplify_end
);

static void
remove_distances(
	Pairs & planeDists,
	FArray1A_float plane
);

static void
find_closest_planes(
 	FArray2A_float chain,  // coordinates of CA-trace
	int start,
	int end,
	int const nres,         // number of residues in chain
	Pairs & planes,
	bool const by_area      // is this sorting planes by area or by distance
);

static bool
simplify_chain(
	int const nres,
	FArray2A_float chain, // modified on output
	Pairs & planes,
	int const simplify_start,
	int const simplify_end,
	bool const by_area,    //flag for sorting planes by area instead of distance
	int const max_trouble
);

static void
get_plane(
	int & res1, // may be modified
	int & res3, // may be modified
	FArray2A_float chain, // CA trace
	int const nres, // chain length
	FArray1A_float plane // first 4 slots: A, B, C, D
);

static void
find_trouble(
			 FArray1A_float plane,
			 FArray2A_float chain,
			 int const nres, // chain length
			 FArray2A_float trouble_above, // pairs of points describing
			 FArray2A_float trouble_below, // line segments that may intersect the plane
			 float tolerance, // closest allowed approach
			 int & counter, // logical size of trouble_above,trouble_below
			 int max_trouble // physical size of trouble_above,trouble_below
);

static void
is_part_of_plane(
				 int & point1,
				 int & point2,
				 FArray1A_float plane,
				 bool & in_plane
);

static bool
any_cross_plane(
				FArray1A_float plane,
				FArray2A_float chain,
				int const nres, // logical size of chain
				int & ntrouble,
				FArray2A_float t_above,
				FArray2A_float t_below,
				float tolerance, // closest allowed approach
				bool const check_all // if true, trouble vectors and modify t_start,t_end,ntrouble
);

static bool
inside_triangle(
				FArray1A_float plane, // 1-4 = equation of inifinite plane  Ax + By + Cz + D
				FArray2A_float ca_chain,
				int const nres,
				FArray1A_float i_point, // point lying on infinite plane
				float tolerance // closest allowed approach
);

static bool
segment_intersects_triangle(
							FArray1A_float x,
							FArray1A_float y,
							FArray1A_float a,
							FArray1A_float b,
							FArray1A_float c,
							float tolerance // closest allowed approach of segments
);



////////////////////////////////////////////////////////////////////////////////
/// @begin knotfind
///
/// @brief looks for knots in a chain
///
/// @detailed
///
/// @param[in]   chain - in -  c-alpha coordinates; removed residue indicator
/// @param[in]   nres - in - number of residues
/// @param[in]   simplify_start - in - first res of region to simplify
/// @param[in]   simplify_end - in - last res of region to simplify
/// @param[in]   end_trim - in - number of residues to ignore on chain ends
///
/// @return  true if chain contains at least one knot
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors car 9/5/03
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
bool
knotfind(
	FArray2A_float chain, // coordinates of CA-trace, chain(4,i)=1 indicates
	int const nres, // number of residues in chain
	int simplify_start, //can be changed if ends are trimmed
	int simplify_end,   //can be changed if ends are trimmed
	int const end_trim
)
{
	chain.dimension( 4, nres );

	bool knotfind = true; // Return value
	bool by_area = false; //flag to run different version of knotfind using area

// make sure that end_trim*2 is not greater than N-4

  if ((end_trim*2) >= nres-4) {
    std::cout << "END_TRIM set too high, KNOTFIND CANCELLED" << std::endl;
    knotfind = false;
    return knotfind;
  }


// set all 4th rows equal to 0 (denotes "point not removed")
	for ( int i = 1; i <= nres; ++i ) {
		chain(4,i) = 0;
	}
// trim ends
	for ( int i = 1; i <= end_trim; ++i ) {
		chain(4,i) = 1;        //set 4th row equal to 1 (denotes "point removed)
		simplify_start = i+1;  //new start of chain
		chain(4,nres+1-i) = 1;       //set 4th row equal to 1
		simplify_end = nres-i; //new end of chain
	}

	//calculate distances between 1st & 3rd residues in all the planes
	Pairs planes;
	find_closest_planes(chain, simplify_start, simplify_end, nres, planes, by_area);

	if (simplify_chain(nres,chain,planes,simplify_start,simplify_end,by_area,nres)) {
		knotfind = false;  //if chain was simplified, then there is no knot
	}
	else { //chain might be knotted
		std::cout << "chain might have a knot, double-checking..." << std::endl;

		//then run knotfind again using the areas of the triangles to calc. the distances
		by_area = true; //now when sorting which planes to check first, use area of plane
		for ( int i = 1; i <= nres; ++i ) {
			chain(4,i) = 0;
		}
		for ( int i = 1; i <= end_trim; ++i ) {
			chain(4,i) = 1;
			simplify_start = i+1;
			chain(4,nres+1-i) = 1;
			simplify_end = nres-i;
		}
		Pairs planes2;
		find_closest_planes(chain, simplify_start, simplify_end, nres, planes2, by_area);
		if (simplify_chain(nres,chain,planes2,simplify_start,simplify_end,by_area,nres)) {
			knotfind = false;  //if chain was simplified, then there is no knot
		}
		else {
			for ( int j = 1; j <= nres; ++j ) { //see what residues did not get simplified
				if (chain(4,j) == 0) { //then point was not removed, report it
					std::cout << "Residue" << SS(j) << " remains" << std::endl;
				}
			}
		}
	}
	return knotfind;

}

////////////////////////////////////////////////////////////////////////////////
/// @begin simplify_chain
///
/// @brief make one pass through chain, removing CAs that do not result
///       in knot ambiguities/conflicts
///
/// @detailed
///car this function carries out one interation of chain simplification
///car If no simplification was made, improved is returned
///car false and no further simplification is possible
///car The presence of clashes when no further simplification can be made is
///car indicative of a knot.
///
/// @param[in]   nres - in - length of chain
/// @param  chain - [in/out]? -  CA coordinates, removed residue marker
/// @param[out]   nclashes - out - number of clashes detected
/// @param[out]   improved - out - was simplification made?
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// private function for knotfind
///
/// @references
///
/// @authors car 9/5/2003
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
bool
simplify_chain(
	int const nres,
	FArray2A_float chain, // modified on output
	Pairs & planeDists, //gets modified
	int const simplify_start,
	int const simplify_end,
	bool const by_area, //is this sorting planes by area or by distance
	int const max_trouble
)
{
	//using namespace param;

	chain.dimension( 4, nres );

//car parameter

	float const tolerance = { 0.0003 }; // closest allowed approach

	//always use tolerance = 0.0003 with gcc34, fdk 9/2005

//car local
 	int res1;
 	int ntrouble, res3;
	FArray1D_float plane( 7 ); // 1-4:  transform pts into plane coord sys
	float dist1 = 0.0;  //changed in calculate_distance
	float dist2 = 0.0;  //changed in calculate_distance
	bool done1, done2;

	// 5-7:  pts comprising plane (integers)
	FArray2D_float t_above( 4, max_trouble );
	FArray2D_float t_end( 4, max_trouble );

	bool improved = false;

	int failed = 0; //this is if the top of the stack can't be simplified
                  //then start at next smallest distance
	int length =  planeDists.size();

	while (failed < length) {

		iterator iter = planeDists.begin() + failed;
		res1 = iter->res1;
		res3 = iter->res2;  //NOTE that in all cases Pair->res2 is actually the 3rd res!! :(

		get_plane(res1, res3, chain, simplify_end, plane);

		find_trouble(plane, chain, nres, t_above, t_end, tolerance, ntrouble, max_trouble);

		if ( any_cross_plane(plane, chain, nres, ntrouble, t_above, t_end,
		 tolerance, false) ) {
			++failed; //next iteration will start at next shortest distance
			improved = false;
		} else {   //then this plane can be simplifed
			chain(4,static_cast< int >(plane(6))) = 1;
			 // eliminate second residue in plane
			improved = true;
			failed = 0; //reset to check the first element of the stack

			planeDists.erase(iter); //erase current distance between 1st and 3rd residues

			//also need to remove all distances involving the second residue
			remove_distances(planeDists, plane);

			//now need to add new distances created by removing this plane
			bool res1start = false; //is res1 the first residue in the chain?
			bool res3end = false;    //is res3 the last residue in the chain?
			int old1st = res1; //save for later, will also be middle residue of a new plane
			int old3rd = res3; //save for later, will also be middle residue of a new plane

			//find previous and next residues that weren't removed
			find_new_ends_of_plane(chain, nres, res1, res3,
														 res1start, res3end, simplify_start, simplify_end);

			//check to see if entire chain has been simplified
			if ((res1 == simplify_start) && (res3 == simplify_end)) {//then stack is empty
				int res2 = res1 + 1;
				while ( chain(4,res2) == 1) {
					++res2;
				}
				if (res2 == res3) { //then all residues between start and end have been removed
					return improved; //chain is simplified
				}
				if ((old1st == simplify_start) || (old3rd == simplify_end)) {
					  //then we still have a residue in between the start and end
					  //so just add this plane to the stack, no need to calc. the distance
   					//this case was 1 2 3 4 and got rid of 1-2-3 or 2-3-4, left with 1-X-4
					planeDists.push_back(Pair(res1,res3,269.00));
				}
				else { //there are 2 residues between the start and the end, but stack is empty
					     //case was 1 2 3 4 5, removed 2-3-4, left with 1-2-4-5 but empty stack

					calculate_distance(chain,nres,res1,old1st,old3rd,dist1,by_area); //dist 1-2-4
					calculate_distance(chain,nres,old1st,old3rd,res3,dist2,by_area); //dist 2-4-5
					add_to_end_in_order(dist1, dist2, res1, res3, old1st, old3rd, planeDists);
				}
			}
			else { // Otherwise calculate distances for the new planes
				done1 = false;
 				done2 = false;
				//do new1st residue and old3rd residue first
				if (!res1start) { //then res1 is not the first residue in the chain
					calculate_distance(chain, nres, res1, old1st, old3rd, dist1, by_area);
				}
				else { //there is no need to insert anything for res1, pretend it has been done
					done1 = true;
				}
				if (!res3end) { //then res3 is not the last residue in the chain
					calculate_distance(chain, nres, old1st, old3rd, res3, dist2, by_area);
				}
				else { //there is no need to insert anything for res3, pretend it has been done
					done2 = true;
				}

				//now insert them into the correct positions in the vector
				insert_into_correct_position(planeDists, dist1, dist2, done1, done2, res1, res3,
																		 old1st, old3rd);
			} //end dist calcs

			length = planeDists.size();  //length of stack has changed, keep it current

		} //end plane simplification

	}
  return improved;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin get_plane
///
/// @brief  given the three points, returns the plane that they create,
/// in the form: @f$ Ax + By + Cz + D = 0 @f$
///
/// @detailed  plane is created using first three non-removed CA atoms starting
///       at current residue in chain; parameter for plane equation
///       are calculated and stored in plane, as are the identity
///       of the points used in creating the plane
///
///
/// @param  res1 - [in/out]? - first residue to use in creating plane
/// @param  res3 - [in/out]? - first residue to use in creating plane
/// @param[in]   chain - in - CA coordinates, removed residue marker
/// @param[out]   nres - out - length of chain
/// @param[out]   plane - out - equation of plane and points used
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// private function for knotfind
///
///  if res1 is marked as removed in chain, the value will
///  be modified, increasing the value to the first un-removed residue in
///  chain
///
/// chain 1-4:  plane equation coefficients; 5-7 (int) points used to define plane
///
/// @references
///
/// @authors car 9/6/2003
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
void
get_plane(
	int & res1, // may be modified
	int & res3, // may be modified
	FArray2A_float chain, // CA trace
	int const nres, // chain length
	FArray1A_float plane // first 4 slots: A, B, C, D
)
{
	chain.dimension( 4, nres );
	plane.dimension( 7 );

// initialize
	plane = 0.0f;
	plane(7) = nres+1; // set high in case we return early

// now extract first residue
	float res1_x = chain(1,res1);
	float res1_y = chain(2,res1);
	float res1_z = chain(3,res1);

// find second residue
	int res2 = res1 + 1;
	if ( res2 >= nres ) return;
	while ( chain(4,res2) == 1 ) {
		++res2;
		if ( res2 >= nres ) return;
	}

// extract second residue
	float res2_x = chain(1,res2);
	float res2_y = chain(2,res2);
	float res2_z = chain(3,res2);

// use third residue
	float res3_x = chain(1,res3);
	float res3_y = chain(2,res3);
	float res3_z = chain(3,res3);

// given the three points, returns the plane that they create,
// in the form:  Ax + By + Cz + D = 0

// this is A
	plane(1) = (res1_y * (res2_z - res3_z)) + (res2_y *
	 (res3_z - res1_z)) + (res3_y * (res1_z - res2_z));

// this is B
	plane(2) = (res1_z * (res2_x - res3_x)) + (res2_z *
	 (res3_x - res1_x)) + (res3_z * (res1_x - res2_x));

// this is C
	plane(3) = (res1_x * (res2_y - res3_y)) + (res2_x *
	 (res3_y - res1_y)) + (res3_x * (res1_y - res2_y));

// this is D
	plane(4) = -((res1_x * ((res2_y * res3_z) - (res3_y * res2_z))) +
	 (res2_x * ((res3_y * res1_z) - (res1_y * res3_z))) +
	 (res3_x * ((res1_y * res2_z) - (res2_y * res1_z))));

// note the three points that make the plane
	plane(5) = res1;
	plane(6) = res2;
	plane(7) = res3;

}

////////////////////////////////////////////////////////////////////////////////
/// @begin find_trouble
///
/// @brief  find all CA-CA bonds in a chain that intersect a plane
///
/// @detailed  for all non-removed residues in chain, identify line segments
///       between sequential residues that intersect the infinite
///       plane defined by a triangle of three residues
///       used to identify potential knot conflicts in a chain, (CA-CA
///       bonds that might intersect the current triangle of interest)
///
///
/// @param[in]   plane - in - current triangle of interest
/// @param[in]   chain - in -  CA trace
/// @param[in]   nres - in - length of CA trace
/// @param[out]   t_start - out - endpoints of line segments intersecting plane
/// @param[out]   t_end - out - endpoints of line segments intersecting plane
/// @param[in]   tolerance - in - how close segment can approach plane before it 'intersects'
/// @param[out]   counter - out - logical length of t_start,t_end
/// @param[in]   max_trouble - in - physical size of t_start,t_end
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// private function for knotfind
///
///  t_start(i) & t_end(i) define a segment that crosses the plane
///  there is no guarantee that all points in t_start lie on the same side
///  of the plane
///
/// @references
///
/// @authors car 9/5/2003
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
void
find_trouble(
	FArray1A_float plane,
	FArray2A_float chain,
	int const nres, // chain length
	FArray2A_float t_start, // pairs of points describing
	FArray2A_float t_end,
	 // line segments that may intersect the plane
	float tolerance, // closest allowed approach
	int & ntrouble, // logical size of t_start,t_end
	int max_trouble // physical size of t_start,t_end
)
{
	plane.dimension( 7 );
	chain.dimension( 4, nres );
	t_start.dimension( 4, max_trouble );
	t_end.dimension( 4, max_trouble );

//car local
	bool in_plane;
	float S_of_1st_point, S_of_2nd_point;

	ntrouble = 0;
	int point1 = 1;
	while ( true ) {
		while ( chain(4,point1) == 1 ) { // point has been
			++point1; // removed
			if ( point1 >= nres ) return;
		}                                  // so ignore it
		if ( point1 >= nres ) return;
		int point2 = point1 + 1;
		while ( chain(4,point2) == 1 ) { // same thing
			++point2; // as above
			if ( point2 > nres ) return;
		}
		if ( point2 > nres ) return;

		is_part_of_plane(point1,point2,plane,in_plane);
		if ( !in_plane ) {
			S_of_1st_point = (plane(1) * chain(1,point1)) +
			 (plane(2) * chain(2,point1)) +
			 (plane(3) * chain(3,point1)) + plane(4);

// S equals Ax+By+Cz+D, this is done to the 1st point and to the 2nd
// using the x,y,z values of point 1 and point 2 respectively


//car S_of_point1 = distance from point 1 to the plane
			S_of_2nd_point = (plane(1) * chain(1,point2)) +
			 (plane(2) * chain(2,point2)) +
			 (plane(3) * chain(3,point2)) + plane(4);

// if both S are on the same side of the triangle plane (ie the xy plane),
//  the segment cant intersect the triangle

			if ( (S_of_1st_point <= tolerance && S_of_2nd_point >=-tolerance ) ||
			 (S_of_1st_point >= tolerance && S_of_2nd_point <= -tolerance )) {

// if the two points create a line segment that potentially cross the triangle
// (ie one above the xy plane, one below)
// mark them as trouble and keep them to check whether they actually cross
// note that t_start does not necessarily contain all points above the
// xy plane
				++ntrouble;
				if ( ntrouble > max_trouble ) {
					std::cout << "STOPPING:: increase max_trouble " << std::endl;
					std::cout << " in function find_trouble, knots.cc" << std::endl;
					utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
				}

				t_start(1,ntrouble) = chain(1,point1);
				t_start(2,ntrouble) = chain(2,point1);
				t_start(3,ntrouble) = chain(3,point1);
				t_start(4,ntrouble) = point1;
				t_end(1,ntrouble) = chain(1,point2);
				t_end(2,ntrouble) = chain(2,point2);
				t_end(3,ntrouble) = chain(3,point2);
				t_end(4,ntrouble) = point2;

			}
		}
		point1 = point2;
	}

}


////////////////////////////////////////////////////////////////////////////////
/// @begin is_part_of_plane
///
/// @brief returns true if either point1 or point2 are one of the points
///       used to define plane
///
/// @detailed
///
/// @param[in]   point1 - in - first point number to check
/// @param[in]   point2 - in - secont point number to check
/// @param[in]   plane - in - plane equation coefficients & point numbers used to define plane
/// @param[out]   in_plane - out - is point1 or point2 a plane-defining point?
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// private function for knotfind
///
/// @references
///
/// @authors car 9/5/2003
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
void
is_part_of_plane(
	int & point1,
	int & point2,
	FArray1A_float plane,
	bool & in_plane
)
{
	plane.dimension( 7 );

	in_plane = false;
	if ( ( point1 == plane(5) ) || ( point1 == plane(6) ) ||
		  ( point1 == plane(7) ) || ( point2 == plane(5) ) ||
		  ( point2 == plane(5) ) || ( point2 == plane(7) ) ) in_plane = true;

}


////////////////////////////////////////////////////////////////////////////////
/// @begin any_cross_plane
///
/// @brief
/// returns true if any segments described by pairs of points in t_start,t_end,
/// intersect with a triangle defined in plane
///
/// @detailed checks line segments that intersect the infinite plane to see if the
///      intersection point lies within the triangle of interest
///
/// @param[in]   plane - in -  plane equation coefficient, and points numbers defining triangle
/// @param[in]   chain - in - CA coordinates and removed residue marker
/// @param[in]   nres - in - length of chain
/// @param[in]   ntrouble - in - logical size of t_start,t_end
/// @param[in]   t_start - in - endpoints of segments intersecting infinite plane
/// @param[in]   t_end - in - endpoints of segments intersecting infinite plane
/// @param[in]   tolerance - in - how close can intersection point and triangle approach?
/// @param[in]   check_all- in - if true, check all trouble vectors and return all clashing vectors
///
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// private function for knotfind
/// modifies ntrouble, t_start,t_end to reflect vectors that have been found
/// to conflict
///
/// @references
///
/// @authors car 9/5/2003
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
bool
any_cross_plane(
	FArray1A_float plane,
	FArray2A_float chain,
	int const nres, // logical size of chain
	int & ntrouble,
	FArray2A_float t_start,
	FArray2A_float t_end,
	float tolerance, // closest allowed approach
	bool const check_all  // check all trouble vectors?
)
{
	plane.dimension( 7 );
	chain.dimension( 4, nres );
	t_start.dimension( 4, ntrouble );
	t_end.dimension( 4, ntrouble );

// local
	FArray1D_float i_point( 3 );
	float denom, mu;

	int ncheck = ntrouble;
	ntrouble = 0;

	bool intersects = false;
	for ( int k = 1; k <= ncheck; ++k ) {
//car find intersection of segment with infinite plane:
		denom =
		 (plane(1) * (t_start(1,k) - t_end(1,k))) +
		 (plane(2) * (t_start(2,k) - t_end(2,k))) +
		 (plane(3) * (t_start(3,k) - t_end(3,k)));

//car check for division by zero--
		if ( std::abs(denom) <= std::max(.0001f,tolerance) ) {
//car zero denom means segment is parallel to plane; we have already
//car excluded pairs that are both above or below plane, so zero denom
//car must mean that segment lies in plane-- must check this segment
//car for intersection (ie segment crosses at least one triangle edge
			if ( segment_intersects_triangle(t_start(1,k),t_end(1,k),
			 chain(1,static_cast< int >(plane(5))),
			 chain(1,static_cast< int >(plane(6))),
			 chain,tolerance ) ) {
				intersects = true;
				++ntrouble;
				t_start(1,ntrouble) = t_start(1,k);
				t_start(2,ntrouble) = t_start(2,k);
				t_start(3,ntrouble) = t_start(3,k);
				t_start(4,ntrouble) = t_start(4,k);
				t_end(1,ntrouble) = t_end(1,k);
				t_end(2,ntrouble) = t_end(2,k);
				t_end(3,ntrouble) = t_end(3,k);
				t_end(4,ntrouble) = t_end(4,k);
				//std::cout << "segment intersects triangle! " << std::endl; //debug
				if (! check_all) return intersects;
			}
			goto L100;
		}

		mu = ( plane(4) + (plane(1) * (t_start(1,k))) +
		 (plane(2) * (t_start(2,k))) + (plane(3) * (t_start(3,k))) ) / denom;

//car if 0 < mu < 1, then the line segment intersects the infinite
//car plane and more checking is necessary. Otherwise ignore it.
		if ( mu < tolerance || mu > 1.0+tolerance ) goto L100;

		i_point(1) = ( t_start(1,k) + ( mu * (t_end(1,k) - t_start(1,k)) ) );
		i_point(2) = ( t_start(2,k) + ( mu * (t_end(2,k) - t_start(2,k)) ) );
		i_point(3) = ( t_start(3,k) + ( mu * (t_end(3,k) - t_start(3,k)) ) );


		if ( inside_triangle(plane, chain, nres, i_point, tolerance )) {
			intersects = true;
			++ntrouble;
			t_start(1,ntrouble) = t_start(1,k);
			t_start(2,ntrouble) = t_start(2,k);
			t_start(3,ntrouble) = t_start(3,k);
			t_start(4,ntrouble) = t_start(4,k);
			t_end(1,ntrouble) = t_end(1,k);
			t_end(2,ntrouble) = t_end(2,k);
			t_end(3,ntrouble) = t_end(3,k);
			t_end(4,ntrouble) = t_end(4,k);
			if (! check_all) return intersects;
		}
L100:;
	}
	return intersects;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin inside_triangle
///
/// @brief
/// returns true  if  i_point lies within the triangular facet
/// defined by plane; i_point is assumed to lie on the inifinite plane
/// described by plane
///
/// @detailed
///
/// @param[in]   plane - in - plane equation coeffiecients ant point numbers of triangle
/// @param[in]   ca_chain - in - ca coordinates and removed residue marker
/// @param[in]   nres - in - length of ca_chain
/// @param[in]   i_point - in -  point of interest
/// @param[out]   tolerance - out - closest allowed approach of point and triangle
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// private function for knotfind
///
///car determined using the fact that the sum of internal angles will
///car be 2*pi for points inside the triangle and <2*pi for those outside.
///
/// @references
///
/// @authors car 9/5/2003
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
bool
inside_triangle(
	FArray1A_float plane, // 1-4 = equation of inifinite plane  Ax + By + Cz + D
	FArray2A_float ca_chain,
	int const nres,
	FArray1A_float i_point, // point lying on infinite plane
	float tolerance // closest allowed approach
)
{
	using namespace numeric::constants::f;
	using numeric::sin_cos_range;

	plane.dimension( 7 );
	ca_chain.dimension( 4, nres );
	i_point.dimension( 3 );

//car local
	float norm_factor, value, angle1, angle2, angle3;
	FArray1D_float Pa1( 3 );
	FArray1D_float Pa2( 3 );
	FArray1D_float Pa3( 3 );

	bool in_triangle = false;

// points of triangle // type conversions
	int pa = static_cast< int >( plane(5) );
	int pb = static_cast< int >( plane(6) );
	int pc = static_cast< int >( plane(7) );

// determine the line seqment Pa1

	Pa1(1) = (ca_chain(1,pa) - i_point(1));
	Pa1(2) = (ca_chain(2,pa) - i_point(2));
	Pa1(3) = (ca_chain(3,pa) - i_point(3));
// normalize it
	norm_factor = 1.0 / std::sqrt(
	 ( Pa1(1) * Pa1(1) ) + ( Pa1(2) * Pa1(2) ) + ( Pa1(3) * Pa1(3) ) );
	Pa1(1) = Pa1(1) * norm_factor;
	Pa1(2) = Pa1(2) * norm_factor;
	Pa1(3) = Pa1(3) * norm_factor;

// determine the line seqment Pa2

	Pa2(1) = ca_chain(1,pb) - i_point(1);
	Pa2(2) = ca_chain(2,pb) - i_point(2);
	Pa2(3) = ca_chain(3,pb) - i_point(3);
// normalize it
	norm_factor = 1.0 / std::sqrt(
	 ( Pa2(1) * Pa2(1) ) + ( Pa2(2) * Pa2(2) ) + ( Pa2(3) * Pa2(3) ) );
	Pa2(1) = Pa2(1) * norm_factor;
	Pa2(2) = Pa2(2) * norm_factor;
	Pa2(3) = Pa2(3) * norm_factor;

// determine the line seqment Pa3

	Pa3(1) = ca_chain(1,pc) - i_point(1);
	Pa3(2) = ca_chain(2,pc) - i_point(2);
	Pa3(3) = ca_chain(3,pc) - i_point(3);
// normalize it
	norm_factor = 1.0 / std::sqrt(
	 ( Pa3(1) * Pa3(1) ) + ( Pa3(2) * Pa3(2) ) + ( Pa3(3) * Pa3(3) ) );
	Pa3(1) = Pa3(1) * norm_factor;
	Pa3(2) = Pa3(2) * norm_factor;
	Pa3(3) = Pa3(3) * norm_factor;

// determine 3 angles formed within the triangle by the 3 lines

	angle1 = ( Pa1(1) * Pa2(1) ) + ( Pa1(2) * Pa2(2) ) + ( Pa1(3) * Pa2(3) );
	angle2 = ( Pa2(1) * Pa3(1) ) + ( Pa2(2) * Pa3(2) ) + ( Pa2(3) * Pa3(3) );
	angle3 = ( Pa3(1) * Pa1(1) ) + ( Pa3(2) * Pa1(2) ) + ( Pa3(3) * Pa1(3) );

// determine, from total angle, if line segment intersects facet

	value =
	 std::acos( sin_cos_range( angle1 ) ) +
	 std::acos( sin_cos_range( angle2 ) ) +
	 std::acos( sin_cos_range( angle3 ) );

// if value equals 2 pi, then point lies inside facet

	if ( std::abs( value - pi_2 ) < std::max(.0001f,tolerance) ) in_triangle = true;
	return in_triangle;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin segment_intersects_triangle
///
/// @brief
///car does segment x-y intersect triangle a,b,c?
///car x,y,a,b,c are all co-planar
///
/// @detailed
///
/// @param[in]   x - in - one of two points defining a line segment
/// @param[in]   y - in - one of two points defining a line segment
/// @param[in]   a - in - one of three points defining a triangle
/// @param[in]   b - in - one of three points defining a triangle
/// @param[in]   c - in - one of three points defining a triangle
/// @param[in]   tolerance - in - closest allowed approach of segment and triangle
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors  car 9/4/2003
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
bool
segment_intersects_triangle(
	FArray1A_float x,
	FArray1A_float y,
	FArray1A_float a,
	FArray1A_float b,
	FArray1A_float c,
	float tolerance // closest allowed approach of segments
)
{
	x.dimension( 3 );
	y.dimension( 3 );
	a.dimension( 3 );
	b.dimension( 3 );
	c.dimension( 3 );

//car local
	FArray1D_float xy( 3 ); // xy vector
	FArray1D_float ab( 3 ); // ab vector
	FArray1D_float x1( 3 );
	FArray1D_float y1( 3 );
	FArray1D_float a1( 3 );
	FArray1D_float b1( 3 );
	FArray1D_float c1( 3 );
	FArray2D_float mat( 3, 3 );
	float norm;
	float upper, lower, limit;
	float intercept;

//car set up coordinate system with x at origin, xy-> x-axis, ab in xy plane

	xy(1) = y(1)-x(1);
	xy(2) = y(2)-x(2);
	xy(3) = y(3)-x(3);

	ab(1) = b(1)-a(1);
	ab(2) = b(2)-a(2);
	ab(3) = b(3)-a(3);

//  ab cross xy, normalize -> z axis

	mat(1,3) = ab(2)*xy(3) - ab(3)*xy(2);
	mat(2,3) = ab(3)*xy(1) - ab(1)*xy(3);
	mat(3,3) = ab(1)*xy(2) - ab(2)*xy(1);

	norm = std::sqrt(mat(1,3)*mat(1,3)+mat(2,3)*mat(2,3)+mat(3,3)*mat(3,3));
	if ( norm == 0.0 ) norm = 0.0001;
	norm = 1.0/norm;
	mat(1,3) *= norm;
	mat(2,3) *= norm;
	mat(3,3) *= norm;

//car normalize x axis
	norm = std::sqrt(xy(1)*xy(1)+xy(2)*xy(2)+xy(3)*xy(3));
	if ( norm == 0.0 ) norm = 0.0001;
	norm = 1.0/norm;
	mat(1,1) = norm*xy(1);
	mat(2,1) = norm*xy(2);
	mat(3,1) = norm*xy(3);


//car z cross x -> y axis
	mat(1,2) = mat(2,3)*mat(3,1) - mat(3,3)*mat(2,1);
	mat(2,2) = mat(3,3)*mat(1,1) - mat(1,3)*mat(3,1);
	mat(3,2) = mat(1,3)*mat(2,1) - mat(2,3)*mat(1,1);


//car rotate points
	for ( int i = 1; i <= 3; ++i ) {
		x1(i) = mat(i,1)*x(1)+mat(i,2)*x(2)+mat(i,3)*x(3);
		y1(i) = mat(i,1)*y(1)+mat(i,2)*y(2)+mat(i,3)*y(3);
		a1(i) = mat(i,1)*a(1)+mat(i,2)*a(2)+mat(i,3)*a(3);
		b1(i) = mat(i,1)*b(1)+mat(i,2)*b(2)+mat(i,3)*b(3);
		c1(i) = mat(i,1)*c(1)+mat(i,2)*c(2)+mat(i,3)*c(3);
	}
//car offset points
	for ( int i = 1; i <= 3; ++i ) {
		a1(i) -= x1(i);
		b1(i) -= x1(i);
		c1(i) -= x1(i);
		y1(i) -= x1(i);
	}

//car x1 at 0,0,0; y1 at y1(1),0,0

	bool intersects = true;

	lower = 0.0 - tolerance;
	upper = y1(1) + tolerance;
	limit = std::max(0.0001f,tolerance);
//car do vertices lie on segment?
	if ( std::abs(a1(2)) <= limit && a1(1) >= lower && a1(1) <= upper ) return intersects;

	if ( std::abs(b1(2)) <= limit && b1(1) >= lower && b1(1) <= upper ) return intersects;

	if ( std::abs(c1(2)) <= limit && c1(1) >= lower && c1(1) <= upper ) return intersects;

//car do edges cross segment?
	if ( ( a1(2) >= -tolerance && b1(2) <= tolerance ) ||
	 ( a1(2) <= tolerance && b1(2) >= -tolerance ) ) {
		intercept = a1(1) - a1(2)*(a1(1)-b1(1))/(a1(2)-b1(2));
		if ( intercept >= lower && intercept <= upper ) return intersects;
	}
	if ( ( b1(2) >= -tolerance && c1(2) <= tolerance ) ||
	 ( b1(2) <= tolerance && c1(2) >= -tolerance ) ) {
		intercept = b1(1) - b1(2)*(b1(1)-c1(1))/(b1(2)-c1(2));
		if ( intercept >= lower && intercept <= upper ) return intersects;
	}
	if ( ( c1(2) >= -tolerance && a1(2) <= tolerance ) ||
	 ( c1(2) <= tolerance && a1(2) >= -tolerance ) ) {
		intercept = c1(1) - c1(2)*(c1(1)-a1(1))/(c1(2)-a1(2));
		if ( intercept >= lower && intercept <= upper ) return intersects;
	}

	intersects = false;
	return intersects;
}

//////////////////////////////////////////////////////////////////////////////////
/// @begin find_closest_planes
///
/// @brief go through entire chain and calculate all distances and sort vector
///
/// @detailed iterate through entire chain, three consecutive residues at a time,
///           calculate all distances between the first residues and the third res
///           or (if the flag "by_area" is true)
///           calculate the area of all the triangles
///
/// @param[in]   chain - in - CA coordinates, removed residue marker
/// @param[in]   start - in - first residue of the chain
/// @param[in]   end - in - last residue of the chain
/// @param[in]   nres - in - length of chain
/// @param[out]  planeDists - out - vector residue pairs and their distances
/// @param[in]   by_area - in - flag to denote which way to calculate the distance
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// private function for knotfind
///
/// distance in this function could be the distance between the first and third
/// residues OR the area of the triangle made up of the first, second, and third
/// residues (depending on which version of knotfind is being run: flag = by_area)
///
/// distance is actually always squared:
///          either the distance between first_res and third_res squared
///          or the area of the triangle squared
///
/// @references
///
/// @authors fdk 10/3/2005
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////

void
find_closest_planes(
 	FArray2A_float chain,
	int start,
	int end,
	int const nres,
	Pairs & planes,
	bool const by_area //is this finding planes by area or by distance
)
{

	chain.dimension( 4, nres );
	float distance;
	int i,j,k;
	float a,b,c,s;

	for (i = start; i <= (end - 2); ++i) { //go through chain and calc distances
		j = i + 2;  //third residue
		if (!by_area) { //then we are not calculating distances by the area of the plane
			//distance is from the first residue to the third residue
			distance = ( square( chain(1,i) - chain(1,j) ) +
									 square( chain(2,i) - chain(2,j) ) +
									 square( chain(3,i) - chain(3,j) ) );
		}
		else { //then we are calculating the distance as the area of the triangle
			//calculate the area between res1->res2->res3 using Heron's Formula
			k = i + 1; //second residue
			//the distance between res1 and res2:
			a = std::sqrt( square( chain(1,i) - chain(1,k) ) +
										 square( chain(2,i) - chain(2,k) ) +
										 square( chain(3,i) - chain(3,k) ) );
			//the distance between res2 and res3:
			b = std::sqrt( square( chain(1,k) - chain(1,j) ) +
										 square( chain(2,k) - chain(2,j) ) +
										 square( chain(3,k) - chain(3,j) ) );
			//the distance between res1 and res3:
			c = std::sqrt( square( chain(1,i) - chain(1,j) ) +
										 square( chain(2,i) - chain(2,j) ) +
										 square( chain(3,i) - chain(3,j) ) );
			s = ((a+b+c) / 2 );
			distance = (s * (s-a) * (s-b) * (s-c));  //this is the area_squared
			assert (distance >= 0.0); //make sure the area squared is positive
			//don't bother taking the square root of the area
		}
			planes.push_back(Pair(i,j,distance)); //push it onto the vector
	}

	std::sort(planes.begin(),planes.end(),compare_Pair);

}

//////////////////////////////////////////////////////////////////////////////////
/// @begin insert_into_correct_position
///
/// @brief insert two new distances into the stack in the correct spots
///
/// @detailed iterates through the vector to insert two new distances
///           into their correct positions
///
/// @param       planeDists - [in/out] - vector residue pairs and their distances
/// @param[in]   dist1 - in - the distance used for res1 and old3rd
/// @param[in]   dist2 - in - the distance used for old1st and res3
/// @param[in]   done1 - in - res1 is the 1st res in the chain, no need to insert
/// @param[in]   done2 - in - res3 is the last res in the chain, no need to insert
/// @param[in]   res1 - in - first residue to use in creating plane
/// @param[in]   res3 - in - third residue to use in creating plane
/// @param[in]   old1st - in - previously the first residue
/// @param[in]   old3rd - in - previously the third residue
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// private function for knotfind
///
/// distance in this function could be the distance between the first and third
/// residues OR the area of the triangle made up of the first, second, and third
/// residues (depending on which version of knotfind is being run: flag = by_area)
///
/// @references
///
/// @authors fdk 10/3/2005
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////

void
insert_into_correct_position(
	Pairs & planeDists,
	float const dist1,
	float const dist2,
	bool done1,
	bool done2,
	int const res1,
	int const res3,
	int const old1st,
	int const old3rd
)
{
	for ( iterator iterate = planeDists.begin(); iterate != planeDists.end(); ++iterate ) {
		if ( ( iterate->dist > dist1 ) && !done1 ) {
			planeDists.insert( iterate, Pair(res1,old3rd,dist1) );
			done1 = true;
			break;
		}
	}
	for ( iterator iterate = planeDists.begin(); iterate != planeDists.end(); ++iterate ) {
		if ( ( iterate->dist > dist2 ) && !done2 ) {
			planeDists.insert( iterate, Pair(old1st,res3,dist2) );
			done2 = true;
			break;
		}
	}
	if ( !( done1 && done2 ) ) { // then either one or both distances were not added
		//	the distance to insert is bigger than the last distance on the stack?
		if ( !done1 && !done2 ) { // then both distances need to be added
			add_to_end_in_order( dist1, dist2, res1, res3, old1st, old3rd, planeDists );
		} else { // only 1 distance needs to be added
			if ( !done1 ) { // dist1 needs to be added
				planeDists.push_back( Pair(res1,old3rd,dist1) );
			} else { // dist2 needs to be added
				planeDists.push_back( Pair(old1st,res3,dist2) );
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////
/// @begin add_to_end_in_order
///
/// @brief insert the two new distances to the end of the stack
///
/// @detailed inserts two new distances the the back of the vector
///           in the correct order
///
/// @param[in]   dist1 - in - the distance used for res1 and old3rd
/// @param[in]   dist2 - in - the distance used for old1st and res3
/// @param[in]   res1 - in - first residue to use in creating plane
/// @param[in]   res3 - in - third residue to use in creating plane
/// @param[in]   old1st - in - previously the first residue
/// @param[in]   old3rd - in - previously the third residue
/// @param       planeDists - [in/out] - vector residue pairs and their distances
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// private function for knotfind
///
/// distance in this function could be the distance between the first and third
/// residues OR the area of the triangle made up of the first, second, and third
/// residues (depending on which version of knotfind is being run: flag = by_area)
///
/// @references
///
/// @authors fdk 10/3/2005
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////

void
add_to_end_in_order(
	float const dist1,
	float const dist2,
	int const res1,
	int const res3,
	int const old1st,
	int const old3rd,
	Pairs & planeDists
	)
{

	if (dist1 < dist2) {
		planeDists.push_back(Pair(res1,old3rd,dist1));
		planeDists.push_back(Pair(old1st,res3,dist2));
	}
	else {
		planeDists.push_back(Pair(old1st,res3,dist2));
		planeDists.push_back(Pair(res1,old3rd,dist1));
	}

}

/////////////////////////////////////////////////////////////////////////////////
/// @begin calculate_distance
///
/// @brief calculate the correct distance to be used for 1st residue and 3rd residue
///
/// @detailed calculate the distance between the first residue and the third residue
///           or (if the flag "by_area" is true)
///           calculate the area of triangle first_res, second_res, third_res
///
/// @param[in]   chain - in - CA coordinates, removed residue marker
/// @param[in]   nres - in - length of chain
/// @param[in]   first_res - in - first residue to use in creating plane
/// @param[in]   second_res - in - second residue to use in creating plane
/// @param[in]   third_res - in - third residue to use in creating plane
/// @param[out]  dist - out - the distance used for first_res and third_res
/// @param[in]   by_area - in - flag to denote which way to calculate the distance
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// private function for knotfind
///
/// distance in this function could be the distance between the first and third
/// residues OR the area of the triangle made up of the first, second, and third
/// residues (depending on which version of knotfind is being run: flag = by_area)
///
/// distance is actually always squared:
///          either the distance between first_res and third_res squared
///          or the area of the triangle squared
///
/// @references
///
/// @authors fdk 10/3/2005
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////

void
calculate_distance(
	FArray2A_float chain,
	int const nres,      // logical size of chain
	int const first_res,
	int const second_res,  //only used if by_area is true
	int const third_res,
	float & dist,        // will be changed on output
	bool const by_area   // is this sorting planes by area or by distance
)
{

	chain.dimension( 4, nres );
	float a,b,c,s;

	if (!by_area) { //then just calculate the distance between first_res and third_res
		dist = ( square( chain(1,first_res) - chain(1,third_res) ) +
						 square( chain(2,first_res) - chain(2,third_res) ) +
						 square( chain(3,first_res) - chain(3,third_res) ) );
	}
	else { //then calculate the area between first_res, second_res, and third_res
		a = std::sqrt( square( chain(1,first_res) - chain(1,second_res) ) +
									 square( chain(2,first_res) - chain(2,second_res) ) +
									 square( chain(3,first_res) - chain(3,second_res) ) );
		b = std::sqrt( square( chain(1,second_res) - chain(1,third_res) ) +
									 square( chain(2,second_res) - chain(2,third_res) ) +
									 square( chain(3,second_res) - chain(3,third_res) ) );
		c = std::sqrt( square( chain(1,first_res) - chain(1,third_res) ) +
									 square( chain(2,first_res) - chain(2,third_res) ) +
									 square( chain(3,first_res) - chain(3,third_res) ) );
		s = ((a+b+c) / 2 );
		dist = (s * (s-a) * (s-b) * (s-c)); //dist is the area_squared
		assert (dist >= 0.0); //make sure area_squared is positive
    //don't waste time taking root
	}

}

//////////////////////////////////////////////////////////////////////////////////
/// @begin find_new_ends_of_plane
///
/// @brief find previous & next residues that weren't removed
///
/// @detailed after the second residue has been removed, need to find a new first
///           residue as well as a new third residue
///
/// @param[in]   chain - in - CA coordinates, removed residue marker
/// @param[in]   nres - in - length of chain
/// @param  res1 - [in/out] - first residue to use in creating plane
/// @param  res3 - [in/out] - third residue to use in creating plane
/// @param[out]  res1start - out - flag that res1 is the first residue in the chain
/// @param[out]  res3end - out - flag that res3 is the last residue in the chain
/// @param[in]   simplify_start - in - the first residue of the chain being examined
/// @param[in]   simplify_end - in - the last residue of the chain being examined
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// private function for knotfind
///
/// example: previously had residues 1 2 3 4 5 7 (6 was removed already)
///          just removed residue 4 from triangle 3-4-5
///          now need to find triangle 2-3-5 and 3-5-7
///          res1 used to be 3, will now be 2. old1st will still be 3.
///          res3 was 5, will now be 7.        old3rd will still be 5.
///
/// @references
///
/// @authors fdk 10/3/2005
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////

void
find_new_ends_of_plane(
	FArray2A_float chain,
	int const nres,    // logical size of chain
	int & res1,        // may be modified
	int & res3,        // may be modified
	bool & res1start,  // may be modified
	bool & res3end,    // may be modified

	int const simplify_start,
  int const simplify_end
)
{
	chain.dimension( 4, nres );

	if (res1 != simplify_start) { //res1 does not change if it is already the 1st res
		--res1; //go back to find previous residue that wasn't removed
		while ( chain(4,res1) == 1) {
			--res1; //will always stop at simplify_start (which is never == 1)
		}
	}
	else { //res1 is the first residue in the chain
		res1start = true;
	}

	//now do the same for the 3rd residue

	if (res3 != simplify_end) { //res3 does not change if it is already the last res
		++res3; //find next residue that hasn't been removed already
		while ( chain(4,res3) == 1) {
			++res3; //will always stop at simplify_end (which is never == 1)
		}
	}
	else { //res3 is the last residue in the chain
		res3end = true;
	}
}

//////////////////////////////////////////////////////////////////////////////////
/// @begin remove_distances
///
/// @brief remove any distances involving the second residue
///
/// @detailed need to go through the vector of distances and remove all distances
///           involving the second residue
///
/// @param    planeDists - [in/out] - vector residue pairs and their distances
/// @param[in]   plane - in - equation of plane and points used
///
/// @global_read
///
/// @global_write
///
/// @remarks
/// private function for knotfind
///
/// note that in the Pairs vector res2 is actually the third residue
///
/// @references
///
/// @authors fdk 10/3/2005
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////

void
remove_distances(
  Pairs & planeDists,   //priority queue of plane distances
  FArray1A_float plane //equation of plane and 3 residues, plane(6) is the 2nd residue
)
{
	plane.dimension( 7 );

	for (iterator cycle = planeDists.begin(); cycle != planeDists.end(); cycle++) {
		if (cycle->res1 == plane(6)) {
			planeDists.erase(cycle); //erase case where res1 = second residue just removed
			break;
		}
	}
	//need to do seperately or else problems occur if they are consecutive in the stack
	//ie used to have one check (if res1=plane6 || res2==plane6) wigs out if consec.
	for (iterator cycle = planeDists.begin(); cycle != planeDists.end(); cycle++) {
		if (cycle->res2 == plane(6)) {  //NOTE that cycle->res2 is actually the 3rd res!!
			planeDists.erase(cycle); //erase case where 3rdRes = second residue just removed
			break;
		}

	}
}



bool
knotfind(
		 const core::pose::Pose & pose,
		 core::Size start,
		 core::Size end
		 )
{
	if (start > end || start < 1 || start > pose.total_residue()) {
		utility_exit_with_message("Bad indices to knotfind");
	}

	core::Size nres = end - start + 1;
	
	ObjexxFCL::FArray2D_float chain(4, nres);
	for (core::Size ii = start; ii <= end; ++ ii) {
		const core::Vector & pos = pose.residue(ii).nbr_atom_xyz();
		core::Size chain_ii = ii - start + 1;
		
		chain(1, chain_ii) = pos.x();
		chain(2, chain_ii) = pos.y();
		chain(3, chain_ii) = pos.z();
		chain(4, chain_ii) = 0.0;
	}
	
	return knotfind(chain, nres, start, end, 0);
}



} /* filter */
} /* util */
} /* interactive */
