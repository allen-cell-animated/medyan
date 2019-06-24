/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

A Coords structure that holds coordinates, bead indices and a number of auxiliary variables.
*/

#ifndef DIST_COORDS
#define DIST_COORDS

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <initializer_list>

#include "dist_moduleV2/dist_common.h"

namespace dist {

	// This structure accepts coordinates and indicies (or could mock generate them for testing purposes)
	struct Coords {
		// --- Variables ------------------	
		// 
		// public

		std::vector<int> indices;
		std::vector<float> x, y, z;
	
		// --- Methods ------------------
	
		template <typename T, typename U>
		void init_coords(const std::vector<T> &xx, const std::vector<T> &yy, const std::vector<T> &zz, const std::vector<U> &indx)		
		{
			uint N = xx.size();
			this->resize(N);

			copy(xx.begin(),xx.end(),x.begin());
			copy(yy.begin(),yy.end(),y.begin());
			copy(zz.begin(),zz.end(),z.begin());
			copy(indx.begin(),indx.end(),indices.begin());
		}
		
		void _init_coords_mock(uint N)
		{
			this->resize(N);	
			create_mock_values();
			uint i=0;
			std::generate(indices.begin(),indices.end(),[&i](){return i++;});
		}
					
		Coords() = default;

		Coords(uint N)
		{
			_init_coords_mock(N);
		}

		template <typename T, typename U>
		Coords(const std::vector<T> &xx, const std::vector<T> &yy, const std::vector<T> &zz, const std::vector<U> &indx)
		{
			init_coords(xx,yy,zz,indx);
		}
	
		void resize(uint N){
			x.resize(N);
			y.resize(N);
			z.resize(N);
			indices.resize(N);
		}
			
		uint size() const {return x.size();}
	
		void create_mock_values()
		{	
			std::random_device rd;
			std::mt19937 mt(rd());
		    std::uniform_real_distribution<float> dist_d(1.0, 10.0);
		
			uint N = x.size();
	
			for(uint i=0; i<N; ++i){
				x[i] = dist_d(mt);
				y[i] = dist_d(mt);
				z[i] = dist_d(mt);
			}
		}	
	};

} // end-of-namespace dist

#endif // DIST_COORDS
	
