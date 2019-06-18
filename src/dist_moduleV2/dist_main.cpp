#include <iostream>

#include "dist_moduleV2/dist_driver.h"

using namespace std;

int main(int argc, char* argv[]){
	
	const uint N = 6000; //test
	
	// cout << "\n\nStarting tests...\n" << endl;
	// bool tests_passed = dist::test_dist_module(N,N);
	// cout << "All tests_passed? " << boolalpha << tests_passed << endl;

	cout << "\n\nStarting benchmarks...\n" << endl;
	dist::benchmark_dist_module(N);
	
	// cout << "\n\nAn example where coordinates and indices are external...\n" << endl;
	// dist::examples_dist_module();
	
	// dist::dist_example_avx_veccoord();
		
	return 0;
}

// Garyk's Macbook Pro          Deepthought II (login node)
// 1 CPU Timings
//
// SERIAL1: 60 ms              90 ms
// SERIAL2: 70 ms              105 ms
// SIMD-AVX-F1: 16 ms          35 ms
// SIMD-AVX-F2: 27 ms          50 ms
// SIMD-AVX2-F1: 12 ms
// SIMD-AVX2-F2: 19 ms
// SIMD-AVX2-I1: 9 ms
// SIMD-AVX2-I2: 15 ns
//
//
// 4 CPU Timings
//
// SIMD-AVX-F2:  8 ms           16 ms
// SIMD-AVX2-F2: 5 ms
// SIMD-AVX2-I2: 4 ns
//


