/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

A simple testing suite. Mainly checks the SIMD version results again the serial ones.
*/

#include <random>
#include <iostream>
#include <chrono>

#include "dist_moduleV2/dist_coords.h"
#include "dist_moduleV2/dist_driver.h"
#include "dist_moduleV2/dist_out.h"

namespace dist {

	using namespace std;

	int pack_two_ints(int high, int low)
	{
		return (high << 16) | low; // pack lower 16 bits of each integer into a new integer; low, high < 65536 (i.e N < 65536)
	}

	int int_low(int i){
		int mask_l = (1<<16)-1;
		return i & mask_l;
	}

	void print_int_low(int i)
	{
		cout << int_low(i) << " ";
	}

	int int_high(int i)
	{
		int mask_h = ((1<<16)-1) << 16;
		 return ((i & mask_h) >> 16);
	}
	
	void print_int_high(int i)
	{
		cout << int_high(i) << " ";
	}

	template <uint D, bool SELF>
	vector<int> convert_dout_unique(dOut<D,SELF> &out, uint d=0){

		uint ncontacts = out.counter[d];
		// cout << "ncontacts=" << ncontacts << endl;
	
		vector<int> vout(ncontacts);
		for(uint k=0; k<ncontacts; ++k){
			vout[k] = pack_two_ints(out.dout[2*d][k], out.dout[2*d+1][k]);
		}
	
		// int mask_l = (1<<16)-1;
		// int mask_h = ((1<<16)-1) << 16;
		//
		// print_bits(mask_l,32);
		// print_bits(mask_h,32);
	
		// for(uint k=1000; k<1010; ++k){
		// 	cout << "i, j = [" << out.dout[0][k] << ":" << out.dout[1][k] << "], vout=" << vout[k] << endl;
		// 	cout << "i from vout: " << (vout[k] & mask_l) << endl;
		// 	cout << "j from vout: " << ((vout[k] & mask_h) >> 16) << endl;
		// }
	
		sort(vout.begin(), vout.end());
		return vout;
	}

	double test_compare_douts(const vector<int> &vout1, const vector<int> &vout2)
	{
		bool success = false;
	
		size_t n_max = max(vout1.size(),vout2.size());
		size_t res_size = vout1.size() + vout2.size();
		vector<int> res(res_size);
		auto newend = std::set_symmetric_difference(vout1.begin(),vout1.end(),
									   vout2.begin(),vout2.end(),
									   res.begin());
	    size_t n_different = std::distance(res.begin(), newend);
		cout << "test_compare_douts(vout1, vout2): different " << n_different << " elements out of max " << n_max << " elements" << endl;
		double fraction_failure = (double)n_different/n_max;
		cout << "fraction_failure=" << fraction_failure << endl;
		// if(fraction_failure>1e-6){
		// 	vector<int> diffvec1(n_different);
		// 	auto end_diffvec1 = std::set_difference(vout1.begin(),vout1.end(),
		// 								           vout2.begin(),vout2.end(),
		// 								           diffvec1.begin());
		// 	for(auto x = diffvec1.begin(); x < end_diffvec1; ++x){
		// 		cout << "Non-Serial: ";
		// 		print_int_high(*x); print_int_low(*x);
		// 		cout << endl;
		// 	}
		//
		// 	vector<int> diffvec2(n_different);
		// 	auto end_diffvec2 = std::set_difference(vout2.begin(),vout2.end(),
		// 								           vout1.begin(),vout1.end(),
		// 								           diffvec2.begin());
		// 	for(auto x = diffvec2.begin(); x < end_diffvec2; ++x){
		// 		cout << "Serial: ";
		// 		print_int_high(*x); print_int_low(*x);
		// 		cout << endl;
		// 	}
		// }
		return 1-fraction_failure;
	
	}

	template <typename TAG>
	bool test_algo1(Coords &coords, TAG &tag, string algo_name)
	{
		cout << "\ntest_algo1(...), with " << algo_name << endl;
		
		
		double min_fract_succ = 0.99, fraction_success;
		float d_l{5.0f}, d_h{12.0f};
	
		uint N = coords.size();
	
		dOut<1> out_serial(N, {d_l,d_h});
		find_distances(out_serial,coords,tag_simd<simd_no,float>());	
		vector<int> vout_serial = convert_dout_unique(out_serial);
		
		dOut<1> out(N, {d_l,d_h});
		find_distances(out,coords,tag);	
		vector<int> vout = convert_dout_unique(out);
		fraction_success = test_compare_douts(vout,vout_serial);
		bool success = fraction_success > min_fract_succ;
	
		if(success){
			cout << "test_algo(...): " << algo_name << " successful? " << boolalpha << true << endl;
			cout << "serial ncontacts=" << out_serial.counter[0] << " " << algo_name << " ncontacts=" << out.counter[0] << endl;
		}
		else{
			cout << "test_algo(...): " << algo_name << " did not pass. fraction_success=" << fraction_success << endl;
			cout << "serial ncontacts=" << out_serial.counter[0] << " " << algo_name << " ncontacts=" << out.counter[0] << endl;
			// for(uint k=0; k<out_serial.counter[0]; ++k){
			for(uint k=0; k<20; ++k){
				// if(int_high(vout_serial[k])!=0 && int_high(vout[k])!=0)
				// 	continue;
				cout << "k=" << k << " serial i and j: ";
				print_int_high(vout_serial[k]); print_int_low(vout_serial[k]);
				cout << " and " << algo_name << " i and j: ";
				print_int_high(vout[k]); print_int_low(vout[k]);
				cout << endl;
				// cout << "vout_serial=" << vout_serial[k] << endl;
				// cout << "vout=" << vout[k] << endl;
			}
		}

		return success;
	}
	
	template <typename TAG>
	bool test_algo2(Coords &coords, TAG &tag, string algo_name)
	{
		cout << "\ntest_algo2(...), with " << algo_name << endl;
		
		double min_fract_succ = 0.99, fraction_success;
	
		uint N = coords.size();
	
		dOut<2> out_serial(N, {5.0f, 12.0f, 6.0f, 13.0f});
		find_distances(out_serial,coords,tag_simd<simd_no,float>());
		dOut<2> out(N, {5.0f, 12.0f, 6.0f, 13.0f});
		find_distances(out,coords,tag);	
	
		bool success = true;	
		for(uint d=0; d<2; ++d){
			vector<int> vout_serial = convert_dout_unique(out_serial,d);
	
			vector<int> vout = convert_dout_unique(out,d);
			fraction_success = test_compare_douts(vout,vout_serial);
			success = fraction_success > min_fract_succ;
	
			if(success)
				cout << "test_algo2(...): " << algo_name << " d=" << d << " successful? " << boolalpha << true << endl;
			else{
				cout << "test_algo2(...): " << algo_name << " d=" << d << " did not pass. fraction_success=" << fraction_success << endl;
				cout << "serial ncontacts=" << out_serial.counter[0] << " " << algo_name << " ncontacts=" << out.counter[0] << endl;
				for(uint k=0; k<20; ++k){
					cout << "serial i and j: ";
					print_int_low(vout_serial[k]); print_int_high(vout_serial[k]);
					cout << " and " << algo_name << " i and j: ";
					print_int_low(vout[k]); print_int_high(vout[k]);
					cout << endl;
					// cout << "vout_serial=" << vout_serial[k] << endl;
					// cout << "vout=" << vout[k] << endl;
				}
			}
		}

		return success;
	}


	template <typename TAG>
	bool test_algo1_betwn_comps(Coords &c1, Coords &c2, TAG &tag, string algo_name)
	{
		cout << "\ntest_algo1_betwn_comps(...), with " << algo_name << endl;
		double min_fract_succ = 0.999, fraction_success;
		float d_l{5.0f}, d_h{12.0f};
	
		uint N1 = c1.size();
		uint N2 = c2.size();
		
		// uint N = std::max({N1,N2});
	
		dOut<1,false> out_serial(N1, N2, {d_l,d_h});
		find_distances(out_serial,c1,c2,tag_simd<simd_no,float>());	
		vector<int> vout_serial = convert_dout_unique(out_serial);
		
		dOut<1,false> out(N1, N2, {d_l,d_h});
		find_distances(out,c1,c2,tag);	
		vector<int> vout = convert_dout_unique(out);
		fraction_success = test_compare_douts(vout,vout_serial);
		bool success = fraction_success > min_fract_succ;
	
		if(success)
			cout << "test_algo1_betwn_comps(...): " << algo_name << " successful? " << boolalpha << true << endl;
		else{
			cout << "test_algo1_betwn_comps(...): " << algo_name << " did not pass. fraction_success=" << fraction_success << endl;
			cout << "serial ncontacts=" << out_serial.counter[0] << " " << algo_name << " ncontacts=" << out.counter[0] << endl;
			for(uint k=0; k<20; ++k){
				cout << "serial i and j: ";
				print_int_low(vout_serial[k]); print_int_high(vout_serial[k]);
				cout << " and " << algo_name << " i and j: ";
				print_int_low(vout[k]); print_int_high(vout[k]);
				cout << endl;
				// cout << "vout_serial=" << vout_serial[k] << endl;
				// cout << "vout=" << vout[k] << endl;
			}
		}

		return success;
	}
	
	
	template <typename TAG>
	bool test_algo2_betwn_comps(Coords &c1, Coords &c2, TAG &tag, string algo_name)
	{
		cout << "\ntest_algo2_betwn_comps(...), with " << algo_name << endl;
		
		double min_fract_succ = 0.999, fraction_success;
	
		uint N1 = c1.size();
		uint N2 = c2.size();		
	
		dOut<2,false> out_serial(N1, N2, {5.0f, 12.0f, 6.0f, 13.0f});
		find_distances(out_serial,c1,c2,tag_simd<simd_no,float>());
		dOut<2,false> out(N1, N2, {5.0f, 12.0f, 6.0f, 13.0f});
		find_distances(out,c1,c2,tag);	
	
		bool success = true;	
		for(uint d=0; d<2; ++d){
			vector<int> vout_serial = convert_dout_unique(out_serial,d);
	
			vector<int> vout = convert_dout_unique(out,d);
			fraction_success = test_compare_douts(vout,vout_serial);
			success = fraction_success > min_fract_succ;
	
			if(success){
				cout << "test_algo2_betwn_comps(...): " << algo_name << " d=" << d << " successful? " << boolalpha << true << endl;
				cout << "serial ncontacts=" << out_serial.counter[d] << " " << algo_name << " ncontacts=" << out.counter[d] << endl;
			}
			else{
				cout << "test_algo2_betwn_comps(...): " << algo_name << " d=" << d << " did not pass. fraction_success=" << fraction_success << endl;
				cout << "serial ncontacts=" << out_serial.counter[d] << " " << algo_name << " ncontacts=" << out.counter[d] << endl;
				for(uint k=0; k<20; ++k){
					cout << "serial i and j: ";
					print_int_low(vout_serial[k]); print_int_high(vout_serial[k]);
					cout << " and " << algo_name << " i and j: ";
					print_int_low(vout[k]); print_int_high(vout[k]);
					cout << endl;
					// cout << "vout_serial=" << vout_serial[k] << endl;
					// cout << "vout=" << vout[k] << endl;
				}
			}
		}

		return success;
	}
	

	bool test_dist_module(uint N1, uint N2)
	{
		bool success = true;
	
		Coords c1(N1), c2(N2);
		// scale_to_int16(c1);
	
		tag_simd<simd_no,  float>           t_serial;
		tag_simd<simd_avx, float>           t_avx;		
		tag_simd<dist::simd_avx_par,float>  t_avx_par;

#ifdef __CUDACC__
		tag_simd<cuda,     float>           t_cuda;
#endif
	
		cout << "\nTesting single comparisions" << endl;
		success =  success && test_algo1(c1, t_avx, "AVX");
		
		success =  success && test_algo1(c1, t_avx_par, "AVX-PARALLEL");


#ifdef __CUDACC__
		success =  success && test_algo1(c1, t_cuda, "CUDA");
#endif


		cout << "\nNow testing two comparisions" << endl;

		success =  success && test_algo2(c1, t_avx, "AVX");
		
		success =  success && test_algo2(c1, t_avx_par, "AVX-PARALLEL");

#ifdef __CUDACC__
		success =  success && test_algo2(c1, t_cuda, "CUDA");
#endif


		cout << "\nNow testing two-compartment functions:" << endl;

#ifdef __CUDACC__
		success =  success && success && test_algo1_betwn_comps(c1, c2, t_cuda, "CUDA");
#endif

		success =  success && test_algo1_betwn_comps(c1, c2, t_avx, "AVX");
		
		success =  success && test_algo1_betwn_comps(c1, c2, t_avx_par, "AVX-PARALLEL");
		

		cout << "\nNow testing two-compartment functions with two comparisions:" << endl;

		success =  success && test_algo2_betwn_comps(c1, c2, t_avx, "AVX");

		success =  success && test_algo2_betwn_comps(c1, c2, t_avx_par, "AVX-PARALLEL");
		

#ifdef __CUDACC__
		success =  success && test_algo2_betwn_comps(c1, c2, t_cuda, "CUDA");
#endif

		
		return success;
	}

} // end-of-namespace dist
