
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

/*
 Cylinder 1: x1------x2
 Cylinder 2: y1------y2


 x = x2-x1
 y = y2-y1
 z = x1-y1

 a = x.x;
 b = y.y;
 c = z.z;
 d = x.y;
 e = x.z;
 f = y.z;

 */

#include "CylinderExclVolRepulsion.h"
#include "CylinderExclVolRepulsionCUDA.h"
#include "CylinderExclVolume.h"

#include "Bead.h"
#include "Cylinder.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>
#include <CGMethod.h>

#ifdef CUDAACCL
#include "nvToolsExt.h"
#include <cuda.h>
#include <cuda_runtime.h>
#endif
typedef std::numeric_limits< floatingpoint > dbl;

using namespace mathfunc;
#ifdef CUDAACCL
//struct unaryfn: std::unary_function<size_t, unsigned long> {
//    int operator()(unsigned long i) const { return 12* i *sizeof(floatingpoint); }
//
//};
void CylinderExclVolRepulsion::deallocate(){
//    if(!(CUDAcommon::getCUDAvars().conservestreams))
//        CUDAcommon::handleerror(cudaStreamDestroy(stream),"cuda stream",
// "CylinderExclVolumeRepulsion.cu");
    CUDAcommon::handleerror(cudaFree(gU_i),"cudaFree", "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaFree(gU_sum),"cudaFree", "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaFree(gFF),"cudaFree", "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaFree(ginteraction),"cudaFree", "CylinderExclVolume.cu");
}
void CylinderExclVolRepulsion::optimalblocksnthreads( int nint, cudaStream_t stream_pass) {
//    //CUDA stream create
//    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
//        CUDAcommon::handleerror(cudaStreamCreate(&stream),"cuda stream", "CylinderExclVolume.cu");
    //
    stream = stream_pass;
    blocksnthreadse.clear();
    blocksnthreadsez.clear();
    blocksnthreadsf.clear();
    if(nint>0){
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the
    // maximum occupancy for a full device launch
    int gridSize;    // The actual grid size needed, based on input size
//    unaryfn::argument_type blksize;
//    unaryfn::result_type result;
//    unaryfn ufn;

    CUDAcommon::handleerror(cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                            CUDAExclVolRepulsionenergy, blockToSmem, 0),"cuda occupancy", "CylinderExclVolume.cu");
//    std::cout<<(nint +blockSize -1) / blockSize<<" "<<blockSize<<endl;
//
//    cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,
//                                        CUDAExclVolRepulsionenergy, 0, 0);
    blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
    blocksnthreadse.push_back(blockSize);
//    std::cout<<(nint +blockSize -1) / blockSize<<" "<<blockSize<<endl;
    blockSize = 0;

        CUDAcommon::handleerror(cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                CUDAExclVolRepulsionenergyz, blockToSmem, 0),"cuda occupancy", "CylinderExclVolume.cu");
    blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
    blocksnthreadsez.push_back(blockSize);
    blockSize = 0;

        CUDAcommon::handleerror(cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                CUDAExclVolRepulsionforce, blockToSmem, 0),"cuda occupancy", "CylinderExclVolume.cu");
    blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
    blocksnthreadsf.push_back(blockSize);
//get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMemsetAsync(gU_i, 0, bntaddvec2.at(0) * sizeof(floatingpoint), stream),
                                "cuda data transfer",
                                "CylinderExclVolume.cu");
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(floatingpoint)),"cuda data transfer",
                                "CylinderExclVolume.cu");
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)),"cuda data transfer",
                                "CylinderExclVolume.cu");
        char a[] = "Excluded Volume";
        char b[] =  "Cylinder Excluded Volume";
        CUDAcommon::handleerror(cudaMalloc((void **) &gFF, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMalloc((void **) &ginteraction, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMemcpyAsync(gFF, a, 100 * sizeof(char),
                                            cudaMemcpyHostToDevice, stream));
        CUDAcommon::handleerror(cudaMemcpyAsync(ginteraction, b, 100 * sizeof(char),
                                            cudaMemcpyHostToDevice, stream));


}
    else{
        blocksnthreadse.push_back(0);
        blocksnthreadse.push_back(0);
        blocksnthreadsez.push_back(0);
        blocksnthreadsez.push_back(0);
        blocksnthreadsf.push_back(0);
        blocksnthreadsf.push_back(0);
    }

}
floatingpoint* CylinderExclVolRepulsion::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                        floatingpoint *krep, int *params) {
//    if(blocksnthreadse[1]>0) {
//

//        CUDAExclVolRepulsionenergy << < blocksnthreadse[0], blocksnthreadse[1],
//                (12 * blocksnthreadse[1]) * sizeof(floatingpoint), stream >> >(coord, f, beadSet, krep, params, gU_i,
//                CUDAcommon::getCUDAvars().gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
////        CUDAcommon::handleerror(cudaEventRecord(event, stream));


//        CUDAcommon::handleerror(cudaGetLastError(),"CUDAExclVolRepulsionenergy", "CylinderExclVolumeRepulsion.cu");

//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0, stream>>>(gU_i,params, gU_sum, gpu_Utot);

//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;

//        CUDAcommon::handleerror( cudaGetLastError() ,"CUDAExclVolRepulsionenergy", "CylinderExclVolumeRepulsion.cu");

//
//        return gU_sum;
//    }
//    else return NULL;
}

floatingpoint* CylinderExclVolRepulsion::energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, floatingpoint *z, int *params) {

//    if(blocksnthreadse[1]>0) {
//        CUDAExclVolRepulsionenergy << < blocksnthreadse[0], blocksnthreadse[1],
//                (12 * blocksnthreadse[1]) * sizeof(floatingpoint), stream >> >(coord, f, beadSet, krep, params, gU_i, z,
//                CUDAcommon::getCUDAvars().gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        CUDAcommon::handleerror(cudaGetLastError(),"CUDAExclVolRepulsionenergy", "CylinderExclVolumeRepulsion.cu");
//    }
    if(blocksnthreadsez[1]>0) {
        auto boolvarvec = CUDAcommon::cudavars.backtrackbools;
        CUDAExclVolRepulsionenergyz << < blocksnthreadsez[0], blocksnthreadsez[1],
                12 * blocksnthreadsez[1] * sizeof(floatingpoint),stream >> > (coord, f, beadSet,
                krep, params, gU_i, CUDAcommon::cudavars.gpu_energyvec, z,
                CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction, boolvarvec.at(0),
                boolvarvec.at(1));
    }
    if(blocksnthreadse[1]<=0 && blocksnthreadsez[1]<=0)
        return NULL;
    else{
#ifdef CUDA_INDIVIDUAL_ESUM
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
        resetfloatingpointvariableCUDA<<<1,1,0,stream>>>(gU_sum);
        addvectorred2<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint),stream>>>(gU_i,
                params, gU_sum, gpu_Utot);
#endif
                CUDAcommon::handleerror( cudaGetLastError() ,"CUDAExclVolRepulsionenergy",
"CylinderExclVolumeRepulsion.cu");
        return gU_sum;
    }
}

void CylinderExclVolRepulsion::forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, int *params) {
//    cudaEvent_t start, stop;
//    CUDAcommon::handleerror(cudaEventCreate( &start));
//    CUDAcommon::handleerror(cudaEventCreate( &stop));
//    CUDAcommon::handleerror(cudaEventRecord( start, 0));

    if(blocksnthreadsf[1]>0) {
//        floatingpoint *gU_ii;
//        floatingpoint *gf1, *gf2, *gf3, *gf4, *gf5;
//        floatingpoint *gc1, *gc2, *gcheckU;
//        floatingpoint U_ii[blocksnthreadsf[0] * blocksnthreadsf[1]];
//        floatingpoint c1[3 * blocksnthreadsf[0] * blocksnthreadsf[1]], c2[3 * blocksnthreadsf[0] * blocksnthreadsf[1]];
//        floatingpoint F_i[3 * blocksnthreadsf[0] * blocksnthreadsf[1]];
//        floatingpoint checkU[blocksnthreadsf[1]];

//        std::cout << "CEVF Number of Blocks: " << blocksnthreadsf[0] << endl;
//        std::cout << "Threads per block: " << blocksnthreadsf[1] << endl;

        //TODO  since the number of threads needed is constant through out the minimization, consider storing the pointer.
//        CUDAcommon::handleerror(cudaMalloc((void **) &gf1, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof
// (floatingpoint)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gf2, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof
// (floatingpoint)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gf3, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof
// (floatingpoint)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gf4, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof
// (floatingpoint)));
//        CUDAcommon::handleerror(
//                cudaMalloc((void **) &gf5, 45 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof(floatingpoint)));

//            size_t freeMem, totalMem;
//
//    cudaMemGetInfo(&freeMem, &totalMem);
//
//    std::cout<<"Memory "<<freeMem<<" "<<totalMem<<endl;
//        struct cudaDeviceProp properties;
//        cudaGetDeviceProperties(&properties, 0);
//        cout << "using " << properties.multiProcessorCount << " multiprocessors" << endl;
//        cout << "max threads per processor: " << properties.maxThreadsPerMultiProcessor << endl;
//        std::cout << 12 *blocksnthreadsf[1] * sizeof(floatingpoint) << endl;
//        int blockSize;   // The launch configurator returned block size
//        int minGridSize; // The minimum grid size needed to achieve the
//        // maximum occupancy for a full device launch
//        int gridSize;    // The actual grid size needed, based on input size
//        cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,
//                                            CUDAExclVolRepulsionforce, 0, 0);
//        gridSize = (blocksnthreadsf[0] * blocksnthreadsf[1] + blockSize -1) / blockSize;
//
//        std::cout<<gridSize<<" "<<blockSize<<endl;
//        size_t my_kernel_sm_size, count;

//        cudaOccupancyMaxPotentialBlockSizeVariableSMem( &minGridSize, &blockSize, CUDAExclVolRepulsionforce,
//                                                        my_kernel_sm_size, count);
//        std::cout<<minGridSize<<" "<<blockSize<<" "<<my_kernel_sm_size<<" "<<count<<endl;
//        int numblocks;
//        unsigned int flag;
//        cudaOccupancyMaxActiveBlocksPerMultiprocessorWithFlags(&numblocks,  CUDAExclVolRepulsionforce,blockSize,
//                                                               my_kernel_sm_size,flag);
//        std::cout<<numblocks<<" "<<blockSize<<" "<<my_kernel_sm_size<<endl;

        CUDAExclVolRepulsionforce << < blocksnthreadsf[0],blocksnthreadsf[1],
                12 *blocksnthreadsf[1] * sizeof(floatingpoint),stream >> >
                (coord, f, beadSet, krep, params);

        CUDAcommon::handleerror(cudaGetLastError(),"CUDAExclVolRepulsionforce", "CylinderExclVolumeRepulsion.cu");
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
//    CUDAcommon::handleerror( cudaPeekAtLastError() );
       //CUDAcommon::handleerror(cudaDeviceSynchronize());

//    floatingpoint f1 [3 * blocksnthreadsf[0]*blocksnthreadsf[1]];
//    floatingpoint f2 [3 * blocksnthreadsf[0]*blocksnthreadsf[1]];
//    floatingpoint f3 [3 * blocksnthreadsf[0]*blocksnthreadsf[1]];
//    floatingpoint f4 [3 * blocksnthreadsf[0]*blocksnthreadsf[1]];
//    floatingpoint f5 [45 * blocksnthreadsf[0]*blocksnthreadsf[1]];
//    cudaMemcpy(f1, gf1, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof(floatingpoint), cudaMemcpyDeviceToHost);
//    cudaMemcpy(f2, gf2, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof(floatingpoint), cudaMemcpyDeviceToHost);
//    cudaMemcpy(f3, gf3, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof(floatingpoint), cudaMemcpyDeviceToHost);
//    cudaMemcpy(f4, gf4, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof(floatingpoint), cudaMemcpyDeviceToHost);
//    cudaMemcpy(f5, gf5, 45 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof(floatingpoint), cudaMemcpyDeviceToHost);
//
//    for(auto i=0;i < blocksnthreadsf[0] * blocksnthreadsf[1]; i++) {
//        for(auto j=0;j<44;j++)
//            std::cout<<f5[44*i+j]<<" ";
//        std::cout << f1[3 * i] << " " << f1[3 * i + 1] << " " << f1[3 * i + 2] << " " << f2[3 * i] << " "
//                  << f2[3 * i + 1] << " " << f2[3 * i + 2] << " "
//                  << f3[3 * i] << " " << f3[3 * i + 1] << " " << f3[3 * i + 2] << " " << f4[3 * i] << " "
//                  << f4[3 * i + 1] << " " << f4[3 * i + 2] << endl;
//    }

//        CUDAcommon::handleerror(cudaFree(gf1));
//        CUDAcommon::handleerror(cudaFree(gf2));
//        CUDAcommon::handleerror(cudaFree(gf3));
//        CUDAcommon::handleerror(cudaFree(gf4));
//        CUDAcommon::handleerror(cudaFree(gf5));
    }
//    CUDAcommon::handleerror(cudaEventRecord( stop, 0));
//    CUDAcommon::handleerror(cudaEventSynchronize(stop));
//    float elaspedtime;
//    CUDAcommon::handleerror(cudaEventElapsedTime(&elaspedtime, start, stop));
//    CUDAvars cvars=CUDAcommon::getCUDAvars();
//    cvars.Ccforce += elaspedtime;
//    std::cout<<"C CFE "<<elaspedtime<<endl;
//    CUDAcommon::cudavars=cvars;
//    CUDAcommon::handleerror(cudaEventDestroy(start));
//    CUDAcommon::handleerror(cudaEventDestroy(stop));
}

void CylinderExclVolRepulsion::checkforculprit() {
    CUDAcommon::printculprit("Excluded Volume", "Cylinder Excluded Volume");
    int i = 0;
    cout<<"Printing culprit cylinders.."<<endl;
    for (auto cyl: Cylinder::getCylinders()) {
            auto id1 = cyl->getFirstBead()->_dbIndex;
            auto id2 = cyl->getSecondBead()->_dbIndex;
            if(id1 == CUDAcommon::getCUDAvars().culpritID[0] && id2 == CUDAcommon::getCUDAvars().culpritID[1])
                cyl->printSelf();
            else if(id1 == CUDAcommon::getCUDAvars().culpritID[2] && id2 == CUDAcommon::getCUDAvars().culpritID[3])
                cyl->printSelf();
        }
    exit(EXIT_FAILURE);
}

#endif
floatingpoint CylinderExclVolRepulsion::energy(floatingpoint *coord, floatingpoint *force,
                                               int *beadSet, floatingpoint *krep) {
	floatingpoint *c1, *c2, *c3, *c2temp, *c4, *newc2, d;

	doubleprecision a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
	doubleprecision ATG1, ATG2, ATG3, ATG4;

	//Additional paramteres
	doubleprecision sqmag_D, sqmag_E;
	doubleprecision g,h,I;
	doubleprecision *cp = new doubleprecision[3];
	floatingpoint *vec_A = new floatingpoint[3];
	floatingpoint *vec_B = new floatingpoint[3];
	floatingpoint *vec_C = new floatingpoint[3];

	int nint = CylinderExclVolume<CylinderExclVolRepulsion>::numInteractions;
	int n = CylinderExclVolume<CylinderExclVolRepulsion>::n;

	floatingpoint U_i = 0.0;
	floatingpoint U = 0.0;
	newc2 = new floatingpoint[3];
//    std::cout<<"SERL ecvol nint "<<nint<<endl;
	for (int i = 0; i < nint; i++) {

		c1 = &coord[3 * beadSet[n * i]];
		c2 = &coord[3 * beadSet[n * i + 1]];
		c2temp = &coord[3 * beadSet[n * i + 1]];
		c3 = &coord[3 * beadSet[n * i + 2]];
		c4 = &coord[3 * beadSet[n * i + 3]];

		floatingpoint trialvec[4]={0.01, 0.1, 1.0, 2.7};
		uint trial = 0;
//
//        for(uint trial = 0; trial < 5; trial++) {
		//trial = 0 signifies usual calculaton. trial = 1 & 2 signify calculation
		// after moving a bead. Trial 2 will decide if the calculation has to end
		// with error.

		//Calculate energy
		//@{
		a = scalarProduct(c1, c2, c1, c2);//always positive
		b = scalarProduct(c3, c4, c3, c4);//always positive
		c = scalarProduct(c3, c1, c3, c1);//always positive
		//If c1, c2, c3, and c4 are coplanar, e = d + f
		d = scalarProduct(c1, c2, c3, c4);
		e = scalarProduct(c1, c2, c3, c1);
		F = scalarProduct(c3, c4, c3, c1);

		//Additional parameters added to minimize loss of signifiance calculations.
		g = scalarProduct(c1, c2, c1, c4);
		h = scalarProduct(c3, c4, c3, c2);
		I = scalarProduct(c3, c4, c1, c4);

		//Preprocessing required to calculate AA through JJ.
		//@{
		sqmag_D = sqmagnitude(c1, c4);
		sqmag_E = sqmagnitude(c3, c2);

		for(auto dim = 0; dim < 3 ; dim++) {
			vec_A[dim] = c1[dim] - c2[dim];
			vec_B[dim] = c3[dim] - c4[dim];
			vec_C[dim] = c3[dim] - c1[dim];
		}
		//@}

		doubleprecision ac = sqrt(a*c);
		doubleprecision bc = sqrt(b*c);
		doubleprecision aD = sqrt(a*sqmag_D);
		doubleprecision bE = sqrt(b*sqmag_E);

		AA = sqrt(max(0.0,(ac + e)*(ac - e)));
		BB = sqrt(max(0.0,(bc + F)*(bc - F)));

		CC = d*e - a*F;
		DD = b*e - d*F;

		EE = sqrt(max(0.0,(aD + g)*(aD - g)));
		FF = sqrt(max(0.0,(bE + h)*(bE - h)));

		GG = d*g - a*I;
		HH = CC + GG - DD;
		crossProduct(cp, vec_A, vec_B);
		JJ = scalarProduct(cp,vec_C);
		JJ = - JJ*JJ;

//            cout<<AA<<" "<<BB<<" "<<CC<<" "<<DD<<" "<<EE<<" "<<FF<<" "<<GG<<" "<<HH<<" "<<JJ<<endl;

		/*AA = sqrt(a*c - e*e);//always positive
		BB = sqrt(b*c - F*F);//always positive

		CC = d*e - a*F;
		DD = b*e - d*F;

		EE = sqrt( a*(b + c - 2*F) - (d - e)*(d - e) );
		FF = sqrt( b*(a + c + 2*e) - (d + F)*(d + F) );

		GG = d*d - a*b - CC;
		HH = CC + GG - DD;
		JJ = c*(GG + CC) + e*DD - F*CC;*/

//            cout<<AA<<" "<<BB<<" "<<CC<<" "<<DD<<" "<<EE<<" "<<FF<<" "<<GG<<" "<<HH<<" "<<JJ<<endl;

		ATG1 = atan( (a + e)/AA) - atan(e/AA);
		ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
		ATG3 = atan((F)/BB) - atan((F - b)/BB);
		ATG4 = atan((d + F)/FF) - atan((d + F - b)/FF);

		U_i = 0.5 * krep[i]/ JJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);

		//@@}
		/* cout<<"Analytical energy "<<U_i<<endl;
		 doubleprecision temp = 	energyN(coord, force, beadSet, krep, i);*/

		//Check if energy is acceptable
		if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
		   || U_i != U_i || U_i < -1.0) {
			//evaluate numerically
			U_i = energyN(coord, force, beadSet, krep, i);
			if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
			   || U_i != U_i || U_i < -1.0) {
				short found = 0;
				for (auto cyl:Cylinder::getCylinders()) {
					auto dbIndex1 = cyl->getFirstBead()->_dbIndex;
					auto dbIndex2 = cyl->getSecondBead()->_dbIndex;
					if (dbIndex1 == beadSet[n * i] &&
					    dbIndex2 == beadSet[n * i + 1]) {
						CylinderVolumeInteractions::_cylinderCulprit1 = cyl;
						found++;
						if (found >= 2)
							break;
					} else if (dbIndex1 == beadSet[n * i + 2] &&
					           dbIndex2 == beadSet[n * i +
					                               3]) {
						CylinderVolumeInteractions::_cylinderCulprit2 = cyl;
						found++;
						if (found >= 2)
							break;
					}
				}
				cout << "Printing relevant coordinates " << endl;
				cout << "c1 " << c1[0] << " " << c1[1] << " " << c1[2] << endl;
				cout << "c2 " << c2[0] << " " << c2[1] << " " << c2[2] << endl;
				cout << "c3 " << c3[0] << " " << c3[1] << " " << c3[2] << endl;
				cout << "c4 " << c4[0] << " " << c4[1] << " " << c4[2] << endl;

				cout << "Printing infinite energy contributions " << endl;
				cout << "a " << a << " b " << b << " c " << c << " d " << d << " e "
				     << e << " f " << F << endl;
				cout << "JJ " << JJ << " CC " << CC << " AA " << AA << " ATG1 "
				     << ATG1 << " GG " << GG << " EE "
				     << EE << " ATG2 " << ATG2 << " DD " << DD << " BB " << BB
				     << " ATG3 " << ATG3 << " HH " << HH << " FF "
				     << FF << " ATG4 " << ATG4 << " U_i " << U_i << endl;
				AA = sqrt(a * c - e * e);
				BB = sqrt(b * c - F * F);

				CC = d * e - a * F;
				DD = b * e - d * F;

				EE = sqrt(a * (b + c - 2 * F) - (d - e) * (d - e));
				FF = sqrt(b * (a + c + 2 * e) - (d + F) * (d + F));

				GG = d * d - a * b - CC;
				HH = CC + GG - DD;
				JJ = c * (GG + CC) + e * DD - F * CC;
				cout << "JJ " << JJ << " CC " << CC << " AA " << AA << " ATG1 "
				     << ATG1 << " GG " << GG << " EE "
				     << EE << " ATG2 " << ATG2 << " DD " << DD << " BB " << BB
				     << " ATG3 " << ATG3 << " HH " << HH << " FF "
				     << FF << " ATG4 " << ATG4 << endl;
				return -1;
			}
			else {
				U += U_i;
				continue;
			}

			if (false) {

				//if in same plane, recalculate after moving beads.
				if (trial < 4) {
					auto xx = areInPlane(c1, c2, c3, c4);
					cout << "areinPlane " << xx << endl;
					if (xx) {
						cout << "moving out of plane" << endl;
						//slightly move point
						movePointOutOfPlane(c1, c2temp, c3, c4, newc2, 2,
						                    trialvec[trial]);
						// SysParams::exvolcounter[1] += 1;
						c2 = newc2;
					}
				}
					//If all trials are done, set Culprit and return
				else {
					short found = 0;
					for (auto cyl:Cylinder::getCylinders()) {
						auto dbIndex1 = cyl->getFirstBead()->_dbIndex;
						auto dbIndex2 = cyl->getSecondBead()->_dbIndex;
						if (dbIndex1 == beadSet[n * i] &&
						    dbIndex2 == beadSet[n * i + 1]) {
							CylinderVolumeInteractions::_cylinderCulprit1 = cyl;
							found++;
							if (found >= 2)
								break;
						} else if (dbIndex1 == beadSet[n * i + 2] &&
						           dbIndex2 == beadSet[n * i +
						                               3]) {
							CylinderVolumeInteractions::_cylinderCulprit2 = cyl;
							found++;
							if (found >= 2)
								break;
						}
					}
					cout << "Printing relevant coordinates " << endl;
					cout << "c1 " << c1[0] << " " << c1[1] << " " << c1[2] << endl;
					cout << "c2 " << c2[0] << " " << c2[1] << " " << c2[2] << endl;
					cout << "c3 " << c3[0] << " " << c3[1] << " " << c3[2] << endl;
					cout << "c4 " << c4[0] << " " << c4[1] << " " << c4[2] << endl;

					cout << "Printing infinite energy contributions " << endl;
					cout << "a " << a << " b " << b << " c " << c << " d " << d << " e "
					     << e << " f " << F << endl;
					cout << "JJ " << JJ << " CC " << CC << " AA " << AA << " ATG1 "
					     << ATG1 << " GG " << GG << " EE "
					     << EE << " ATG2 " << ATG2 << " DD " << DD << " BB " << BB
					     << " ATG3 " << ATG3 << " HH " << HH << " FF "
					     << FF << " ATG4 " << ATG4 << " U_i " << U_i << endl;
					AA = sqrt(a * c - e * e);
					BB = sqrt(b * c - F * F);

					CC = d * e - a * F;
					DD = b * e - d * F;

					EE = sqrt(a * (b + c - 2 * F) - (d - e) * (d - e));
					FF = sqrt(b * (a + c + 2 * e) - (d + F) * (d + F));

					GG = d * d - a * b - CC;
					HH = CC + GG - DD;
					JJ = c * (GG + CC) + e * DD - F * CC;
					cout << "JJ " << JJ << " CC " << CC << " AA " << AA << " ATG1 "
					     << ATG1 << " GG " << GG << " EE "
					     << EE << " ATG2 " << ATG2 << " DD " << DD << " BB " << BB
					     << " ATG3 " << ATG3 << " HH " << HH << " FF "
					     << FF << " ATG4 " << ATG4 << endl;
					return -1;
				}
			}
		}
			//add energy to total energy and move on to the next interaction.
		else
			U += U_i;

//        }
	}
	delete [] newc2;
	delete [] cp;
	delete [] vec_A;
	delete [] vec_B;
	delete [] vec_C;
	return U;
}


floatingpoint CylinderExclVolRepulsion::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                               floatingpoint *krep, floatingpoint z) {


	floatingpoint d;
	doubleprecision invDSquare;
	floatingpoint *f1, *f2, *f3, *f4;
	doubleprecision a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
	doubleprecision ATG1, ATG2, ATG3, ATG4;
	floatingpoint *c1us, *c2us, *c3us, *c4us;
	floatingpoint *c1 = new floatingpoint[3];//stretched
	floatingpoint *c2 = new floatingpoint[3];
	floatingpoint *c3 = new floatingpoint[3];
	floatingpoint *c4 = new floatingpoint[3];
	floatingpoint *newc2 = new floatingpoint[3];

	int nint = CylinderExclVolume<CylinderExclVolRepulsion>::numInteractions;
	int n = CylinderExclVolume<CylinderExclVolRepulsion>::n;

	floatingpoint U_i = 0.0;
	floatingpoint U = 0.0;

	for (int i = 0; i < nint; i++) {

		c1us = &coord[3 * beadSet[n * i]];
		c2us = &coord[3 * beadSet[n * i +1]];
		c3us = &coord[3 * beadSet[n * i +2]];
		c4us = &coord[3 * beadSet[n * i +3]];

		//stretch coords
		f1 = &f[3 * beadSet[n * i]];
		f2 = &f[3 * beadSet[n * i + 1]];
		f3 = &f[3 * beadSet[n * i + 2]];
		f4 = &f[3 * beadSet[n * i + 3]];

		c1[0] = c1us[0] + z * f1[0];
		c1[1] = c1us[1] + z * f1[1];
		c1[2] = c1us[2] + z * f1[2];
		c2[0] = c2us[0] + z * f2[0];
		c2[1] = c2us[1] + z * f2[1];
		c2[2] = c2us[2] + z * f2[2];
		c3[0] = c3us[0] + z * f3[0];
		c3[1] = c3us[1] + z * f3[1];
		c3[2] = c3us[2] + z * f3[2];
		c4[0] = c4us[0] + z * f4[0];
		c4[1] = c4us[1] + z * f4[1];
		c4[2] = c4us[2] + z * f4[2];

		//check if parallel
		if(areParallel(c1, c2, c3, c4)) {
//            SysParams::exvolcounterz[0] += 1;
			d = twoPointDistance(c1, c3);
			invDSquare =  1 / (d * d);
			U_i = krep[i] * invDSquare * invDSquare;
//            std::cout<<"P Energy"<<U_i<<endl;
			if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
			   || U_i != U_i || U_i < -1.0) {
				short found = 0;
				for(auto cyl:Cylinder::getCylinders()){
					auto dbIndex1 = cyl->getFirstBead()->_dbIndex;
					auto dbIndex2 = cyl->getSecondBead()->_dbIndex;
					if(dbIndex1 == beadSet[n * i] && dbIndex2 == beadSet[n * i + 1]) {
						CylinderVolumeInteractions::_cylinderCulprit1 = cyl;
						found++;
						if(found>=2)
							break;
					}
					else if(dbIndex1 == beadSet[n * i + 2] && dbIndex2 == beadSet[n * i +
					                                                              3]){
						CylinderVolumeInteractions::_cylinderCulprit2 = cyl;
						found++;
						if(found>=2)
							break;
					}
				}
				return -1;
			}
			U += U_i;
			continue;
		}

		//check if in same plane
		if(areInPlane(c1, c2, c3, c4)) {
//            SysParams::exvolcounterz[1] += 1;
			//slightly move point
			movePointOutOfPlane(c1, c2, c3, c4, newc2, 2, 0.01);
			c2 = newc2;
//            std::cout<<"move"<<endl;
//            std::cout<<i<<" 2.0 "<<c1[0]<<" "<<c1[1]<<" "<<c1[2]<<" "<<c2[0]<<" "<<c2[1]<<" "<<c2[2]<<" "<<c3[0]<<" "
//                    ""<<c3[1]<<" "
//                             ""<<c3[2]<<" "<<c4[0]<<" "<<c4[1]<<" "<<c4[2]<<" "<<U_i<<endl;
		}


		a = scalarProduct(c1, c2, c1, c2);
		b = scalarProduct(c3, c4, c3, c4);
		c = scalarProduct(c3, c1, c3, c1);
		d = scalarProduct(c1, c2, c3, c4);
		e = scalarProduct(c1, c2, c3, c1);
		F = scalarProduct(c3, c4, c3, c1);
//    std::cout<<i<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<F<<endl;
		AA = sqrt(a*c - e*e);
		BB = sqrt(b*c - F*F);

		CC = d*e - a*F;
		DD = b*e - d*F;

		EE = sqrt( a*(b + c - 2*F) - (d - e)*(d - e) );
		FF = sqrt( b*(a + c + 2*e) - (d + F)*(d + F) );

		GG = d*d - a*b - CC;
		HH = CC + GG - DD;
		JJ = c*(GG + CC) + e*DD - F*CC;


		ATG1 = atan( (a + e)/AA) - atan(e/AA);
		ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
		ATG3 = atan((F)/BB) - atan((F - b)/BB);
		ATG4 = atan((d + F)/FF) - atan((d + F - b)/FF);

		U_i = 0.5 * krep[i]/ JJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);

//        std::cout<<"analytical energy "<<U_i<<endl;

		if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
		   || U_i != U_i || U_i < -1.0) {
			short found = 0;
			for(auto cyl:Cylinder::getCylinders()){
				auto dbIndex1 = cyl->getFirstBead()->_dbIndex;
				auto dbIndex2 = cyl->getSecondBead()->_dbIndex;
				if(dbIndex1 == beadSet[n * i] && dbIndex2 == beadSet[n * i + 1]) {
					CylinderVolumeInteractions::_cylinderCulprit1 = cyl;
					found++;
					if(found>=2)
						break;
				}
				else if(dbIndex1 == beadSet[n * i + 2] && dbIndex2 == beadSet[n * i +
				                                                              3]){
					CylinderVolumeInteractions::_cylinderCulprit2 = cyl;
					found++;
					if(found>=2)
						break;
				}
			}

			return -1.0;
		}
/*        std::cout<<i<<" "<<U_i<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<F<<" "<<AA<<" "<<BB<<" "<<CC<<" "<<DD<<" "
                ""<<EE<<" "
                ""<<FF<<""
                " "<<GG<<" "<<HH<<" "<<JJ<<" "<<ATG1<<" "<<ATG2<<" "<<ATG3<<" "<<ATG4<<endl;
        std::cout<<"Energy "<<U_i<<endl;*/
		U += U_i;
	}
//    std::cout<<"Total energy serial "<<U<<endl;
	delete [] c1;
	delete [] c2;
	delete [] c3;
	delete [] c4;
	delete [] newc2;
//    delete [] c1us;
//    delete [] c2us;
//    delete [] c3us;
//    delete [] c4us;
	return U;
}

void CylinderExclVolRepulsion::forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep) {

	floatingpoint *c1, *c2, *c3, *c4, *newc2, d, U;
	floatingpoint *f1, *f2, *f3, *f4;

	doubleprecision a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ, invJJ;
	doubleprecision ATG1, ATG2, ATG3, ATG4;
	doubleprecision A1, A2, E1, E2, B1, B2, F1, F2, A11, A12, A13, A14;
	doubleprecision E11, E12, E13, E14, B11, B12, B13, B14, F11, F12, F13, F14;

	//Additional paramteres
	doubleprecision sqmag_D, sqmag_E;
	doubleprecision g,h,I;
	doubleprecision *cp = new doubleprecision[3];
	floatingpoint *vec_A = new floatingpoint[3];
	floatingpoint *vec_B = new floatingpoint[3];
	floatingpoint *vec_C = new floatingpoint[3];

	int nint = CylinderExclVolume<CylinderExclVolRepulsion>::numInteractions;
	int n = CylinderExclVolume<CylinderExclVolRepulsion>::n;

	newc2 = new floatingpoint[3];
	for (int i = 0; i < nint; i++) {

		c1 = &coord[3 * beadSet[n * i]];
		c2 = &coord[3 * beadSet[n * i + 1]];
		c3 = &coord[3 * beadSet[n * i + 2]];
		c4 = &coord[3 * beadSet[n * i + 3]];

		//stretch coords
		f1 = &f[3 * beadSet[n * i]];
		f2 = &f[3 * beadSet[n * i + 1]];
		f3 = &f[3 * beadSet[n * i + 2]];
		f4 = &f[3 * beadSet[n * i + 3]];

		//check if in same plane
		if(areInPlane(c1, c2, c3, c4)) {

			//slightly move point
			movePointOutOfPlane(c1, c2, c3, c4, newc2, 2, 0.01);
			c2 = newc2;

#ifdef DETAILEDOUTPUT
			std::cout<<"Mv"<<c1[0]<<" "<<c1[1]<<" "<<c1[2]<<" "<<
                     c2[0]<<" "<<c2[1]<<" "<<c2[2]<<" "<<
                     c3[0]<<" "<<c3[1]<<" "<<c3[2]<<" "<<
                     c4[0]<<" "<<c4[1]<<" "<<c4[2]<<endl;
            std::cout<<"M ";
#endif
		}
#ifdef DETAILEDOUTPUT
		else{
            std::cout<<"N ";
        }
#endif
		a = scalarProduct(c1, c2, c1, c2);
		b = scalarProduct(c3, c4, c3, c4);
		c = scalarProduct(c3, c1, c3, c1);
		d = scalarProduct(c1, c2, c3, c4);
		e = scalarProduct(c1, c2, c3, c1);
		F = scalarProduct(c3, c4, c3, c1);

		//Additional parameters added to minimize loss of signifiance calculations.
		g = scalarProduct(c1, c2, c1, c4);
		h = scalarProduct(c3, c4, c3, c2);
		I = scalarProduct(c3, c4, c1, c4);

		//Preprocessing required to calculate AA through JJ.
		//@{
		sqmag_D = sqmagnitude(c1, c4);
		sqmag_E = sqmagnitude(c3, c2);

		for(auto dim = 0; dim < 3 ; dim++) {
			vec_A[dim] = c1[dim] - c2[dim];
			vec_B[dim] = c3[dim] - c4[dim];
			vec_C[dim] = c3[dim] - c1[dim];
		}
		//@}

		doubleprecision ac = sqrt(a*c);
		doubleprecision bc = sqrt(b*c);
		doubleprecision aD = sqrt(a*sqmag_D);
		doubleprecision bE = sqrt(b*sqmag_E);

		AA = sqrt(max(0.0,(ac + e)*(ac - e)));
		BB = sqrt(max(0.0,(bc + F)*(bc - F)));

		CC = d*e - a*F;
		DD = b*e - d*F;

		EE = sqrt(max(0.0,(aD + g)*(aD - g)));
		FF = sqrt(max(0.0,(bE + h)*(bE - h)));

		GG = d*g - a*I;
		HH = CC + GG - DD;
		crossProduct(cp, vec_A, vec_B);
		JJ = scalarProduct(cp,vec_C);
		JJ = - JJ*JJ;
		invJJ = 1/JJ;

//        std::cout<<"N "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<F<<endl;
		/* AA = sqrt(a*c - e*e);
		 BB = sqrt(b*c - F*F);

		 CC = d*e - a*F;
		 DD = b*e - d*F;

		 EE = sqrt( a*(b + c - 2*F) - (d - e)*(d - e) );
		 FF = sqrt( b*(a + c + 2*e) - (d + F)*(d + F) );

		 GG = d*d - a*b - CC;
		 HH = CC + GG - DD;
		 JJ = c*(GG + CC) + e*DD - F*CC;*/
//        std::cout<<"N2 "<<AA<<" "<<BB<<" "<<CC<<" "<<DD<<" "<<EE<<" "<<FF<<" "<<GG<<" "<<HH<<" "<<JJ<<endl;


		ATG1 = atan( (a + e)/AA) - atan(e/AA);
		ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
		ATG3 = atan((F)/BB) - atan((F - b)/BB);
		ATG4 = atan((d + F)/FF) - atan((d + F - b)/FF);

#ifdef DETAILEDOUTPUT
		std::cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<F<<" "<<AA<<" "<<BB<<" "<<CC<<" "
                ""<<DD<<" "<<EE<<" "<<FF<<" "<<GG<<" "<<HH<<" "<<JJ<<" "<<ATG1<<" "<<ATG2<<" "
                         ""<<ATG3<<" "<<ATG4<<" "<<U<<" "<<krep[i]<<endl;
#endif
		U = 0.5 * krep[i]*invJJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);

		A1 = AA*AA/(AA*AA + e*e);
		A2 = AA*AA/(AA*AA + (a + e)*(a + e));

		E1 = EE*EE/(EE*EE + (a + e - d)*(a + e - d));
		E2 = EE*EE/(EE*EE + (e - d)*(e - d));

		B1 = BB*BB/(BB*BB + (F - b)*(F - b));
		B2 = BB*BB/(BB*BB + F*F);

		F1 = FF*FF/(FF*FF + (d + F - b)*(d + F - b));
		F2 = FF*FF/(FF*FF + (d + F)*(d + F));

		A11 = ATG1/AA;
		A12 = -((ATG1*CC)/(AA*AA)) + (A1*CC*e)/(AA*AA*AA) -
		      (A2*CC*(a + e))/(AA*AA*AA);
		A13 = -((A1*CC)/(AA*AA)) + (A2*CC)/(AA*AA);
		A14 = (A2*CC)/(AA*AA);

		E11 = ATG2/EE;
		E12 = (E2*(-a + d - e)*GG)/(EE*EE*EE) + (E1*(-d + e)*GG)/(EE*EE*EE) -
		      (ATG2*GG)/(EE*EE);
		E13 = -((E1*GG)/(EE*EE)) + (E2*GG)/(EE*EE);
		E14 = (E2*GG)/(EE*EE);

		B11 = ATG3/BB;
		B12 = -((ATG3*DD)/(BB*BB)) - (B2*DD*F)/(BB*BB*BB) +
		      (B1*DD*(-b + F))/(BB*BB*BB);
		B13 = -((B1*DD)/(BB*BB)) + (B2*DD)/(BB*BB);
		B14 = (B1*DD)/(BB*BB);

		F11 = ATG4/FF;
		F12 = (F2*(-d - F)*HH)/(FF*FF*FF) + (F1*(-b + d + F)*HH)/(FF*FF*FF) -
		      (ATG4*HH)/(FF*FF);
		F13 = -((F1*HH)/(FF*FF)) + (F2*HH)/(FF*FF);
		F14 = (F1*HH)/(FF*FF);

		doubleprecision f1l[3], f2l[3], f3l[3], f4l[3];

		f1l[0] =  - 0.5*invJJ*( (c2[0] - c1[0] ) *( A13 + E13 + B11*b - F11*b + A11*d -
		                                            E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F)
		                                            + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA)
		                                                              + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[0] - c3[0] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[0] - c3[0] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );

		f1l[1] =  - 0.5*invJJ*( (c2[1] - c1[1] ) *( A13 + E13 + B11*b - F11*b + A11*d -
		                                            E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[1] - c3[1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[1] - c3[1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );

		f1l[2] =  - 0.5*invJJ*( (c2[2] - c1[2] ) *( A13 + E13 + B11*b - F11*b + A11*d -
		                                            E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[2] - c3[2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[2] - c3[2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );


		f2l[0] =  - invJJ*( (c2[0] - c1[0] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)
		                                       +(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) )
		                    + 0.5*(c4[0] - c3[0])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e
		                                           - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[0] - c3[0] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );

		f2l[1] = - invJJ*( (c2[1] - c1[1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+
		                                      (E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[1] - c3[1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[1] - c3[1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );

		f2l[2] = - invJJ*( (c2[2] - c1[2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+
		                                      (E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[2] - c3[2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[2] - c3[2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );

		f3l[0] =  - 0.5*invJJ*( (c2[0] - c1[0] )*(-A13 - F13 - B11*b + F11*b - A11*d -
		                                          E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[0] - c3[0] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[0] - c3[0] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );

		f3l[1] =  - 0.5*invJJ*( (c2[1] - c1[1] )*(-A13 - F13 - B11*b + F11*b - A11*d -
		                                          E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[1] - c3[1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[1] - c3[1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) ) ;

		f3l[2] =  - 0.5*invJJ*( (c2[2] - c1[2] )*(-A13 - F13 - B11*b + F11*b - A11*d -
		                                          E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[2] - c3[2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[2] - c3[2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );


		f4l[0] =  - invJJ*( 0.5*(c2[0] - c1[0] )*( -E13 + F13 + 2*E11*d + 2*F11*d -
		                                           4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) +
		                    (c4[0] - c3[0])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e
		                                     - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))
		                    + 0.5*(c1[0] - c3[0] )* (B13 + F13 - A11*a
		                                             + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) )  ;

		f4l[1] =  - invJJ*( 0.5*(c2[1] - c1[1] )*( -E13 + F13 + 2*E11*d + 2*F11*d -
		                                           4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[1] - c3[1])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[1] - c3[1] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) ) ;

		f4l[2] =  - invJJ*( 0.5*(c2[2] - c1[2] )*( -E13 + F13 + 2*E11*d + 2*F11*d -
		                                           4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[2] - c3[2])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[2] - c3[2] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) ) ;


		if(checkNaN_INF(f1l, 0, 2)||checkNaN_INF(f2l,0,2)||checkNaN_INF(f3l, 0, 2)
		   ||checkNaN_INF(f4l,0,2)){
			forceN(coord, f, beadSet, krep, i);
		}
		else{
			for(int dim = 0; dim<3; dim++) {
				f1[dim] += f1l[dim];
				f2[dim] += f2l[dim];
				f3[dim] += f3l[dim];
				f4[dim] += f4l[dim];
			}
		}

#ifdef DETAILEDOUTPUT
		std::cout<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<" "<<
                 f2[0]<<" "<<f2[1]<<" "<<f2[2]<<" "<<
                 f3[0]<<" "<<f3[1]<<" "<<f3[2]<<" "<<
                 f4[0]<<" "<<f4[1]<<" "<<f4[2]<<endl;
#endif
	}
	delete [] newc2;
	delete [] cp;
	delete [] vec_A;
	delete [] vec_B;
	delete [] vec_C;
}

floatingpoint CylinderExclVolRepulsion::energyN(floatingpoint *coord, floatingpoint *force,
                                                int *beadSet, floatingpoint *krep, int i) {
	floatingpoint *c1, *c2, *c3, *c4;

	doubleprecision a, b, c, d, e, F;

	doubleprecision U_i = 0.0;

	int n = CylinderExclVolume<CylinderExclVolRepulsion>::n;

	c1 = &coord[3 * beadSet[n * i]];
	c2 = &coord[3 * beadSet[n * i + 1]];
	c3 = &coord[3 * beadSet[n * i + 2]];
	c4 = &coord[3 * beadSet[n * i + 3]];

	a = scalarProduct(c1, c2, c1, c2);//always positive
	b = scalarProduct(c3, c4, c3, c4);//always positive
	c = scalarProduct(c3, c1, c3, c1);//always positive
	//If c1, c2, c3, and c4 are coplanar, e = d + f
	d = scalarProduct(c1, c2, c3, c4);
	e = scalarProduct(c1, c2, c3, c1);
	F = scalarProduct(c3, c4, c3, c1);

	//Calculate energy

	//Version 2, Simpson's composite rule in 2D.
	int N1  = 400;
	int N2 = 400;
	doubleprecision deltas = 1.0/(doubleprecision(N1));
	doubleprecision deltat = 1.0/(doubleprecision(N2));

	//terms that are independent of s and t.
	doubleprecision Termset1 = getenergyintegrand(a, b, c, d ,e, F, 0.0,0.0) +
	                           getenergyintegrand(a, b, c, d ,e, F, 0.0,1.0) +
	                           getenergyintegrand(a, b, c, d ,e, F, 1.0,0.0) +
	                           getenergyintegrand(a, b, c, d ,e, F, 1.0,1.0);

	//Terms that are summed over s alone
	doubleprecision Termset2 = 0.0;
	//Terms that are summed over t alone
	doubleprecision Termset3 = 0.0;
	//Terms that are summed over both s and t
	doubleprecision Termset4 = 0.0;

	doubleprecision s = 0.0;
	doubleprecision t = 0.0;

	doubleprecision factors = 0.0;
	doubleprecision factort = 0.0;

	for(int i = 0; i <= N1; i++){ //c1-c2
		t = 0.0;
		s += deltas;

		factors = 4.0;
		if(i%2 == 0 ) //even
			factors = 2.0;

		if(i>0 && i < N1){
			Termset2 += factors * (getenergyintegrand(a, b, c, d ,e, F, s, 0.0) + getenergyintegrand(a,
			                                                                                         b, c, d ,e, F, s, 1.0));
		}

		for(int j = 0; j <= N2; j++){ //c3-c4
			t += deltat;

			factort = 4.0;
			if(j%2 == 0 ) //even
				factort = 2.0;

			if(j>0 && j < N2){
				if(i == 0){
					Termset3 += factort * (getenergyintegrand(a, b, c, d ,e, F, 0.0, t) +
					                       getenergyintegrand(a, b, c, d ,e, F, 1.0, t));
				}
				else if (i < N1){
					Termset4 += factors * factort * getenergyintegrand(a, b, c, d ,e, F, s, t);
				}
			}
		}
	}

	U_i = krep[i] * (deltas * deltat/9.0) * (Termset1 + Termset2 + Termset3 + Termset4);

	cout<<"Numerical result_2 "<<U_i<<endl;
	return U_i;
}

void CylinderExclVolRepulsion::forceN(floatingpoint *coord, floatingpoint *f,
                                      int *beadSet, floatingpoint *krep, int i) {
	floatingpoint *c1, *c2, *c3, *c4;
	floatingpoint *f1, *f2, *f3, *f4;
	doubleprecision a, b, c, d, e, F;

	vector<doubleprecision> vecA, vecB, vecC;

	doubleprecision *integrandarray;

	integrandarray = new doubleprecision[6];

	int n = CylinderExclVolume<CylinderExclVolRepulsion>::n;

	c1 = &coord[3 * beadSet[n * i]];
	c2 = &coord[3 * beadSet[n * i + 1]];
	c3 = &coord[3 * beadSet[n * i + 2]];
	c4 = &coord[3 * beadSet[n * i + 3]];

	f1 = &f[3 * beadSet[n * i]];
	f2 = &f[3 * beadSet[n * i + 1]];
	f3 = &f[3 * beadSet[n * i + 2]];
	f4 = &f[3 * beadSet[n * i + 3]];

	a = scalarProduct(c1, c2, c1, c2);//always positive
	b = scalarProduct(c3, c4, c3, c4);//always positive
	c = scalarProduct(c3, c1, c3, c1);//always positive
	//If c1, c2, c3, and c4 are coplanar, e = d + f
	d = scalarProduct(c1, c2, c3, c4);
	e = scalarProduct(c1, c2, c3, c1);
	F = scalarProduct(c3, c4, c3, c1);

	for(int dim = 0; dim < 3; dim++){
		vecA.push_back(c1[dim]-c2[dim]);
		vecB.push_back(c3[dim]-c4[dim]);
		vecC.push_back(c3[dim]-c1[dim]);
	}

	//Version 2, Simpson's composite rule in 2D.
	int N1  = 10000;
	int N2 = 10000;
	doubleprecision deltas = 1.0/(doubleprecision(N1));
	doubleprecision deltat = 1.0/(doubleprecision(N2));

	doubleprecision Termset1[6] = {0.0};
	//Terms that are summed over s alone
	doubleprecision Termset2[6] = {0.0};
	//Terms that are summed over t alone
	doubleprecision Termset3[6] = {0.0};
	//Terms that are summed over both s and t
	doubleprecision Termset4[6] = {0.0};

	doubleprecision limits[2] = {0.0, 1.0};
	doubleprecision s0, t0;

	for(int i = 0; i < 2; i++){
		s0 = limits[i];
		for(int j = 0; j < 2; j++){
			t0 = limits[j];
			getforceintegrand(a, b, c, d, e, F, s0, t0, integrandarray);
			for(int t = 0; t<6; t++){
				Termset1[t] += integrandarray[t];
			}
		}
	}

	doubleprecision s = 0.0;
	doubleprecision t = 0.0;

	doubleprecision factors = 0.0;
	doubleprecision factort = 0.0;

	for(int i = 0; i <= N1; i++){ //c1-c2
		t = 0.0;
		s += deltas;

		factors = 4.0;
		if(i%2 == 0 ) //even
			factors = 2.0;

		if(i>0 && i < N1){
			for(int l = 0; l<2; l++){
				getforceintegrand(a, b, c, d, e, F, s, limits[l], integrandarray);
				for(int t = 0; t<6; t++)
					Termset2[t] += factors * integrandarray[t];
			}
		}

		for(int j = 0; j <= N2; j++){ //c3-c4
			t += deltat;

			factort = 4.0;
			if(j%2 == 0 ) //even
				factort = 2.0;

			if(j>0 && j < N2){

				if(i == 0){
					for(int l = 0; l<2; l++){
						getforceintegrand(a, b, c, d, e, F, limits[l], t, integrandarray);
						for(int t = 0; t<6; t++)
							Termset3[t] += factort * integrandarray[t];
					}
				}
				else if (i < N1){
					getforceintegrand(a, b, c, d, e, F, s, t, integrandarray);
					for(int t = 0; t<6; t++)
						Termset4[t] += factors * factort * integrandarray[t];
				}
			}
		}
	}

	//ID vs integrand
	//0     4/r^6
	//1     4s/r^6
	//2     4s*s/r^6
	//3     4*s*t/r^6
	//4     4*t/r^6
	//5     4*t*t/r^6

	doubleprecision integration[6];
	for(int t = 0; t<6; t++)
		integration[t] = (deltas * deltat/9.0) *(Termset1[t] + Termset2[t] +
		                                         Termset3[t] + Termset4[t]);

	doubleprecision f1l[3], f2l[3], f3l[3], f4l[3];

	for(int dim=0; dim<3; dim++){
		f1l[dim] =4.0*krep[i] * ( -integration[0]*vecC[dim] + integration[1]*(vecC[dim]- vecA[dim])
		                          + integration[2]*vecA[dim] + integration[4] * vecB[dim] - integration[3]
		                                                                                    *vecB[dim]);
		f2l[dim] = 4.0*krep[i] * ( -integration[1]*vecC[dim] - integration[2]*vecA[dim] +
		                           integration[3]*vecB[dim]);
		f3l[dim] =4.0*krep[i] * (  integration[0]*vecC[dim] + integration[1]*vecA[dim] -
		                           integration[4]* (vecB[dim]+vecC[dim]) - integration[3]*vecA[dim] +
		                           integration[5]*vecB[dim]);
		f4l[dim] = 4.0*krep[i] * ( integration[4]*vecC[dim] + integration[3]*vecA[dim]
		                           -integration[5]*vecB[dim]);
	}

//	cout<<"printing numerical forces"<<endl;
/*	cout<<f1l[0]<<" "<<f1l[1]<<" "<<f1l[2]<<endl;
	cout<<f2l[0]<<" "<<f2l[1]<<" "<<f2l[2]<<endl;
	cout<<f3l[0]<<" "<<f3l[1]<<" "<<f3l[2]<<endl;
	cout<<f4l[0]<<" "<<f4l[1]<<" "<<f4l[2]<<endl;*/

	delete [] integrandarray;

	if(checkNaN_INF(f1l, 0, 2)||checkNaN_INF(f2l,0,2)||checkNaN_INF(f3l, 0, 2)
	   ||checkNaN_INF(f4l,0,2)){

		cout<<"Cylinder Exclusion Force becomes infinite. Printing data "<<endl;

		short found = 0;
		Cylinder *cyl1, *cyl2;
		for(auto cyl:Cylinder::getCylinders()){
			auto dbIndex1 = cyl->getFirstBead()->_dbIndex;
			auto dbIndex2 = cyl->getSecondBead()->_dbIndex;
			if(dbIndex1 == beadSet[n * i] && dbIndex2 == beadSet[n * i + 1]) {
				cyl1 = cyl;
				found++;
				if(found>=2)
					break;
			}
			else if(dbIndex1 == beadSet[n * i + 2] && dbIndex2 == beadSet[n * i + 3]){
				cyl2 = cyl;
				found++;
				if(found>=2)
					break;
			}
		}
		cout<<"Cylinder IDs "<<cyl1->getId()<<" "<<cyl2->getId()<<" with cIndex "
		    <<cyl1->_dcIndex<<" "<<cyl2->_dcIndex<<" and bIndex "
		    <<cyl1->getFirstBead()->_dbIndex<<" "
		    <<cyl1->getSecondBead()->_dbIndex<<" "
		    <<cyl2->getFirstBead()->_dbIndex<<" "
		    <<cyl2->getSecondBead()->_dbIndex<<endl;

		cout<<"Printing coords"<<endl;
		cout<<c1[0]<<" "<<c1[1]<<" "<<c1[2]<<endl;
		cout<<c2[0]<<" "<<c2[1]<<" "<<c2[2]<<endl;
		cout<<c3[0]<<" "<<c3[1]<<" "<<c3[2]<<endl;
		cout<<c4[0]<<" "<<c4[1]<<" "<<c4[2]<<endl;
		cout<<"Printing force"<<endl;
		cout<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<endl;
		cout<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<endl;
		cout<<f3[0]<<" "<<f3[1]<<" "<<f3[2]<<endl;
		cout<<f4[0]<<" "<<f4[1]<<" "<<f4[2]<<endl;
		cout<<"Printing binary Coords"<<endl;
		printvariablebinary(c1,0,2);
		printvariablebinary(c2,0,2);
		printvariablebinary(c3,0,2);
		printvariablebinary(c4,0,2);
		cout<<"Printing binary Force"<<endl;
		printvariablebinary(f1,0,2);
		printvariablebinary(f2,0,2);
		printvariablebinary(f3,0,2);
		printvariablebinary(f4,0,2);

		cout<<"Printing infinite energy contributions "<<endl;
		cout<<"a "<<a<<" b "<<b<<" c "<<c<<" d "<<d<<" e "<<e<<" f "<<F<<endl;

		exit(EXIT_FAILURE);
	}
	else{
		for(int dim = 0; dim<3; dim++) {
			f1[dim] += f1l[dim];
			f2[dim] += f2l[dim];
			f3[dim] += f3l[dim];
			f4[dim] += f4l[dim];
		}
	}
}