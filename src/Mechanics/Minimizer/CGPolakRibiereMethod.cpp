

//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "CGPolakRibiereMethod.h"

#include "ForceFieldManager.h"
#include "Composite.h"
#include "Output.h"
#include "cross_check.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "Structure/Bead.h"

MinimizationResult PolakRibiere::minimize(ForceFieldManager &FFM, floatingpoint GRADTOL,
                            floatingpoint MAXDIST, floatingpoint LAMBDAMAX,
                            floatingpoint LAMBDARUNNINGAVERAGEPROBABILITY,
                            string _LINESEARCHALGORITHM,
                            bool steplimit) {

	#ifdef TRACKDIDNOTMINIMIZE
	SysParams::Mininimization().beta.clear();
	SysParams::Mininimization().Lambda.clear();
	SysParams::Mininimization().Energyvec.clear();
	SysParams::Mininimization().TotalE.clear();
	SysParams::Mininimization().maxF.clear();
	SysParams::Mininimization().safeModeORnot.clear();
	SysParams::Mininimization().tempEnergyvec.clear();
	SysParams::Mininimization().gradientvec.clear();
	#endif

	MinimizationResult result;
#ifdef CUDATIMETRACK
	chrono::high_resolution_clock::time_point tbeginTot, tendTot;
	chrono::high_resolution_clock::time_point tbeginII, tendII;
	tbeginTot = chrono::high_resolution_clock::now();
	tbeginII = chrono::high_resolution_clock::now();
#endif
	//number of steps
	int N;
	if (steplimit) {
		int beadMaxStep = 3 * Bead::numBeads();
		N = (beadMaxStep > _MINNUMSTEPS ? beadMaxStep : _MINNUMSTEPS);
	} else
		N = numeric_limits<int>::max();

	//@@@{ STEP 1: Start minimization
	tbegin = chrono::high_resolution_clock::now();
	startMinimization();//TODO needs to be hostallocdefault and MemCpyAsync followed by CudaStreamSynchronize

#ifdef ALLSYNC
	cudaDeviceSynchronize();
#endif
	FFM.vectorizeAllForceFields();//each forcefield needs to use hostallocdefault and MemCpyAsync followed by CudaStreamSynchronize

#ifdef ALLSYNC
	cudaDeviceSynchronize();
#endif
	tend = chrono::high_resolution_clock::now();
	chrono::duration<floatingpoint> elapsed_vec(tend - tbegin);
	CUDAcommon::tmin.vectorize += elapsed_vec.count();
	//@@@} Start Minimization
#ifdef OPTIMOUT
	floatingpoint lambdatime = 0.0;
	int safestatuscount = 0;
#endif

#ifdef CUDAACCL
	cross_checkclass::Aux=false;
	auto cvars = CUDAcommon::getCUDAvars();
	cvars.streamvec.clear();
	CUDAcommon::cudavars = cvars;
#endif
#ifdef ALLSYNC
	cudaDeviceSynchronize();
#endif
	//@@@{ STEP 2: COMPUTE FORCES
	tbegin = chrono::high_resolution_clock::now();
	FFM.computeForces(Bead::getDbData().coords.data(),
	                  Bead::getDbData().forces.data()); //split and synchronize in the end
	tend = chrono::high_resolution_clock::now();
	chrono::duration<floatingpoint> elapsed_force(tend - tbegin);
	CUDAcommon::tmin.computeforces += elapsed_force.count();
	//@@@}

	//Output energy
//	cout<<"Energy before minimization"<<endl;
//	FFM.computeEnergy(coord, force, 0.0, true);
#ifdef SERIAL // SERIAL
	//@@@{ STEP 3: COPY FORCES
	tbegin = chrono::high_resolution_clock::now();
	Bead::getDbData().forcesAux = Bead::getDbData().forces;
	Bead::getDbData().forcesAuxP = Bead::getDbData().forces;
	auto maxForce = maxF();

	result.energiesBefore = FFM.computeEnergyHRMD(Bead::getDbData().coords.data());

	tend = chrono::high_resolution_clock::now();
	chrono::duration<floatingpoint> elapsed_copy(tend - tbegin);
	CUDAcommon::tmin.copyforces += elapsed_copy.count();
	//@@@}
#endif
	//M as the first letter in variables signifies that it is used by minimizer
	// (as opposed to finding lambda)
	bool Ms_isminimizationstate, Ms_issafestate;
	int numIter = 0;
	floatingpoint lambda = 0.0;
#ifdef CUDAACCL
	volatile bool *Mc_isminimizationstate;
	volatile bool *Mc_issafestate;
	Ms_isminimizationstate = false;
	Ms_issafestate = false;
#endif
#ifdef SERIAL
	//TODO Comment during SERIAL_CUDACROSSCHECK @{
	bool *Mc_isminimizationstate;
	bool *Mc_issafestate;
	bool *M_ETolstate;
//    @}

	Mc_isminimizationstate = new bool[1];
	Mc_issafestate = new bool[1];
	M_ETolstate = new bool[1];
	Mc_isminimizationstate[0] = false;//points to address of Mmh_stop
	Mc_issafestate[0] = false;//points to address of Msh_stop
	M_ETolstate[0] = false;
#endif
#ifdef CUDATIMETRACK
	chrono::high_resolution_clock::time_point tbegin, tend;
#endif
#ifdef  CUDAACCL
#ifdef CUDATIMETRACK
	//    chrono::high_resolution_clock::time_point tbegin, tend;
		tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		for(auto strm:CUDAcommon::getCUDAvars().streamvec)
			CUDAcommon::handleerror(cudaStreamSynchronize(*strm));

#ifdef CUDATIMETRACK
		tend= chrono::high_resolution_clock::now();
		chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
		CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run.count());
		CUDAcommon::cudatime.TcomputeF += elapsed_run.count();
#endif
#ifdef CUDATIMETRACK

		//Reset lambda time tracker.
		CUDAcommon::cudatime.Tlambda = 0.0;
#endif

		if(!(CUDAcommon::getCUDAvars().conservestreams) || stream1 == NULL)
			CUDAcommon::handleerror(cudaStreamCreate(&stream1));
		if(!(CUDAcommon::getCUDAvars().conservestreams) || stream1 == NULL)
			CUDAcommon::handleerror(cudaStreamCreate(&stream2));
		if(!(CUDAcommon::getCUDAvars().conservestreams) || stream1 == NULL)
			CUDAcommon::handleerror(cudaStreamCreate(&stream3));
		FFM.CUDAcopyForces(stream1, CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_force);//pass a
		// stream
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		FFM.CUDAcopyForces(stream2, CUDAcommon::getCUDAvars().gpu_forceAuxP,CUDAcommon::getCUDAvars().gpu_force);//pass a
		// stream
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		floatingpoint *gpu_GRADTOL;
		floatingpoint gradtol[1];
		gradtol[0]= GRADTOL;
		CUDAcommon::handleerror(cudaMalloc((void **) &gpu_GRADTOL, sizeof(floatingpoint)));
		CUDAcommon::handleerror(cudaMemcpy(gpu_GRADTOL, gradtol, sizeof(floatingpoint), cudaMemcpyHostToDevice));
		CGMethod::CUDAallFDotF(stream3);//curGrad //pass a stream

		//synchronize streams
		CUDAcommon::handleerror(cudaStreamSynchronize(stream1));
		CUDAcommon::handleerror(cudaStreamSynchronize(stream2));
		CUDAcommon::handleerror(cudaStreamSynchronize(stream3));
		if(!(CUDAcommon::getCUDAvars().conservestreams)) {
			CUDAcommon::handleerror(cudaStreamDestroy(stream1));
			CUDAcommon::handleerror(cudaStreamDestroy(stream2));
			CUDAcommon::handleerror(cudaStreamDestroy(stream3));
		}
	//PING PONG
		bool  *Mmh_stop, *Mmg_stop1, *Mmg_stop2, *Mmg_s1, *Mmg_s2, *Mmg_ss;//minimization state
		bool  *Msh_stop, *Msg_stop1, *Msg_stop2, *Msg_s1, *Msg_s2, *Msg_ss;//safe state

		//PING PONG
		//minimization state
		cudaMalloc(&Mmg_stop1, sizeof(bool));
		cudaMalloc(&Mmg_stop2, sizeof(bool));
		cudaHostAlloc(&Mmh_stop, sizeof(bool), cudaHostAllocDefault);
		//safe state
		cudaMalloc(&Msg_stop1, sizeof(bool));
		cudaMalloc(&Msg_stop2, sizeof(bool));
		cudaHostAlloc(&Msh_stop, sizeof(bool), cudaHostAllocDefault);
		//@

		if(!(CUDAcommon::getCUDAvars().conservestreams) || Ms1 == NULL)
			CUDAcommon::handleerror(cudaStreamCreate(&Ms1));
		if(!(CUDAcommon::getCUDAvars().conservestreams) || Ms2 == NULL)
			CUDAcommon::handleerror(cudaStreamCreate(&Ms2));
		if(!(CUDAcommon::getCUDAvars().conservestreams) || Ms3 == NULL)
			CUDAcommon::handleerror(cudaStreamCreate(&Ms3));
		if(!(CUDAcommon::getCUDAvars().conservestreams) || Ms4 == NULL)
			CUDAcommon::handleerror(cudaStreamCreate(&Ms4));
		if(!(CUDAcommon::getCUDAvars().conservestreams) || Me1 == NULL)
			CUDAcommon::handleerror(cudaEventCreate(&Me1));
		if(!(CUDAcommon::getCUDAvars().conservestreams) || Me2 == NULL)
			CUDAcommon::handleerror(cudaEventCreate(&Me2));
		if(!(CUDAcommon::getCUDAvars().conservestreams) || stream_shiftsafe == NULL)
			CUDAcommon::handleerror(cudaStreamCreate(&stream_shiftsafe));
		if(!(CUDAcommon::getCUDAvars().conservestreams) || stream_dotcopy == NULL)
			CUDAcommon::handleerror(cudaStreamCreate(&stream_dotcopy));
		if(!(CUDAcommon::getCUDAvars().conservestreams) || event_safe == NULL)
			CUDAcommon::handleerror(cudaEventCreate(&event_safe));
		if(!(CUDAcommon::getCUDAvars().conservestreams) || event_dot == NULL)
			CUDAcommon::handleerror(cudaEventCreate(&event_dot));

		Mmh_stop[0] = true; //Minimizationstate //Yes = Minimize. No = Don't minimize.
		Msh_stop[0] = false; //safe state
		Mc_isminimizationstate = Mmh_stop;//points to address of Mmh_stop
		Mc_issafestate = Msh_stop;//points to address of Msh_stop

		Msp1 = &Ms1;
		Msp2 = &Ms2;
		Mep1 = &Me1;
		Mep2 = &Me2;
		Mmg_s1 = Mmg_stop1;
		Mmg_s2 = Mmg_stop2;
		Msg_s1 = Msg_stop1;
		Msg_s2 = Msg_stop2;
	//set Mmg_stop1, Mmg_stop2 to true and Msg_stop1, Msg_stop2 to false.

	// stick to single stream going forward.
	//@CUDA Get minimizaton state{
		//calculate MAXF and allFDotFA
		CGMethod::CUDAinitializePolak(*Msp1, Mmg_s1, Mmg_s2, Msg_s1, Msg_s2);
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
		CUDAcommon::handleerror(cudaGetLastError(),"CUDAgetPolakvars", "CGPolakRibiereMethod.cu");
#ifdef SERIAL_CUDACROSSCHECK
		CUDAcommon::handleerror(cudaDeviceSynchronize());
		std::cout<<"FMAX SERL "<<maxF()<<endl;
#endif
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		//Copy to host
		CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
		cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
	//@}
#endif

#ifdef SERIAL //SERIAL
	//@@@{ STEP 4 OTHER
	//FIND MAXIMUM ERROR BETWEEN CUDA AND VECTORIZED FORCES{
	//VECTORIZED. Prep for Polak{
	tbegin = chrono::high_resolution_clock::now();
	floatingpoint curGrad = CGMethod::allFDotF();
	Ms_isminimizationstate = true;
	Ms_issafestate = false;
	Ms_isminimizationstate = maxForce > GRADTOL;
	bool ETOLexittstatus = false;
	//
#ifdef DETAILEDOUTPUT
	std::cout<<"printing beads & forces"<<endl;
	long i = 0;
	long index = 0;
	for(auto b:Bead::getBeads()){
		index = 3 * b->getStableIndex();
		std::cout<<b->getId()<<" "<< b->coordinate() <<" "
				"" << b->force() <<endl;
	}
	std::cout<<"printed beads & forces"<<endl;
#endif
	//
#endif
#ifdef CUDATIMETRACK_MACRO
	CUDAcommon::cudatime.Tlambda = 0.0;
	CUDAcommon::cudatime.Ecount = 0;
	CUDAcommon::cudatime.TcomputeE = 0.0;
	CUDAcommon::cudatime.TcomputeF = 0.0;
	floatingpoint c1 = 0.0;
	floatingpoint c2,c3,c4,c5;
	c2 = 0.0;c3 = 0.0;c4 = 0.0;c5=0.0;
#endif
#ifdef CUDATIMETRACK
	tendII= chrono::high_resolution_clock::now();
	chrono::duration<floatingpoint> elapsed_runslice1(tendII - tbeginII);
	std::cout<<"Slice time "<<elapsed_runslice1.count()<<endl;
	tbeginII = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
	while (/* Iteration criterion */  numIter < N &&
		   /* Gradient tolerance  */  (Mc_isminimizationstate[0])) {
//#ifdef CUDATIMETRACK_MACRO
//        chrono::high_resolution_clock::time_point tbeginiter, tenditer;
//        tbeginiter = chrono::high_resolution_clock::now();
//#endif
#ifdef CUDATIMETRACK
		chrono::high_resolution_clock::time_point tbegin, tend;
#endif

#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
//PING PONG SWAP
//        CUDAcommon::handleerror(cudaStreamWaitEvent(*Msp2, *Mep1, 0));

#ifdef CUDATIMETRACK
		tbegin = chrono::high_resolution_clock::now();
#endif
		CUDAcommon::handleerror(cudaStreamSynchronize(*Msp1));
		CUDAcommon::handleerror(cudaStreamSynchronize(stream_shiftsafe));
#ifdef CUDATIMETRACK
		tend = chrono::high_resolution_clock::now();
		chrono::duration<floatingpoint> elapsed_run1(tend - tbegin);
		c5 += elapsed_run1.count();
#endif

		Msps = Msp1;
		Msp1 = Msp2;
		Msp2 = Msps;
		Meps = Mep1;
		Mep1 = Mep2;
		Mep2 = Meps;
		Mmg_ss = Mmg_s1;
		Mmg_s1 = Mmg_s2;
		Mmg_s2 = Mmg_ss;
		Msg_ss = Msg_s1;
		Msg_s1 = Msg_s2;
		Msg_s2 = Msg_ss;
//PING ENDS
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		numIter++;
#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_LAMBDA)
		std::cout<<"SL safestate "<<_safeMode<<endl;
#endif
		if(Mc_issafestate[0]) {
			_safeMode = false;
		}
#ifdef CUDATIMETRACK_MACRO
		CUDAcommon::serltime.TcomputeEiter = 0.0;
		CUDAcommon::cudatime.TcomputeEiter = 0.0;
		chrono::high_resolution_clock::time_point tLbegin, tLend;
		tLbegin = chrono::high_resolution_clock::now();
#endif
//        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp2));//make sure previous iteration is done.
		//find lambda by line search, move beads
		lambda = backtrackingLineSearchCUDA(FFM, MAXDIST, LAMBDAMAX, Msg_s1);
#ifdef CUDATIMETRACK_MACRO
		tLend= chrono::high_resolution_clock::now();
		chrono::duration<floatingpoint> elapsed_runiter(tLend - tLbegin);
		std::cout<<" CUDA Iter "<<numIter<<" Lambda Time taken (s) "<<elapsed_runiter
				.count() - CUDAcommon::serltime.TcomputeEiter <<endl;
#endif
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif

#ifdef CUDATIMETRACK
		tbegin = chrono::high_resolution_clock::now();
#endif
		CUDAmoveBeads(*Msp1, Mmg_s1);
//        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));//This seems unnecessary.
		//wait for movebeads to finish before calculating forces1

#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		//Synchronize to ensure forces are calculated on moved beads.
		CUDAcommon::handleerror(cudaStreamSynchronize(*Msp1));
//        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp2));
		//@CUDA Get minimizaton state
		// @{
//        CGMethod::CUDAgetPolakvars(true, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);
		CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
		//@}
#ifdef CUDATIMETRACK
		tend = chrono::high_resolution_clock::now();
		chrono::duration<floatingpoint> elapsed_run1b(tend - tbegin);
		c1 += elapsed_run1b.count();
#endif

#if defined(CROSSCHECK) || defined(CUDAACCL)
		cross_checkclass::Aux=true;
#endif
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif

		cvars = CUDAcommon::getCUDAvars();
		cvars.streamvec.clear();
		CUDAcommon::cudavars = cvars;
#ifdef CUDATIMETRACK
		tbegin = chrono::high_resolution_clock::now();
#endif
		CUDAcommon::handleerror(cudaStreamSynchronize(stream_dotcopy));
#ifdef CUDATIMETRACK
		tend = chrono::high_resolution_clock::now();
		chrono::duration<floatingpoint> elapsed_run2(tend - tbegin);
		c2 += elapsed_run2.count();
#endif

		//compute new forces
		FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forcesAux.data());//split and synchronize
#ifdef DETAILEDOUTPUT
		std::cout<<"MB printing beads & forces L "<<lambda<<endl;
		long i = 0;
		long index = 0;
		for(auto b:Bead::getBeads()){
			index = 3 * b->getStableIndex();

			std::cout<<b->getId()<<" "<<coord[index]<<" "<<coord[index + 1]<<" "
					""<<coord[index + 2]<<" "
					""<<forceAux[index]<<" "
					""<<forceAux[index + 1]<<" "<<forceAux[index + 2]<<" "<<3 *
					b->getStableIndex()<<endl;
		}
		std::cout<<"MB printed beads & forces"<<endl;
#endif

#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
#ifdef CUDATIMETRACK
		tbegin = chrono::high_resolution_clock::now();
#endif
		//wait for forces to be calculated
		for(auto strm:CUDAcommon::getCUDAvars().streamvec)
			CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
#ifdef CUDATIMETRACK
		tend= chrono::high_resolution_clock::now();
		chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
		CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run.count());
		CUDAcommon::cudatime.TcomputeF += elapsed_run.count();
#endif
#ifdef CUDATIMETRACK
		tbegin = chrono::high_resolution_clock::now();
#endif
//Copying forces back to Host was here.
#ifdef CUDATIMETRACK
		tend = chrono::high_resolution_clock::now();
		chrono::duration<floatingpoint> elapsed_run1c(tend - tbegin);
		c1 += elapsed_run1c.count();
#endif

#ifdef CUDATIMETRACK
		tbegin = chrono::high_resolution_clock::now();
#endif
		//compute direction CUDA
		CGMethod::CUDAallFADotFA(stream_dotcopy); //newGrad
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		CGMethod::CUDAallFADotFAP(stream_dotcopy); //prevGrad
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		CUDAcommon::handleerror(cudaEventRecord(event_dot,stream_dotcopy));
		CUDAcommon::handleerror(cudaStreamWaitEvent(stream_shiftsafe, event_dot,0));
//Copy forces
		FFM.CUDAcopyForces(stream_dotcopy, CUDAcommon::getCUDAvars().gpu_forceAuxP,CUDAcommon::getCUDAvars().gpu_forceAux);
#ifdef CUDATIMETRACK
		tend = chrono::high_resolution_clock::now();
		chrono::duration<floatingpoint> elapsed_run2c(tend - tbegin);
		c2 += elapsed_run2c.count();
#endif
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
#ifdef CUDATIMETRACK
		tbegin = chrono::high_resolution_clock::now();
#endif
		//Polak-Ribieri update beta & shift gradient
		CUDAshiftGradient(stream_shiftsafe, Mmg_s1);
#ifdef SERIAL_CUDACROSSCHECK
		CUDAcommon::handleerror(cudaStreamSynchronize(stream_shiftsafe));
#endif
/*        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp2));
		//@CUDA Get minimizaton state{
//        CGMethod::CUDAgetPolakvars(true, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);
		CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));*/
		CGMethod::CUDAgetPolakvars2(stream_shiftsafe, Msg_s2);
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
		CUDAshiftGradientifSafe(stream_shiftsafe, Mmg_s1, Msg_s1);
#ifdef ALLSYNC
		cudaDeviceSynchronize();
#endif
#ifdef CUDATIMETRACK
		tend = chrono::high_resolution_clock::now();
		chrono::duration<floatingpoint> elapsed_run3(tend - tbegin);
		c3 += elapsed_run3.count();
#endif
#ifdef CUDATIMETRACK
		tbegin = chrono::high_resolution_clock::now();
#endif
		if(Mc_isminimizationstate[0]  == true){
			//Copy to host
			CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
//            CUDAcommon::handleerror(cudaStreamWaitEvent(Ms4, event_safe, 0));//event_safe is not attached to any event
#ifdef ALLSYNC
			cudaDeviceSynchronize();
#endif
			cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
#ifdef ALLSYNC
			cudaDeviceSynchronize();
#endif
			//TODO remove later... June 7, 2018. Removed.
//            std::cout<<"safe state copy"<<endl;
//            cudaMemcpyAsync(Msh_stop, Msg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms4);

		}
#ifdef CUDATIMETRACK
		tend = chrono::high_resolution_clock::now();
		chrono::duration<floatingpoint> elapsed_run4(tend - tbegin);
		c4 += elapsed_run4.count();
#endif
	}

#ifdef CUDATIMETRACK_MACRO
	std::cout<<"CUDA start min time taken (s) "<<CUDAcommon::cudatime.Tstartmin<<endl;
	std::cout<<"CUDA Energy time taken (s) "<<CUDAcommon::cudatime.TcomputeE<<"for total "
			"iters "<<CUDAcommon::cudatime.Ecount<<endl;
	std::cout<<"CUDA Force time taken (s) "<<CUDAcommon::cudatime.TcomputeF<<endl;
	std::cout<<"CUDA Energy time per iter (s/iter) "<<CUDAcommon::cudatime.TcomputeE/
			(floatingpoint(CUDAcommon::cudatime.Ecount))<<endl;
	std::cout<<"CUDA Split Times Iter "<<numIter<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "
			""<<c5<<endl;
	std::cout<<"CUDA Add "<<c5+c4+c3+c2+c1<<endl;
#endif
	std::cout<<"CUDA Total number of iterations "<<numIter<<endl;
#endif //CUDAACCL

#ifdef CUDATIMETRACK_MACRO
	floatingpoint s1 = 0.0;
	floatingpoint s2,s3,s4;
	s2 = 0.0;s3 = 0.0;s4 = 0.0;
	CUDAcommon::serltime.Tlambda = 0.0;
	CUDAcommon::serltime.Ecount = 0;
	CUDAcommon::serltime.TcomputeE = 0.0;
	CUDAcommon::serltime.TcomputeF = 0.0;
#endif
	_safeMode = false;
#ifdef OPTIMOUT
	std::cout << "----------------------------------------" << endl;
	std::cout << "maxF " << maxForce << endl;
#endif

	#ifdef TRACKDIDNOTMINIMIZE
	SysParams::Mininimization().maxF.push_back(maxForce);
	#endif

	tend = chrono::high_resolution_clock::now();
	chrono::duration<floatingpoint> elapsed_other(tend - tbegin);
	CUDAcommon::tmin.tother += elapsed_other.count();

/*	for(auto cyl:Cylinder::getCylinders()){
		cout<<"before while minimization  Cylinder ID = "<<cyl->getId()<<endl;
		cyl->printSelf();
	}*/
			//@@@} STEP 4 OTHER
#ifdef SERIAL
			while (/* Iteration criterion */  numIter < N &&
			/* Gradient tolerance  */  (Ms_isminimizationstate) &&
			                                  !ETOLexittstatus) {

//#ifdef CUDATIMETRACK_MACRO
//        chrono::high_resolution_clock::time_point tbeginiter, tenditer;
//        tbeginiter = chrono::high_resolution_clock::now();
//#endif

				//@@@{ STEP 5 OTHER
				tbegin = chrono::high_resolution_clock::now();
				if (std::is_same<floatingpoint, float>::value) {
					//set the floor of lambda (lowest lambda allowed based on maxf
					int maxForder = static_cast<int>(floor(log10(maxForce)));
					if (maxForder < 0) maxForder--;
					CGMethod::setLAMBDATOL(maxForder);
				}

				tend = chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_other2(tend - tbegin);
				CUDAcommon::tmin.tother += elapsed_other2.count();
				//@@@} OTHER

				double beta, newGrad, prevGrad;
//        std::cout<<"SERL maxF "<<maxF()<<endl;

				numIter++;
#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_LAMBDA)
				std::cout<<"SL safestate "<<_safeMode<<endl;
#endif
#ifdef CUDATIMETRACK_MACRO
				CUDAcommon::serltime.TcomputeEiter = 0.0;
				CUDAcommon::cudatime.TcomputeEiter = 0.0;
				chrono::high_resolution_clock::time_point tLbegin, tLend;
				tLbegin = chrono::high_resolution_clock::now();
#endif
#ifdef TRACKDIDNOTMINIMIZE
				SysParams::Mininimization().safeModeORnot.push_back(_safeMode);
#endif

				//@@@{ STEP 6 FIND LAMBDA
				prevlambda = lambda;
				tbegin = chrono::high_resolution_clock::now();
				bool copysafeMode = _safeMode;
				bool *dummy = nullptr;
				if (_LINESEARCHALGORITHM == "BACKTRACKING") {
					lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, maxForce,
					                                                LAMBDAMAX, dummy,
					                                                M_ETolstate)
					                   : backtrackingLineSearch(FFM, MAXDIST, maxForce,
					                                            LAMBDAMAX,
					                                            LAMBDARUNNINGAVERAGEPROBABILITY,
					                                            dummy, M_ETolstate);
				} else if (_LINESEARCHALGORITHM == "QUADRATIC") {
					lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, maxForce,
					                                                LAMBDAMAX, dummy,
					                                                M_ETolstate)
					                   : quadraticLineSearchV2(FFM, MAXDIST, maxForce,
					                                           LAMBDAMAX,
					                                           LAMBDARUNNINGAVERAGEPROBABILITY,
					                                           dummy, M_ETolstate);
				} else {
					lambda = _safeMode ? safeBacktrackingLineSearchV2(FFM, MAXDIST,
					                                                  maxForce,
					                                                  LAMBDAMAX, dummy,
					                                                  M_ETolstate)
					                   : quadraticLineSearchV2(FFM, MAXDIST, maxForce,
					                                           LAMBDAMAX,
					                                           LAMBDARUNNINGAVERAGEPROBABILITY,
					                                           dummy, M_ETolstate);
				}

				tend = chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_lambda(tend - tbegin);
				CUDAcommon::tmin.findlambda += elapsed_lambda.count();
				#ifdef OPTIMOUT
				lambdatime += elapsed_lambda.count();
				#endif
				////@@@@} FIND LAMBDA

#ifdef TRACKDIDNOTMINIMIZE
				SysParams::Mininimization().Lambda.push_back(lambda);
#endif

#ifdef CUDATIMETRACK_MACRO
				tLend= chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_runiter(tLend - tLbegin);
				std::cout<<"SERL Iter "<<numIter<<" Lambda Time taken (s) "<<elapsed_runiter.count() - CUDAcommon::cudatime.TcomputeEiter<<endl;
#endif

#ifdef SERIAL_CUDACROSSCHECK
				CUDAcommon::handleerror(cudaDeviceSynchronize());
				floatingpoint cuda_lambda[1];
				CUDAcommon::handleerror(cudaMemcpy(cuda_lambda, CUDAcommon::cudavars.gpu_lambda,  sizeof(floatingpoint),
												   cudaMemcpyDeviceToHost));
				std::cout<<"Lambda CUDA "<<cuda_lambda[0]<<" SERL "<<lambda<<endl;
#endif
#ifdef CUDATIMETRACK
				tbegin = chrono::high_resolution_clock::now();
#endif
				//@@@{ STEP7 OTHER
				if (Ms_isminimizationstate) {
#if defined(TRACKDIDNOTMINIMIZE) || defined(EVSALPHA)
					//Backup coordinate
					/*const std::size_t num = Bead::getDbData().coords.size_raw();
					Bead::getDbData().coords_bckup.resize(num);
					Bead::getDbData().forces_bckup.resize(num);

					for (size_t i = 0; i < num; ++i) {
						Bead::getDbData().coords_bckup.value[i] = Bead::getDbData().coords.value[i];
						Bead::getDbData().forces_bckup.value[i] = Bead::getDbData().forces.value[i];
					}
					calculateEvsalpha(FFM, lambda, LAMBDAMAX, allFDotFA());
					cout<<endl;*/

#endif
					tbegin = chrono::high_resolution_clock::now();
					//SERIAL VERSION
					moveBeads(lambda);
					tend = chrono::high_resolution_clock::now();
					chrono::duration<floatingpoint> elapsed_other3(tend - tbegin);
					CUDAcommon::tmin.tother += elapsed_other3.count();
				}
				//@@@} OTHER
#if defined(CROSSCHECK) || defined(CUDAACCL)
				cross_checkclass::Aux=true;
#endif
#ifdef CUDATIMETRACK
				tend = chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_runs1(tend - tbegin);
				s1 += elapsed_runs1.count();
#endif
				///@@@{ STEP 8 compute new forces
				tbegin = chrono::high_resolution_clock::now();
				FFM.computeForces(Bead::getDbData().coords.data(),
				                  Bead::getDbData().forcesAux.data());//split and synchronize
				tend = chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_force(tend - tbegin);
				CUDAcommon::tmin.computeforces += elapsed_force.count();

				maxForce = maxF();

				if (M_ETolstate[0] && maxForce <= 2.5 * GRADTOL) {
					ETOLexittstatus = true;
				} else
					M_ETolstate[0] = false;
				///@@@}
#ifdef DETAILEDOUTPUT
				std::cout<<"MB printing beads & forces L "<<lambda<<endl;
				long i = 0;
				long index = 0;
				for(auto b:Bead::getBeads()){
					index = 3 * b->getStableIndex();

					std::cout<<b->getId()<<" "<<coord[index]<<" "<<coord[index + 1]<<" "
							""<<coord[index + 2]<<" "
							""<<forceAux[index]<<" "
							""<<forceAux[index + 1]<<" "<<forceAux[index + 2]<<" "<<3 *
							b->getStableIndex()<<endl;
				}
				std::cout<<"MB printed beads & forces"<<endl;
#endif
#ifdef CUDATIMETRACK
				tbegin = chrono::high_resolution_clock::now();
#endif
				//@@@{ STEP 9 OTHER
				tbegin = chrono::high_resolution_clock::now();
				//compute direction
//        std::cout<<"serial"<<endl;
				newGrad = CGMethod::allFADotFA();
				prevGrad = CGMethod::allFADotFAP();
#ifdef CUDATIMETRACK
				tend = chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_runs2a(tend - tbegin);
				s2 += elapsed_runs2a.count();
#endif
#ifdef CUDATIMETRACK
				tbegin = chrono::high_resolution_clock::now();
#endif

				//Polak-Ribiere update
				//Max(0,betaPR) allows us to reset the direction under non-ideal circumstances.
				//The direction is reset of steepest descent direction (-gk).
				double betaPR = max<double>((double) 0.0, (newGrad - prevGrad) / curGrad);
				double betaFR = max<double>((double) 0.0, newGrad / curGrad);
				//Efficient hybrid Conjugate gradient techniques, Eq 21
				prevbeta = beta;
				if (betaPR == 0.0)
					beta = betaFR;
				else if (betaPR < 1.25 * betaFR)
					beta = betaPR;
				else
					beta = betaFR;

				//Global convergence properties of conjugate gradient methods for optimization Eq
				// 3.8
				//A SURVEY OF NONLINEAR CONJUGATE GRADIENT METHODS Section 9.
				//Allows for negative beta values.
				/*double betaPR = (newGrad - prevGrad) / curGrad;
				double betaFR = newGrad/ curGrad;
				beta = max<double>(-betaFR, min<double>(betaPR, betaFR));
				cout<<"betaPR "<<betaPR<<" betaFR "<<betaFR<<" beta "<<beta<<endl;*/

//	    cout<<"newGrad "<<newGrad<<" prevGrad "<<prevGrad<<" curGrad "
//	    <<curGrad<<" beta "<<beta<<endl;

//        cout<<"lambda "<<lambda<<" beta "<<beta<<endl;
				if (Ms_isminimizationstate)
					//shift gradient
					shiftGradient(beta);

				tend = chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_other4(tend - tbegin);
				CUDAcommon::tmin.tother += elapsed_other4.count();
				//@@@} OTHER

#ifdef TRACKDIDNOTMINIMIZE
				SysParams::Mininimization().beta.push_back(beta);
				SysParams::Mininimization().maxF.push_back(maxForce);
#endif

#ifdef CUDATIMETRACK
				tend = chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_runs3a(tend - tbegin);
				s3 += elapsed_runs3a.count();
#endif
#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_BETA)
				std::cout<<"Shift Gradient "<<beta<<endl;
				CUDAcommon::handleerror(cudaDeviceSynchronize(),"CGPolakRibiereMethod.cu","CGPolakRibiereMethod.cu");
				std::cout<<"Beta serial "<<beta<<endl;
				std::cout<<"newGrad "<<newGrad<<" prevGrad "<<prevGrad<<" curGrad "<<curGrad<<endl;
#endif
#ifdef CUDATIMETRACK
				tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDATIMETRACK
				tbegin = chrono::high_resolution_clock::now();
#endif

				//@@@{ STEP 10 vectorized copy
				tbegin = chrono::high_resolution_clock::now();
				Bead::getDbData().forcesAuxP = Bead::getDbData().forcesAux;
				tend = chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_copy2(tend - tbegin);
				CUDAcommon::tmin.copyforces += elapsed_copy2.count();
				//@@@}
#ifdef CUDATIMETRACK
				tend = chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_runs2b(tend - tbegin);
				s2 += elapsed_runs2b.count();
#endif
#ifdef CUDATIMETRACK
				tbegin = chrono::high_resolution_clock::now();
#endif
				//direction reset if not downhill, or no progress made
				//For any iteration "k"
				//"force" are the conjugate gradient direction (dk)
				//"forceAux" are the force/steepest descent direction (-gk)
				//<-gk+1,dk+1> < 0 => gk and dk are at acute angles with one another
				/*Note: Gradient and conjugate direction should be at obtuse angles for effective
				 * descent*/
				//curGrad = newGrad => gk+1.gk+1 and gk.gk are equal. Gradient has not
				// changed in magnitude.
				//Note: -grad E = ForceAux = -gk. Descent direction = dk = Force
#ifdef TRACKDIDNOTMINIMIZE
				vector<floatingpoint>gradlocal;
				gradlocal.push_back(CGMethod::allFDotFA());
				gradlocal.push_back(curGrad);
				gradlocal.push_back(newGrad);
				gradlocal.push_back(prevGrad);

				SysParams::Mininimization().gradientvec.push_back(gradlocal);
#endif

				Ms_issafestate = CGMethod::allFDotFA() <= 0 || areEqual(curGrad, newGrad);
//        ||        abs(prevGrad/newGrad)<0.1;
				if (Ms_issafestate && Ms_isminimizationstate) {
					//The direction is reset of steepest descent direction (-gk).
					shiftGradient(0.0);
					_safeMode = true;
					#ifdef OPTIMOUT
					safestatuscount++;
					#endif
#ifdef EVSALPHA

					cout << "newGrad " << newGrad << " prevGrad " << prevGrad << " curGrad "
					     << curGrad << endl;
					cout << "beta " << beta << " prevbeta " << prevbeta << endl;
					cout << "FDotFA<0 " << (CGMethod::allFDotFA() <= 0)
					     << " curGrad==newGrad "
					     <<
					     areEqual(curGrad, newGrad) << " abs(prevGrad/newGrad)<0.1 "
					     << (abs(prevGrad / newGrad) < 0.1) << endl;
					calculateEvsalpha(FFM, lambda, LAMBDAMAX, allFDotFA());
					cout << endl;
#endif
				}
				//Create back up coordinates to go to in case Energy minimization fails at an
				// undeisrable state.
				if (maxForce < 10 * GRADTOL && numIter > N / 2) {
					copycoordsifminimumE(maxForce);
				}

#ifdef CUDATIMETRACK
				tend = chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_runs3b(tend - tbegin);
				s3 += elapsed_runs3b.count();
#endif
#ifdef CUDATIMETRACK
				tbegin = chrono::high_resolution_clock::now();
#endif
				curGrad = newGrad;
				Ms_isminimizationstate = maxForce > GRADTOL;
#ifdef CUDATIMETRACK
				tend = chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_runs1b(tend - tbegin);
				s1 += elapsed_runs1b.count();
#endif
			}// End minimization
#endif //SERIAL

#ifdef OPTIMOUT
		std::cout << "SERL Total number of iterations " <<_LINESEARCHALGORITHM<<" "<<
		numIter << endl;
#endif
#ifdef CUDATIMETRACK_MACRO
		std::cout<<"SERL Energy time taken (s) "<<CUDAcommon::serltime.TcomputeE<<" for total "
				"iters "<<CUDAcommon::serltime.Ecount<<endl;
		std::cout<<"SERL Force time taken (s) "<<CUDAcommon::serltime.TcomputeF<<endl;
		std::cout<<"SERL Energy time per iter (s/iter) "<<CUDAcommon::serltime.TcomputeE/
														  (floatingpoint(CUDAcommon::serltime.Ecount))
				 <<endl;
		std::cout<<"SERL Split Times Iter "<<numIter<<" "<<s1<<" "<<s2<<" "<<s3<<" "<<s4<<endl;
		std::cout<<"SERL Add "<<s1+s2+s3+s4<<endl;
#endif
#ifdef CUDATIMETRACK
		tendII= chrono::high_resolution_clock::now();
		chrono::duration<floatingpoint> elapsed_runslice2(tendII - tbeginII);
		std::cout<<"Slice time "<<elapsed_runslice2.count()<<endl;
		tbeginII = chrono::high_resolution_clock::now();
#endif
		if (M_ETolstate[0]) {
			cout << endl;

			cout << "WARNING: Minimization exited when Energy Tolerance was reached at N = "
			     << N << " steps." << endl;
			cout << "Maximum force in system = " << maxF() << endl;
			cout << "System energy..." << endl;
			FFM.computeEnergy(Bead::getDbData().coords.data(), true);
		}

		if (numIter >= N) {
#ifdef CUDAACCL


#endif
			//TODO think about integrating CUDA version here
			cout << endl;

			cout << "WARNING: Did not minimize in N = " << N << " steps." << endl;
			cout << "Maximum force in system = " << maxF() << endl;

#ifdef CUDAACCL
			auto cvars = CUDAcommon::getCUDAvars();
			cvars.streamvec.clear();
			CUDAcommon::cudavars = cvars;
#endif
			cout << "System energy..." << endl;
			FFM.computeEnergy(Bead::getDbData().coords.data(), true);
#ifdef CUDAACCL
			for(auto strm:CUDAcommon::getCUDAvars().streamvec)
				CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
#endif
			cout << endl;
			//Copy back coordinates that correspond to minimum energy
			copybackupcoordinates();

			cout << "Culprit ..." << endl;
			auto b = maxBead();
			if (b != nullptr) b->getParent()->printSelf();
		}
		// Reset backup coordinates with minimum Energy
		Bead::getDbData().coords_minE.resize(0);



	#ifdef TRACKDIDNOTMINIMIZE
	if(numIter) {
		auto tempparams = SysParams::Mininimization();

		cout << "Obegin maxForce Lambda Beta SafeModestatus FDotFA curGrad NewGrad prevGrad TotalE Evec (";
		auto interactionnames = FFM.getinteractionnames();
		for (auto x:interactionnames)
			cout << x << ", ";
		cout << ")" << endl;
		for (auto i = 0; i < tempparams.maxF.size()-1; i++) {
			cout << tempparams.maxF[i] << " " << tempparams.Lambda[i] << " " << tempparams
					.beta[i] << " " << tempparams.safeModeORnot[i] <<" ";
			for(auto j:tempparams.gradientvec[i])
				cout<< j <<" ";
			cout << tempparams.TotalE[i] << " ";
			for (auto j:tempparams.Energyvec[i]) {
				cout << j << " ";
			}
			cout << endl;
		}
		cout<<"End maxF "<<tempparams.maxF[tempparams.maxF.size()-1]<<endl;
		cout << "Oend ------------------" << endl;
	}
	#endif

	#ifdef TRACKDIDNOTMINIMIZE
	SysParams::Mininimization().beta.clear();
	SysParams::Mininimization().Lambda.clear();
	SysParams::Mininimization().Energyvec.clear();
	SysParams::Mininimization().TotalE.clear();
	SysParams::Mininimization().maxF.clear();
	SysParams::Mininimization().safeModeORnot.clear();
	SysParams::Mininimization().tempEnergyvec.clear();
	SysParams::Mininimization().gradientvec.clear();
	FFM.computeEnergy(Bead::getDbData().coords.data(), false);

	#endif

#if defined(CROSSCHECK) || defined(CUDAACCL)
    cross_checkclass::Aux=false;
#endif
#ifdef CUDAACCL
    cvars = CUDAcommon::getCUDAvars();
    cvars.streamvec.clear();
    CUDAcommon::cudavars = cvars;
#endif

    result.energiesAfter = FFM.computeEnergyHRMD(Bead::getDbData().coords.data());

    //final force calculation
    FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forces.data());
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
    for(auto strm:CUDAcommon::getCUDAvars().streamvec)
        CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
#endif
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run2(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run2.count());
    CUDAcommon::cudatime.TcomputeF += elapsed_run2.count();
#endif

#ifdef CUDAACCL
    FFM.CUDAcopyForces(*Msp1, CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_force);
    //copy back forces and calculate load forces in CPU.
#endif
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
#ifdef SERIAL
    Bead::getDbData().forcesAux = Bead::getDbData().forces;

#endif
#ifdef CUDAACCL
    CUDAcommon::handleerror(cudaFreeHost(Mmh_stop));
    CUDAcommon::handleerror(cudaFree(Mmg_stop1));
    CUDAcommon::handleerror(cudaFree(Mmg_stop2));
    CUDAcommon::handleerror(cudaFree(Msg_stop1));
    CUDAcommon::handleerror(cudaFree(Msg_stop2));
    CUDAcommon::handleerror(cudaFree(gpu_GRADTOL));
    CUDAcommon::handleerror(cudaFreeHost(Msh_stop));
    CUDAcommon::handleerror(cudaStreamSynchronize(Ms1));
    CUDAcommon::handleerror(cudaStreamSynchronize(Ms2));
    CUDAcommon::handleerror(cudaStreamSynchronize(Ms3));
    CUDAcommon::handleerror(cudaStreamSynchronize(Ms4));
    if(!(CUDAcommon::getCUDAvars().conservestreams)) {
        CUDAcommon::handleerror(cudaStreamDestroy(Ms1));
        CUDAcommon::handleerror(cudaStreamDestroy(Ms2));
        CUDAcommon::handleerror(cudaStreamDestroy(Ms3));
        CUDAcommon::handleerror(cudaStreamDestroy(Ms4));
        CUDAcommon::handleerror(cudaEventDestroy(Me1));
        CUDAcommon::handleerror(cudaEventDestroy(Me2));
        CUDAcommon::handleerror(cudaEventDestroy(event_safe));
        CUDAcommon::handleerror(cudaEventDestroy(event_dot));
        CUDAcommon::handleerror(cudaStreamDestroy(stream_dotcopy));
        CUDAcommon::handleerror(cudaStreamDestroy(stream_shiftsafe));
    }
#endif
#ifdef SERIAL
    //TODO Comment during SERIAL_CUDACROSSCHECK @{
    delete [] Mc_isminimizationstate;
    delete [] Mc_issafestate;
    delete [] M_ETolstate;
    //@}
#endif
    //TODO make sure it calculates stretchforce in CUDA.
#ifdef CUDAACCL
    FFM.assignallforcemags();
#endif
    #ifdef OPTIMOUT
		cout<<"Safestatuscount "<<safestatuscount<<endl;
	#endif

    //@ STEP 11 END MINIMIZATION
	tbegin = chrono::high_resolution_clock::now();
    endMinimization();


#ifdef CUDATIMETRACK
    tendII= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runslice4(tendII - tbeginII);
    std::cout<<"Slice time "<<elapsed_runslice4.count()<<endl;
    tbeginII = chrono::high_resolution_clock::now();
#endif
    FFM.computeLoadForces();

	#ifdef OPTIMOUT
	std::cout<<"End Minimization************"<<endl;
	cout<<"Time taken for lambda "<<lambdatime<<endl;
	std::cout << "----------------------------------------" << endl;
	#endif

    // compute the Hessian matrix at this point if the feature is enabled
    if(SysParams::Mechanics().hessTracking){
        int total_DOF = Bead::getDbData().coords.size_raw();
        FFM.computeHessian(Bead::getDbData().coords.data(), Bead::getDbData().forcesAux.data(), total_DOF, SysParams::Mechanics().hessDelta);
    }


    FFM.cleanupAllForceFields();
	tend = chrono::high_resolution_clock::now();
	chrono::duration<floatingpoint> elapsed_end(tend - tbegin);
	CUDAcommon::tmin.endminimization+= elapsed_end.count();
	//@} END MINIMIZTION
#ifdef DETAILEDOUTPUT
    std::cout<<"printing beads & forces"<<endl;
    for(auto b:Bead::getBeads()){
        std::cout<<b->getId()<<" "<< b->coordinate() <<" "
                ""<<b->force() <<endl;
    }
    std::cout<<"printed beads & forces"<<endl;
#endif

#ifdef CUDATIMETRACK
    tendII= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runslice3(tendII - tbeginII);
    std::cout<<"Slice time "<<elapsed_runslice3.count()<<endl;
    tendTot= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runtot(tendTot - tbeginTot);
    std::cout<<"Total Minimization time "<<elapsed_runtot.count()<<endl;
#endif

    return result;
}

void PolakRibiere::calculateEvsalpha(ForceFieldManager &FFM, floatingpoint lambda,
		floatingpoint LAMBDAMAX, floatingpoint FDotFA){
	cout<<"Printing Evslambda information "<<endl;
	cout<<"chosen lambda "<<lambda<<" prev lambda "<<prevlambda<<endl;
//	for(floatingpoint alpha=0.0;alpha<0.1;alpha=alpha+1e-4    ){
//		cout<<alpha<<" ";
//	}
//	cout<<endl;
	floatingpoint energyzerolambda = FFM.computeEnergy(Bead::getDbData().coords_bckup.data());
	cout<<"Energy zero lambda "<<energyzerolambda<<" "<<endl;
	cout<<"Lambda = [";
	for(floatingpoint alpha=LAMBDAMAX;alpha>=1e-4;alpha=alpha*LAMBDAREDUCE){
		cout<<alpha<<" ";
	}
	cout<<"];"<<endl;
	int count = 0;
	bool exityes = false;
	cout<<"Energy = [";
	for(floatingpoint alpha=LAMBDAMAX;alpha>=1e-4;alpha=alpha*LAMBDAREDUCE){
		//moveBeads
		const std::size_t num = Bead::getDbData().coords_bckup.size_raw();
		Bead::getDbData().coordsStr.resize(num);
		for(size_t i = 0; i < num; ++i)
			Bead::getDbData().coordsStr.value[i] = Bead::getDbData().coords_bckup
					.value[i] + alpha * Bead::getDbData().forces_bckup.value[i];
		floatingpoint energy = FFM.computeEnergy(Bead::getDbData().coordsStr.data());
		cout<<energy<<" ";
		if(count > 10)
			exityes = true;
		count++;
	}
	cout<<"];"<<endl;
	cout<<"Armijo = [";
	for(floatingpoint alpha=LAMBDAMAX;alpha>=1e-4;alpha=alpha*LAMBDAREDUCE){
		cout<<energyzerolambda-BACKTRACKSLOPE * alpha * FDotFA<<" ";
	}
	cout<<"];"<<endl;
	if(exityes)
		exit(EXIT_SUCCESS);
//	exit(EXIT_FAILURE);
}
