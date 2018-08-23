function calculateRoC(N,outputfile)
%Calculates radius of curvature from snapshot
RoC_parted10={};
RoC10={};
Length10={};
time_vector10={};
polarity_profile={};
orientation_profile={};
r=runs(N);
for i=1:N
	i
    path=['/lustre/qni/ring/nu/forminKinetics3x/t9-2-3']
    %path=['/lustre/qni/turnover/t7/t7-10']
    %path=['./']
    r(i)=readsnapshot([path,'/snapshot.traj'],[path,'/plusend.traj'],r(i));
    time_vector10=[time_vector10;{r(i).time_vector}];
end
VMDstylesnap(r,outputfile);
save([outputfile,'.mat']);
end
