# Provides class to convert medyan traj files to vmd files

class Converter(object):

    def snapshot2pdb(path, outputfile):
        """Converts snapshot.traj to VMD readable pdb file.

        path - path to snapshot.traj.
        outuputfile - desired outputfile name.
        """

RoC_parted10={};
RoC10={};
Length10={};
time_vector10={};
polarity_profile={};
orientation_profile={};
    r=runs(N);

    for i=1:N
	i
    r(i)=readsnapshot(['snapshot.traj'],r(i));
    time_vector10=[time_vector10;{r(i).time_vector}];
end
VMDstylesnap_polym(r,outputfile);
save([outputfile,'.mat']);
end
