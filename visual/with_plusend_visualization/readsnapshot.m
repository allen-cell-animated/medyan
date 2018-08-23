function r=readsnapshot(filename,filename2,r)
%% READ SNAPSHOT.TRAJ
f1=fopen(filename,'r');
time_step=str2num(fgetl(f1));
       r.time_vector=time_step(2);
       line=fgetl(f1);
       Coord_cell={};
       Coord_myosin_cell=[];
       Coord_linker_cell=[];
       myosin_id=[];
       linker_id=[];
       brancher_id = [];
       Coord_brancher_cell = [];
       s_no=1;
       Coord=[];
       filament_numbeads=[];
       
while(~feof(f1))
%   while(~feof(f1))
    if(length(line)==0)
        time_step=str2num(fgetl(f1));
        r.time_vector=[r.time_vector,time_step(2)];
        s_=snapshot;
        s_.serial=s_no; % frame number;
        s_no=s_no+1;
%       [coord_cell1,coord_cell2]=cluster_filaments(Coord_cell,Coord,filament_numbeads);
        s_.f.coord_cell1=Coord_cell;
        s_.f.coord_cell2={};
        s_.m.coord_cell=Coord_myosin_cell;
        s_.m.id=myosin_id;
        s_.l.id=linker_id;
        s_.b.id=brancher_id;
        s_.l.coord_cell=Coord_linker_cell;
        s_.b.coord_cell = Coord_brancher_cell;
        r=appendsnapshot(r,s_);
        clear s_;
        Coord=[];
        Coord_cell={};
        Coord_myosin_cell=[];
        Coord_linker_cell=[];
        myosin_id=[];
        linker_id=[];
        brancher_id = [];
        Coord_brancher_cell = [];
        filament_numbeads=[];
    elseif(strcmp(line(1),'F')==1)
        dummy=str2num(fgetl(f1));
        Coord_cell=[Coord_cell;{dummy}];
        Coord=[Coord;reshape(dummy,3,[])'];
        filament_numbeads=[filament_numbeads,size(Coord,1)];
        clear dummy;
    elseif(strcmp(line(1),'M')==1)
        Coord_myosin_cell=[Coord_myosin_cell;str2num(fgetl(f1))];
        dummy=strsplit(line,' ');
        myosin_id=[myosin_id,str2double(dummy(2))];       
    elseif(strcmp(line(1),'L')==1)
        Coord_linker_cell=[Coord_linker_cell;str2num(fgetl(f1))];
        dummy=strsplit(line,' ');
        linker_id=[linker_id,str2double(dummy(2))];
    elseif(strcmp(line(1),'B')==1)
        Coord_brancher_cell = [Coord_brancher_cell;str2num(fgetl(f1))];
        dummy=strsplit(line,' ');
        brancher_id = [brancher_id,str2double(dummy(2))];
    end
   line=fgetl(f1);

end
       if(feof(f1))
       s_=snapshot;
       s_.serial=s_no; % frame number;
       s_no=s_no+1;
%       [coord_cell1,coord_cell2]=cluster_filaments(Coord_cell,Coord,filament_numbeads);
       s_.f.coord_cell1=Coord_cell;
       s_.f.coord_cell2={};
       s_.m.coord_cell=Coord_myosin_cell;
       s_.m.id=myosin_id;
       s_.l.id=linker_id;
       s_.l.coord_cell=Coord_linker_cell;
       s_.b.id=brancher_id;
       s_.b.coord_cell=Coord_brancher_cell;
       r=appendsnapshot(r,s_);
       size(r.s)
       else
           r.time_vector=[r.time_vector(1:end-1)];
       end

       clear s_ Coord_cell Coord_myosin_cell myosin_id filament_numbeads Coord_brancher_cell brancher_id;
       
%% READ PLUSEND.TRAJ
fp=fopen(filename2,'r');
time_step2=str2num(fgetl(fp));
       line=fgetl(fp);
       Coord_p_cell=[];
       Coord_c_cell=[];
       pluse_id=[];
       c_id=[];
       p_id=[];
       s_no=1;
       Coord=[];
       emp=0;

while(~feof(fp))
   %while(~feof(fp))
    if(length(line)==0)
       time_step2=str2num(fgetl(fp));
       s_=snapshot;
       s_.serial=s_no; % frame number;
       s_no=s_no+1;
%       [coord_cell1,coord_cell2]=cluster_filaments(Coord_cell,Coord,filament_numbeads);
       s_.p.coord_cell=Coord_p_cell;
       s_.p.id=p_id;
       s_.c.id=c_id;
       s_.c.coord_cell=Coord_c_cell;
       r=appendsnapshot(r,s_);
       clear s_;
       Coord=[];
       Coord_p_cell=[];
       Coord_c_cell=[];
       p_id=[];
       c_id=[];
       
    elseif(strcmp(line(1),' ')==1)
        %in plusend.traj, the end of each frame may contain two lines
       emp=0; 
       
    elseif(strcmp(line(1),'F')==1)
       Coord_cell2=str2num(fgetl(fp)); %store cell coordinates
       dummy=strsplit(line,' ');
       line_plusend = fgetl(fp);
       if(line_plusend(10) == '0') %PA
           Coord_p_cell = [Coord_p_cell;Coord_cell2];
           p_id=[p_id,str2double(dummy(2))]; 
       elseif(line_plusend(10) == '1') %CA
           Coord_c_cell = [Coord_c_cell;Coord_cell2];
           c_id=[c_id,str2double(dummy(2))];    
       end
       clear dummy;
    end
     line=fgetl(fp);
end

       if(feof(fp))
       s_=snapshot;
       s_.serial=s_no; % frame number;
       s_no=s_no+1;
%       [coord_cell1,coord_cell2]=cluster_filaments(Coord_cell,Coord,filament_numbeads);
       s_.p.coord_cell=Coord_p_cell; 
       s_.p.id=p_id;
       s_.c.coord_cell=Coord_c_cell;
       s_.c.id=c_id;
       end

       fclose(f1);
       fclose(fp);
       clear s_ Coord_cell2 Coord_c_cell Coord_p_cell p_id c_id;

end
