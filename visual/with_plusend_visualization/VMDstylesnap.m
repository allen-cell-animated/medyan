function VMDstylesnap(r,outputfile)
%% Eventhough the MATLAB code that gives the class r cannot take more than one type of filament,...
% linker or motor, we create a class to hold more than one type of those for future expansion.
% Edited 10/04/2016

f1=fopen([outputfile,'.pdb'],'w');
count_model=1;
speciesvec=max_copies(1,1,1,1);
chain_limit = 700;

for runs=1
    for snap=1:size(r(runs).s,2)-1
        dummy=r(runs).s(snap).f.coord_cell1;
        %         dummy2=size(r(runs).s(snap).f.coord_cell1,1);
        dummy2=[];
        for it=1:size(dummy,1)
            dummy2=[dummy2,numel(dummy{it})/3];
        end
        if(size(speciesvec.f_max,2)<size(dummy2,2))
            speciesvec.f_max=padarray(speciesvec.f_max,[0 abs(numel(speciesvec.f_max)-numel(dummy2))],'post');
        elseif(size(speciesvec.f_max,2)>size(dummy2,2))
            dummy2=padarray(dummy2,[0 abs(numel(speciesvec.f_max)-numel(dummy2))],'post');
        end
        %showf_max = length(speciesvec.f_max)
        %show_dummy2 = length(dummy2)
        if (length(dummy)>0)
            speciesvec.f_max=max(dummy2,speciesvec.f_max);
            dummy2=size(r(runs).s(snap).l.coord_cell,1);
            speciesvec.l_max(1)=max(dummy2,speciesvec.l_max(1));
            dummy2=size(r(runs).s(snap).m.coord_cell,1);
            speciesvec.m_max(1)=max(dummy2,speciesvec.m_max(1));
            dummy2=size(r(runs).s(snap).b.coord_cell,1);
            speciesvec.b_max(1)=max(dummy2,speciesvec.b_max(1));
        end
    end
    clear dummy dummy2 it;
    check_snap = size(r(runs).s,2)
    if(mod(check_snap,2)) % if it is an odd number
        sizeofruns = (size(r(runs).s,2)-1)/2
    else
        sizeofruns = size(r(runs).s,2)/2
    end
    
    
    for snap=1:sizeofruns
        fprintf(f1,'MODEL     %4i\n',count_model);%can only pring 10k frames.
        count_f=0;
        fils=r(runs).s(snap).f.coord_cell1;
        
        if(size(speciesvec.f_max,2) <= chain_limit)
            for f=1:size(fils,1)
                A=reshape(fils{f},3,[])'./10;% A model that is 10 times smaller.
                b=0.0;
                for i=1:size(A,1)
                    fprintf(f1,'ATOM  %5i  CA  ARG F%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,A(i,1),A(i,2),A(i,3),b);
                    count_f=count_f+1;
                end
                for i=size(A,1)+1:speciesvec.f_max(f)
                    fprintf(f1,'ATOM  %5i  CA  ARG F%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,A(end,1),A(end,2),A(end,3),b);
                    count_f=count_f+1;
                end
                count_f=count_f+1;
            end
            
            for f=size(fils,1)+1:size(speciesvec.f_max,2)
                b=0.0;
                for i=1:speciesvec.f_max(f)
                    fprintf(f1,'ATOM  %5i  CA  ARG F%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,0.0,0.0,0.0,b);
                    count_f=count_f+1;
                end
                count_f=count_f+1;
            end
            
            %        for i=dummy2+1:speciesvec.f_max(1);
            %            fprintf(f1,'ATOM  %5i  CA  ARG F%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,0.0,0.0,0.0,b);
            %            count_f=count_f+1;
            %        end
            fprintf(f1,'TER   %5i      ARG F%4i\n',count_f,count_f);

            if(count_f > 9999)
                disp(size(speciesvec.f_max,2))
                error('one chain is not enough for current speciesvec.f_max, reduce chian_limit')
            end
            % if the number of beads > the maximum allowed in each chain, create new chain N
        elseif((size(speciesvec.f_max,2) > chain_limit) && (size(speciesvec.f_max,2) < 2*chain_limit))
            if(size(fils,1) < chain_limit) % if there less than  9999 beads in this frame
                for f=1:size(fils,1)
                    A=reshape(fils{f},3,[])'./10;% A model that is 10 times smaller.
                    b=0.0;
                    
                    for i=1:size(A,1)
                        fprintf(f1,'ATOM  %5i  CA  ARG F%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,A(i,1),A(i,2),A(i,3),b);
                        count_f=count_f+1;
                    end
                    for i=size(A,1)+1:speciesvec.f_max(f)
                        fprintf(f1,'ATOM  %5i  CA  ARG F%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,A(end,1),A(end,2),A(end,3),b);
                        count_f=count_f+1;
                    end
                    count_f=count_f+1;
                end
                
                
                for f=size(fils,1)+1:chain_limit
                    b=0.0;
                    for i=1:speciesvec.f_max(f)
                        fprintf(f1,'ATOM  %5i  CA  ARG F%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,0.0,0.0,0.0,b);
                        count_f=count_f+1;
                    end
                    count_f=count_f+1;
                end
                count_f_end = count_f;
                fprintf(f1,'TER   %5i      ARG F%4i\n',count_f,count_f);
                
                if(count_f_end > 9999)
                disp(size(speciesvec.f_max,2))
                error('Chain F goes beyond 9999, reduce chain limit')
                end

                for f=(chain_limit + 1):size(speciesvec.f_max,2) % the rest goes to chain N as coordinates = 0.0
                    b=0.0;
                    for i=1:speciesvec.f_max(f)
                        fprintf(f1,'ATOM  %5i  CA  ARG N%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f - count_f_end,count_f - count_f_end,0.0,0.0,0.0,b);
                        count_f=count_f+1;
                    end
                    count_f=count_f+1;
                end
                fprintf(f1,'TER   %5i      ARG N%4i\n',count_f - count_f_end,count_f - count_f_end);
                
                if(count_f - count_f_end > 9999)
                disp(size(speciesvec.f_max,2))
                error('Chain N goes beyond 9999, need to add new chain or increase chain_limit ')
                end
                % if this frame has more than 9999 beads
            else
                % the first 9999 beads go to chain F
                for f=1:chain_limit
                    A=reshape(fils{f},3,[])'./10;% A model that is 10 times smaller.
                    b=0.0;
                    
                    for i=1:size(A,1)
                        fprintf(f1,'ATOM  %5i  CA  ARG F%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,A(i,1),A(i,2),A(i,3),b);
                        count_f=count_f+1;
                    end
                    for i=size(A,1)+1:speciesvec.f_max(f)
                        fprintf(f1,'ATOM  %5i  CA  ARG F%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,A(end,1),A(end,2),A(end,3),b);
                        count_f=count_f+1;
                    end
                    count_f=count_f+1;
                    count_f_end = count_f;
                end
                fprintf(f1,'TER   %5i      ARG F%4i\n',count_f ,count_f);

                if(count_f_end > 9999)
                disp(size(speciesvec.f_max,2))
                error('Chain F goes beyond 9999, need to reduce chain limit')
                end
                % the rest goes to chain N
                for f=(chain_limit + 1):size(fils,1)
                    A=reshape(fils{f},3,[])'./10;% A model that is 10 times smaller.
                    b=0.0;
                    
                    for i=1:size(A,1)
                        fprintf(f1,'ATOM  %5i  CA  ARG N%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f - count_f_end,count_f - count_f_end,A(i,1),A(i,2),A(i,3),b);
                        count_f=count_f+1;
                    end
                    for i=size(A,1)+1:speciesvec.f_max(f)
                        fprintf(f1,'ATOM  %5i  CA  ARG N%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f - count_f_end,count_f - count_f_end,A(end,1),A(end,2),A(end,3),b);
                        count_f=count_f+1;
                    end
                    count_f=count_f+1;
                end
                
                for f=size(fils,1)+1:size(speciesvec.f_max,2) % the rest goes to chain N as coordinates = 0.0
                    b=0.0;
                    for i=1:speciesvec.f_max(f)
                        fprintf(f1,'ATOM  %5i  CA  ARG N%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f - count_f_end,count_f - count_f_end,0.0,0.0,0.0,b);
                        count_f=count_f+1;
                    end
                    count_f=count_f+1;
                end
                fprintf(f1,'TER   %5i      ARG N%4i\n',count_f - count_f_end,count_f - count_f_end);
                
                if(count_f - count_f_end > 9999)
                disp(size(speciesvec.f_max,2))
                error('Chain N goes beyond 9999, need to add new chain or increase chain_limit')
                end
            end
            
        else
            disp(size(speciesvec.f_max,2))
            error('speciesvec.f_max out of bound, need to adjust chain limit or add new chains')
        end
        
        %% LINKER
        count_f=0;
        link=r(runs).s(snap).l.coord_cell./10;
        for i=1:speciesvec.l_max(1)
            if(i<=size(link,1))
                fprintf(f1,'ATOM  %5i  CA  ARG L%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,link(i,1),link(i,2),link(i,3),b);
                count_f=count_f+1;
                fprintf(f1,'ATOM  %5i  CA  ARG L%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,link(i,4),link(i,5),link(i,6),b);
                count_f=count_f+2;
            else
                fprintf(f1,'ATOM  %5i  CA  ARG L%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,0.0,0.0,0.0,b);
                count_f=count_f+1;
                fprintf(f1,'ATOM  %5i  CA  ARG L%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,0.0,0.0,0.0,b);
                count_f=count_f+2;
            end
        end
        fprintf(f1,'TER   %5i      ARG L%4i\n',count_f,count_f);
        clear link;
        
        %% MOTOR
        count_f=0;
        motor=r(runs).s(snap).m.coord_cell./10;
        for i=1:speciesvec.m_max(1)
            if(i<=size(motor,1))
                fprintf(f1,'ATOM  %5i  CA  ARG M%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,motor(i,1),motor(i,2),motor(i,3),b);
                count_f=count_f+1;
                fprintf(f1,'ATOM  %5i  CA  ARG M%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,motor(i,4),motor(i,5),motor(i,6),b);
                count_f=count_f+2;
            else
                fprintf(f1,'ATOM  %5i  CA  ARG M%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,0.0,0.0,0.0,b);
                count_f=count_f+1;
                fprintf(f1,'ATOM  %5i  CA  ARG M%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,0.0,0.0,0.0,b);
                count_f=count_f+2;
            end
        end
        fprintf(f1,'TER   %5i      ARG M%4i\n',count_f,count_f);
        
        %% Brancher
        count_f = 0;
        brancher = r(runs).s(snap).b.coord_cell./10;
        endbr = size(brancher,1);
        disp(speciesvec.b_max(1))
        checkbrlength = length(speciesvec.b_max(1))
        for i=1:speciesvec.b_max(1)
            if(i<=size(brancher,1))
                fprintf(f1,'ATOM  %5i  CA  ARG B%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,brancher(i,1),brancher(i,2),brancher(i,3),b);
                count_f=count_f+1;
            else
                if endbr == 0;
                    fprintf(f1,'ATOM  %5i  CA  ARG B%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,0.0,0.0,0.0,b);
                    count_f=count_f+1;
                else
                    fprintf(f1,'ATOM  %5i  CA  ARG B%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,brancher(endbr,1),brancher(endbr,2),brancher(endbr,3),b);
                    count_f=count_f+1;
                end
            end
        end
        fprintf(f1,'TER   %5i      ARG B%4i\n',count_f,count_f);
        
        %% PLUSEND-PA
        pa=r(runs).s(snap+sizeofruns+1).p.coord_cell./10;
        ca=r(runs).s(snap+sizeofruns+1).c.coord_cell./10;
        endpa=size(pa,1);
        endca=size(ca,1);
        
        count_f=0;
        for i=1:length(speciesvec.f_max)
            if(i<=size(pa,1))
                fprintf(f1,'ATOM  %5i  CA  ARG P%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,pa(i,1),pa(i,2),pa(i,3),b);
                count_f=count_f+1;
            else
                if endpa == 0
                    fprintf(f1,'ATOM  %5i  CA  ARG C%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,0.0,0.0,0.0,b);
                    count_f=count_f+1;
                else
                    fprintf(f1,'ATOM  %5i  CA  ARG P%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,pa(endpa,1),pa(endpa,2),pa(endpa,3),b);
                    %fprintf(f1,'ATOM  %5i  CA  ARG P%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,0.0,0.0,0.0,b);
                    count_f=count_f+1;
                end
            end
        end
        fprintf(f1,'TER   %5i      ARG P%4i\n',count_f,count_f);
        
        %% PLUSEND-CA
        count_f=0;
        for i=1:length(speciesvec.f_max)
            
            if(i<=size(ca,1))
                fprintf(f1,'ATOM  %5i  CA  ARG C%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,ca(i,1),ca(i,2),ca(i,3),b);
                count_f=count_f+1;
            else
                if endca==0
                    %fprintf(f1,'ATOM  %5i  CA  ARG C%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,pa(endpa,1),pa(endpa,2),pa(endpa,3),b);
                    fprintf(f1,'ATOM  %5i  CA  ARG C%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,0.0,0.0,0.0,b);
                    count_f=count_f+1;
                else
                    fprintf(f1,'ATOM  %5i  CA  ARG C%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,count_f,ca(endca,1),ca(endca,2),ca(endca,3),b);
                    count_f=count_f+1;
                end
            end
        end
        fprintf(f1,'TER   %5i      ARG C%4i\n',count_f,count_f);
        fprintf(f1,'ENDMDL\n');
        count_model=count_model+1;
    end
end
fclose(f1);
end