classdef max_copies
    properties
        f_max=0;
        l_max=0;
        m_max=0;
        b_max=0;
    end
    methods
         function obj=max_copies(nf,nl,nm,nb) %number of fil, linker, motor, brancher types
             obj.f_max=nf;
             obj.l_max=nl;
             obj.m_max=nm;
             obj.b_max=nb;
         end
    end
end