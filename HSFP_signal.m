function [dict,lookup,jac] = HSFP_signal(lookup_grid,theta,B1s,TR,limits)
    %crop to relevant pairs
    T1min=limits(1);
    T1max=limits(2);
    T2min=limits(3);
    T2max=limits(4);
    lookup_grid(lookup_grid(:,1)<T1min,:)=[];
    lookup_grid(lookup_grid(:,1)>T1max,:)=[];
    lookup_grid(lookup_grid(:,2)<T2min,:)=[];
    lookup_grid(lookup_grid(:,2)>T2max,:)=[];
    lookup_grid(lookup_grid(:,1)<lookup_grid(:,2),:)=[];
    NB1s=numel(B1s);
    Ngridpts=size(lookup_grid,1);
    Npulses=numel(theta);
    dict=zeros(Ngridpts,Npulses);
    jac=zeros(Ngridpts,Npulses,2);
    for iGrid = 1:Ngridpts
        T1=lookup_grid(iGrid,1);
        T2=lookup_grid(iGrid,2);
        T_grads=9.6E-3;
        InvEff=0.95;
        bc=1-(1+InvEff)*exp(-T_grads/T1);

        Mxy_t_avg=0;
        DF_avg=0;
        for iB1=1:NB1s
            B1=B1s(iB1);
            [Mxy_t,t,DF] = computeDynamics(T1,T2,TR,bc,B1*theta);
            Mxy_t_avg=Mxy_t_avg+Mxy_t/NB1s;
            DF_avg=DF_avg+DF/NB1s;
        end
        dict(iGrid,:)=Mxy_t_avg;
        jac(iGrid,:,:)=DF_avg(:,2:3);%only T1 and T2
    end

 
    lookup=lookup_grid;
end
