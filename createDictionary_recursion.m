function [dict,lookup,jac,P] = createDictionary_recursion(errorthresh,theta,TR,B1s)
    tic
    P.N=2; %number of neighbors

    T1min=0.1;
    T1max=5;
    T2min=0.020;
    T2max=3;
    limits=[T1min T1max T2min T2max];
    T1T2scale=[(T1max-T1min) (T2max-T2min)];
%     T1T2scale=[(T1max-T1min) (T2max-T2min) 1];
    
    split_dir=bsxfun(@times,[1 0;
                             0 1],T1T2scale);%cut in this direction
    stencil=bsxfun(@times,[1 0;
                           1 1;
                           0 0],T1T2scale);%create two new points      
    
    %% now do the adaptive sampling
    % compute new points recursively based on local errors
    lookup=[T1min T2min];
    lookup=round(lookup*1E+6)/1E+6;%stay on grid    

    refine = refine_and_check(lookup,[1 1]);  
    refine = round(refine*1E+6)/1E+6;
    
    new=unique([refine;lookup],'rows');
    [dict,lookup,jac]=HSFP_signal(new,theta,B1s,TR,limits);

    dict(:,1)=[]; %no data from TR/2 interval measured
    jac(:,1,:)=[];%no data from TR/2 interval measured    
    
    fprintf('dict size: %2.2d\n',size(dict,1));
    fprintf('generation time: %2.2f sec\n',toc);
    
    function refine = refine_and_check(x_curr,depth)
        refine=[];
        stencil_scaled=bsxfun(@times,stencil,(1./(2.^(depth-1))));
        new_recur=bsxfun(@plus,x_curr,stencil_scaled);
        idx_inaccess=new_recur(:,1)<new_recur(:,2);
        new_recur(idx_inaccess,2)=new_recur(idx_inaccess,1);
        
        [y_new_recur,lookup_new_recur,jac_new_recur]=HSFP_signal(new_recur,theta,B1s,TR,limits);

        if size(lookup_new_recur,1)>2
            [error,dim_largesterror]=getError(y_new_recur,lookup_new_recur,jac_new_recur,P.N,T1T2scale);
            dimvec=([1 2]==dim_largesterror);
            refine = lookup_new_recur;
            if ~(isempty(error) || any(isnan(error)) || (error<errorthresh))
                split_coor=bsxfun(@plus,x_curr,2.^-depth(dim_largesterror)*split_dir(dim_largesterror,:));
                refine = [refine;...
                          refine_and_check(split_coor,depth+dimvec);...
                          refine_and_check(x_curr,depth+dimvec)];             
            end

        end
    end
    
end