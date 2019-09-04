function [error,dim_largesterror]=getError(dict,lookup,jac,N,T1T2scale)
    idx_N = knnsearch(lookup,lookup(end,:),'k',N+1,'Distance','seuclidean','Scale',T1T2scale);
    x_curr=lookup(end,:);
    y_curr=dict(end,:);
    X_neigh=lookup(idx_N(2:end),:)';
    Y_neigh=dict(idx_N(2:end),:)';
    %synthesize neighboring points with derivative
    A=squeeze(jac(end,:,:));
    Y_neigh_syn =  bsxfun(@plus,A*bsxfun(@plus,X_neigh,-x_curr.'),y_curr.');
    
    relError=rssq(Y_neigh_syn-Y_neigh,1);
    [~,worst_neighbor] = max(relError);
    [~,dim_largesterror] = max(abs(X_neigh(:,worst_neighbor)-x_curr.')./x_curr.');
    
    error = norm(Y_neigh_syn(:)-Y_neigh(:))/norm(Y_neigh(:));
    %[A(iGrid,:,:),~,errors(iGrid)] = getLinearMapping(x_curr',X_neigh,Y_neigh,x_curr',y_curr');
end