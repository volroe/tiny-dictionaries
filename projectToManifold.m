function [param,x_pr, x_match, y_match, PD]...
    =projectToManifold(x_lin,lookup,dictionary,jacobians)
tic
dictionary_matching=bsxfun(@times,dictionary,1./rssq(dictionary,2));

sprod=abs(dictionary_matching*x_lin');
[~,ip_sorted]=max(sprod,[],1);
t_match=toc;  

x_match=lookup(ip_sorted,:);
y_match=dictionary(ip_sorted,:);
PD=sum(single(x_lin).*dictionary_matching(ip_sorted,:),2)./rssq(dictionary(ip_sorted,:),2);
A=jacobians(ip_sorted,:,:);
x_pr=zeros(size(x_lin));
param=zeros(size(x_lin,1),3);
for i = 1:size(A,1)
    a=squeeze(A(i,:,:));
    a_ext=[a(:,1:2) -x_lin(i,:).'];% as defined in eq. (6)
    pinva_ext=pinv(a_ext);
    param(i,:)=pinva_ext*(a*x_match(i,:).'-y_match(i,:).');
    x_pr(i,:)=a(:,1:2)*([param(i,1:2)].'-x_match(i,1:2).')+y_match(i,:).';
    if param(i,end)~=0
        param(i,end)=1./param(i,end);
        x_pr(i,:)=x_pr(i,:)*param(i,end);
    else
        x_pr(i,:)=0*x_pr(i,:);
    end
end

fprintf('matching time: %2.2f sec\n', t_match);
fprintf('projecting time: %2.2f sec\n',toc);

end