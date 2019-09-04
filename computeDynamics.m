function [Mxy_t,t,jac,r_t,jac_r,djac_dtheta] = computeDynamics(T1,T2,TR,bc,theta_t)
%compute dynamics in hybrid space, Asslaender 2018, arXiv (eq. #7)
    N=numel(theta_t);
    switch bc
        case 'INV'
            r_0=-1;
        otherwise
            r_0=bc;
    end    
    
    TR_prep=TR*ones(size(theta_t));
    TR_prep(1)=TR_prep(1)/2;
    a_t=exp(-TR_prep.*cumtrapz(sind(theta_t).^2/T2+cosd(theta_t).^2/T1));
    int=TR_prep.*cumtrapz(cosd(theta_t)./a_t);
    r_t=a_t.*(r_0+1/T1*int);
    Mxy_t=r_t.*sind(theta_t);
    t=TR_prep.*(0.5+(0:N-1));
    
    if nargout > 2
        jac=zeros(N,4);
        jac_r=zeros(N,4);
        djac_dtheta=zeros(N,4,N);%last dim is dtheta
        
        jac_r(:,1)=+r_t;

        da_dT1=a_t.*TR_prep.*cumtrapz(cosd(theta_t).^2/T1^2);
        jac_r(:,2)=(da_dT1.*(r_0+1/T1*int)...
            +a_t.*(-1/T1^2*int...
            -1/T1.*TR_prep.*cumtrapz(cosd(theta_t)./a_t.^2.*da_dT1)));
              
        da_dT2=a_t.*TR_prep.*cumtrapz(sind(theta_t).^2/T2^2);
        jac_r(:,3)=(da_dT2.*(r_0+1/T1*int)...
            +a_t.*(-1/T1.*TR_prep.*cumtrapz(cosd(theta_t)./a_t.^2.*da_dT2)));
        
        %assume theta=B1_eff*theta_nom
        da_dB1eff=pi/180*(-a_t.*TR_prep.*cumtrapz(2*sind(theta_t)...
                      .*cosd(theta_t).*theta_t./T2-2.*sind(theta_t).*cosd(theta_t).*theta_t./T1));
        jac_r(:,4)=   (da_dB1eff.*(r_0+1./T1.*int)...
                      +a_t.*(1./T1.*TR_prep.*cumtrapz(-sind(theta_t)...
                      .*theta_t*pi/180./a_t-cosd(theta_t)./a_t.^2.*da_dB1eff)));
        
        jac(:,1:3)= bsxfun(@times,jac_r(:,1:3),sind(theta_t(:)));
        jac(:,4) = jac_r(:,4).*sind(theta_t(:))...
                  +r_t(:).*cosd(theta_t(:)).*theta_t(:).*pi/180;
    end
end