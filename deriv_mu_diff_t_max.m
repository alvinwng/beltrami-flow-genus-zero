function [mu_t,mu_tt,t_max]=deriv_mu_diff_t_max(pt2,tri,r_vec,s_vec,r_dot,s_dot,outward)
% Alvin Wong, April 10, 2013
% Find the first and second derivatives of the Beltrami differential given
% the current one and the direction of flow (r_dot,s_dot). The Beltrami
% differential and its derivatives are computed as a number on each face.
% Also find the smallest timestep when degeneracy occurs for some face.
I=sqrt(-1);
m=size(tri,1);
%[ux,uy,vx,vy]=deriv_3d_map(pt1,pt2,tri,outward);
ux=ones(m,1);
uy=zeros(m,1);
vx=zeros(m,1);
vy=ones(m,1);
[uxp,uyp,vxp,vyp]=deriv_flow(pt2,tri,r_vec,s_vec,r_dot,s_dot,outward);
% mu=tmp2./tmp1;    % verified, gives the same result as compute_beltrami_differential
% [tux,tuy,tvx,tvy]=deriv_3d_map(pt2,pt2+0.0001*(r_vec.*repmat(r_dot,1,3)+s_vec.*repmat(s_dot,1,3)),tri,outward);
% tmp=(tvy-vy)*10000;
tmp1=ux+I*vx-I*uy+vy;
tmp2=ux+I*vx+I*uy-vy;
% mu2=dtmp2./dtmp1;
dtmp1=uxp+I*vxp-I*uyp+vyp;
dtmp2=uxp+I*vxp+I*uyp-vyp;
A=tmp1.*dtmp2-tmp2.*dtmp1;
mu_t=A./tmp1.^2;
% mu_tt=-2*A.*tmp1.*dtmp1./tmp1.^4;
mu_tt=-2*A.*dtmp1./tmp1.^3;  % simplified from the above line
% now determine max_t from the Jacobian: keep at^2+bt+c > 0
a=uxp.*vyp-uyp.*vxp;
b=ux.*vyp+vy.*uxp-vx.*uyp-uy.*vxp;
c=ux.*vy-uy.*vx;	% definitely positive if nothing goes wrong
tmp=b.^2-4*a.*c;
r1=(-b+sqrt(tmp))./a/2;
r2=(-b-sqrt(tmp))./a/2;
t_max=ones(m,1)*Inf;
t_max(imag(r1)==0 & r1>0)=min(t_max(imag(r1)==0 & r1>0),r1(imag(r1)==0 & r1>0));
t_max(imag(r2)==0 & r2>0)=min(t_max(imag(r2)==0 & r2>0),r2(imag(r2)==0 & r2>0));
% t_max=min(t_max);
end