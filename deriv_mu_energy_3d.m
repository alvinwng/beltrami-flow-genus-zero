function [Et,Ett,t_max]=deriv_mu_energy_3d(pt2,tri,r_vec,s_vec,r_dot,s_dot,M,mu,outward,pt1,mu_tgt)
% Alvin Wong, April 10, 2013
% Find the first and second derivatives of the L2-norm of mu given the map
% and the direction of flow (r_dot,s_dot).
% mu_tgt: optional target mu to achieve, if not provided, assume target is
% 0 mu everywhere
[mu_t,mu_tt,t_max]=deriv_mu_diff_t_max(pt2,tri,r_vec,s_vec,r_dot,s_dot,outward);    % mu_t verified correct numerically
% tmp=compute_beltrami_differential(pt2,pt2+0.0001*(r_vec.*repmat(r_dot,1,3)+s_vec.*repmat(s_dot,1,3)),tri,outward)-compute_beltrami_differential(pt2,pt2,tri,outward);

fz=compute_fz_mesh3(pt1,pt2,tri,outward);
foo1=conj(fz)./fz;
foo2=1-abs(mu).^2;
nu_t=foo1.*foo2.*mu_t;    % VERIFIED CORRECT NUMERICALLY
nu_tt=foo1.*foo2.*mu_tt-2*foo1.^2.*foo2.*conj(mu).*mu_t.^2;     % VERIFIED CORRECT NUMERICALLY

% tmp=compute_beltrami_differential(pt1,pt2+0.0001*(r_vec.*repmat(r_dot,1,3)+s_vec.*repmat(s_dot,1,3)),tri,outward)-compute_beltrami_differential(pt1,pt2,tri,outward);

% tmp=compute_beltrami_differential(pt1,pt2-0.0001*(r_vec.*repmat(r_dot,1,3)+s_vec.*repmat(s_dot,1,3)),tri,outward) ...
%     +compute_beltrami_differential(pt1,pt2+0.0001*(r_vec.*repmat(r_dot,1,3)+s_vec.*repmat(s_dot,1,3)),tri,outward) ...
%     -2*compute_beltrami_differential(pt1,pt2,tri,outward);
% tmp=tmp*10000.^2;

if nargin<11
    mu_tgt=zeros(length(mu),1);
end
x=[real(mu-mu_tgt);imag(mu-mu_tgt)];
%x=[real(mu);imag(mu)];
x_t=[real(nu_t);imag(nu_t)];
x_tt=[real(nu_tt);imag(nu_tt)];
tmp=M+M';
Et=x_t'*tmp*x;
Ett=x_tt'*tmp*x+x_t'*tmp*x_t;
end