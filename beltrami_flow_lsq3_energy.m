function [M,N,C]=beltrami_flow_lsq3_energy(pt,tri,outward,u_vec,v_vec,a,b,A)
% Alvin Wong, April 21, 2013
% Compute the least square energy for the Beltrami flow on a 3D mesh.
% Suppose x=[u;v] is the flow, energy = x'*M*x + N*x +C.
if nargin==7
    A=ones(n,1)/n;
end
m=size(tri,1);
[Mux,Muy,Mvx,Mvy]=deriv_flow_matrix(pt,tri,u_vec,v_vec,outward);
% real part
tmp1=(Mux-Mvy)/2;
% imaginary part
tmp2=(Mvx+Muy)/2;
A2=sparse((1:m)',(1:m)',A);
% now expand (tmp1*x-a)'*A2*(tmp1*x-a)+(tmp2*x-b)'*A2*(tmp2*x-b)
M=tmp1'*A2*tmp1+tmp2'*A2*tmp2;
N=-2*a'*A2*tmp1-2*b'*A2*tmp2;
C=a'*A2*a+b'*A2*b;
end