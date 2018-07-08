function [Mux,Muy,Mvx,Mvy]=deriv_flow_matrix(pt,tri,r_vec,s_vec,outward)
% Alvin Wong, April 8, 2013
% Find the matrices representing the first differentials of a flow on a 3D
% mesh.
% pt: n by 3
% tri: m by 3
% Mux,Muy,Mvx,Mvy: m by 2n sparse
m=size(tri,1);
n=size(pt,1);
if nargin==5
    [Cux,Cuy,Cvx,Cvy]=deriv_flow_coeff(pt,tri,r_vec,s_vec,outward);
else
    [Cux,Cuy,Cvx,Cvy]=deriv_flow_coeff(pt,tri,r_vec,s_vec);
end
I=repmat((1:m)',1,6);
J=[tri n+tri];
Mux=sparse(I(:),J(:),Cux(:),m,2*n);
Muy=sparse(I(:),J(:),Cuy(:),m,2*n);
Mvx=sparse(I(:),J(:),Cvx(:),m,2*n);
Mvy=sparse(I(:),J(:),Cvy(:),m,2*n);
end