function [Cux,Cuy,Cvx,Cvy]=deriv_flow_coeff(pt,tri,r_vec,s_vec,outward)
% Alvin Wong, April 8, 2013
% Compute the coefficients for the vertices on each triangles of the 3D
% mesh to get the first differentials of a function from the mesh.
% u,v,x,y, are local to each triangle. r,s are at each point of the mesh.
m=size(tri,1);
if nargin==5
    [u_vec,v_vec,~,x2,x3,y3]=get_frame_basis_coord(pt,tri,outward);
else
    [u_vec,v_vec,~,x2,x3,y3]=get_frame_basis_coord(pt,tri);
end
r1u=dot(r_vec(tri(:,1),:),u_vec,2);
r2u=dot(r_vec(tri(:,2),:),u_vec,2);
r3u=dot(r_vec(tri(:,3),:),u_vec,2);
s1u=dot(s_vec(tri(:,1),:),u_vec,2);
s2u=dot(s_vec(tri(:,2),:),u_vec,2);
s3u=dot(s_vec(tri(:,3),:),u_vec,2);
r1v=dot(r_vec(tri(:,1),:),v_vec,2);
r2v=dot(r_vec(tri(:,2),:),v_vec,2);
r3v=dot(r_vec(tri(:,3),:),v_vec,2);
s1v=dot(s_vec(tri(:,1),:),v_vec,2);
s2v=dot(s_vec(tri(:,2),:),v_vec,2);
s3v=dot(s_vec(tri(:,3),:),v_vec,2);
tmp=x2.*y3;
Cux=[-r1u.*y3 r2u.*y3 zeros(m,1) -s1u.*y3 s2u.*y3 zeros(m,1)]./repmat(tmp,1,6);
Cuy=[r1u.*x3-r1u.*x2 -r2u.*x3 r3u.*x2 s1u.*x3-s1u.*x2 -s2u.*x3 s3u.*x2]./repmat(tmp,1,6);
Cvx=[-r1v.*y3 r2v.*y3 zeros(m,1) -s1v.*y3 s2v.*y3 zeros(m,1)]./repmat(tmp,1,6);
Cvy=[r1v.*x3-r1v.*x2 -r2v.*x3 r3v.*x2 s1v.*x3-s1v.*x2 -s2v.*x3 s3v.*x2]./repmat(tmp,1,6);
end