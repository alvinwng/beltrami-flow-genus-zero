function [v1,v2,v3,x1,x2,x3]=get_frame_basis_coord(pt,tri,outward)
% Alvin Wong, April 8, 2013
% Get local orthonormal basis spanning R^3 for a triangular mesh in 3D.
% Returns result for each face. All vectors are #face by 1.
% On each face, the vertices of the triangle take coordinates (0,0),
% (x1,0), (x2,x3)
% outward is optional. No flipping to v3 is done if not given.
id1=tri(:,1);
id2=tri(:,2);
id3=tri(:,3);
p1=pt(id1,:);
p2=pt(id2,:);
p3=pt(id3,:);
v1=p2-p1;
x1=sqrt(v1(:,1).^2+v1(:,2).^2+v1(:,3).^2);
v1=v1./[x1 x1 x1];
v2=p3-p1;
x2=dot(v2,v1,2);
v3=cross(v1,v2);
tmp=sqrt(v3(:,1).^2+v3(:,2).^2+v3(:,3).^2);
v3=v3./[tmp tmp tmp];
if nargin==3
    tmp=sign(dot(v3,outward,2));
    v3=v3.*[tmp tmp tmp];
end
v2=cross(v3,v1);
x3=dot(p3-p1,v2,2);
end