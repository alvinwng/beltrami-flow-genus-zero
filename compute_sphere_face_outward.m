function outward=compute_sphere_face_outward(pt,tri,r)
% Alvin Wong, January 10, 2013
% Compute the outward normal of each face of a mesh of a sphere of radius
% r.

% n=size(tri,1);
% outward=zeros(n,3);
% for i=1:n
%     tmp=(pt(tri(i,1),:)+pt(tri(i,2),:)+pt(tri(i,3),:))/3;
%     outward(i,:)=tmp/norm(tmp);
% end

if nargin==2
    r=1;
end
pt=pt/r;

outward=(pt(tri(:,1),:)+pt(tri(:,2),:)+pt(tri(:,3),:))/3;
tmp=sqrt(outward(:,1).^2+outward(:,2).^2+outward(:,3).^2);
outward=outward./[tmp tmp tmp];
end
