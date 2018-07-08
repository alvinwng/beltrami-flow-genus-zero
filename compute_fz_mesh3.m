function f_z=compute_fz_mesh3(pt1,pt2,tri,outward)

% Alvin Wong, November 15, 2012
% slightly modified from compute_beltrami_differential.m

% Alvin Wong, October 10, 2012

% Compute the Beltrami differential of a mesh map.
% pt1: n by 3
% pt2: n by 2 or n by 3
% tri: m by 3
% outward: optional, outward normal for each face of mesh 2

n=size(pt1,1);

if size(pt1,2)==2
    pt1=[pt1 zeros(n,1)];
end

if size(pt2,2)==2
    pt2=[pt2 zeros(n,1)];
end

% m=size(tri,1);
% 
% if nargin==3
%     outward=zeros(m,3);
% end
% 
% f_z=zeros(m,1);
% 
% for i=1:m
%     v1=pt1(tri(i,2),:)-pt1(tri(i,1),:);
%     v2=pt1(tri(i,3),:)-pt1(tri(i,1),:);
%     
%     x1=norm(v1);
%     v1=v1/x1;
%     x2=dot(v2,v1);
%     x3=norm(v2-x2*v1);
%     
%     w1=pt2(tri(i,2),:)-pt2(tri(i,1),:);
%     w2=pt2(tri(i,3),:)-pt2(tri(i,1),:);
%     
%     % compute w3 first
%     w3=cross(w1,w2);
%     w3=w3/norm(w3);
%     if dot(w3,outward(i,:))<0
%         w3=-w3;
%     end
%     
%     y1=norm(w1);
%     w1=w1/y1;
%     y2=dot(w2,w1);
%     w2=cross(w3,w1);
%     y3=dot(pt2(tri(i,3),:)-pt2(tri(i,1),:),w2);
%     
%     % now the map becomes:
%     % (x1,0) |-> (y1,0)
%     % (x2,x3) |-> (y2,y3)
%     
%     % f_zbar=(u_x-v_y)/2+i(v_x+u_y)/2
%     % f_z=(u_x+v_y)/2+i(v_x-u_y)/2
%     M=[y1 y2; 0 y3]/[x1 x2; 0 x3];
%     f_z(i)=M(1,1)+M(2,2)+sqrt(-1)*(M(2,1)-M(1,2));
% end

v1=pt1(tri(:,2),:)-pt1(tri(:,1),:);
v2=pt1(tri(:,3),:)-pt1(tri(:,1),:);

x1=sqrt(v1(:,1).^2+v1(:,2).^2+v1(:,3).^2);
v1=v1./[x1 x1 x1];
x2=dot(v2,v1,2);
tmp=v2-[x2 x2 x2].*v1;
x3=sqrt(tmp(:,1).^2+tmp(:,2).^2+tmp(:,3).^2);

w1=pt2(tri(:,2),:)-pt2(tri(:,1),:);
w2=pt2(tri(:,3),:)-pt2(tri(:,1),:);

w3=cross(w1,w2);
tmp=sqrt(w3(:,1).^2+w3(:,2).^2+w3(:,3).^2);
w3=w3./[tmp tmp tmp];
tmp=sign(dot(w3,outward,2));
w3=w3.*[tmp tmp tmp];

y1=sqrt(w1(:,1).^2+w1(:,2).^2+w1(:,3).^2);
w1=w1./[y1 y1 y1];
y2=dot(w2,w1,2);
w2=cross(w3,w1);
y3=dot(pt2(tri(:,3),:)-pt2(tri(:,1),:),w2,2);

M11=y1.*x3;
M12=-y1.*x2+y2.*x1;
% M13=zeros(m,1);
M22=y3.*x1;

f_z=M11+M22-sqrt(-1)*M12;
f_z=f_z./(x1.*x3)/2;

end










