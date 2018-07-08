function [M,C]=beltrami_flow_lsq3_energy_grad(pt,tri,outward,u_vec,v_vec,a,b,A)
% Alvin Wong, April 21, 2013
% Compute the gradient of the least square energy for the Beltrami flow on
% a 3D mesh.
% Suppose x=[u;v] is the flow, gradient = M*x + C.
n=size(tri,1);
if nargin==4
    A=ones(n,1)/n;
end
[M,C,~]=beltrami_flow_lsq3_energy(pt,tri,outward,u_vec,v_vec,a,b,A);
M=2*M;
C=C';
end

function [M,RHS,C]=beltrami_flow_lsq3_energy_grad_original(pt,tri,outward,u_vec,v_vec,a,b,A)
% now pt is n by 3, the flow depends on the u_vec, and a,b are beltrami
% differentials
% energy=[u v]M[u;v]+[u v]RHS+C
% result=0;
n=size(tri,1);
if nargin==7
    A=ones(n,1)/n;
end
m=size(pt,1);
% M=sparse(2*m,2*m);
% M=spalloc(2*m,2*m,36*m);
% RHS=zeros(2*m,1);
% C=0;
% I=zeros(36*n,1);
% J=zeros(36*n,1);
% S=zeros(36*n,1);
% 
% for i=1:n
%     id1=tri(i,1);
%     id2=tri(i,2);
%     id3=tri(i,3);
%     p1=pt(id1,:);
%     p2=pt(id2,:);
%     p3=pt(id3,:);
%     % change pt(id1,:), pt(id2,:), pt(id3,:) into (0,0), (x1,0), (x2,x3)
%     % and compute v1=x-direction, v2=y-direction, v3=normal of the face
%     v1=p2-p1;
%     x1=norm(v1);
%     v1=v1/x1;
%     v2=p3-p1;
%     x2=dot(v2,v1);
%     v3=cross(v1,v2);
%     v3=v3/norm(v3);
%     if dot(v3,outward(i,:))<0
%         v3=-v3;
%     end
%     v2=cross(v3,v1);
%     x3=dot(p3-p1,v2);
%     [tmp1,tmp2]=xy_derivatives([0 0],[x1 0],[x2 x3]);
%     
%     tmp=[dot(u_vec(id1,:),v1) 0 0 dot(v_vec(id1,:),v1) 0 0; ...
%          0 dot(u_vec(id2,:),v1) 0 0 dot(v_vec(id2,:),v1) 0; ...
%          0 0 dot(u_vec(id3,:),v1) 0 0 dot(v_vec(id3,:),v1); ...
%          dot(u_vec(id1,:),v2) 0 0 dot(v_vec(id1,:),v2) 0 0; ...
%          0 dot(u_vec(id2,:),v2) 0 0 dot(v_vec(id2,:),v2) 0; ...
%          0 0 dot(u_vec(id3,:),v2) 0 0 dot(v_vec(id3,:),v2)];
%     tmp_v1x=tmp1*tmp(1:3,:);
%     tmp_v1y=tmp2*tmp(1:3,:);
%     tmp_v2x=tmp1*tmp(4:6,:);
%     tmp_v2y=tmp2*tmp(4:6,:);
%     
% %     v1_x=tmp_v1x*[u(id1);u(id2);u(id3);v(id1);v(id2);v(id3)];
% %     v1_y=tmp_v1y*[u(id1);u(id2);u(id3);v(id1);v(id2);v(id3)];
% %     v2_x=tmp_v2x*[u(id1);u(id2);u(id3);v(id1);v(id2);v(id3)];
% %     v2_y=tmp_v2y*[u(id1);u(id2);u(id3);v(id1);v(id2);v(id3)];
% %     tmp=(u_x-v_y)^2+(v_x+u_y)^2+4*a(i)^2+4*b(i)^2-4*a(i)*(u_x-v_y)-4*b(i)*(v_x+u_y);
%     
%     tmp_M=(tmp_v1x-tmp_v2y)'*(tmp_v1x-tmp_v2y)+(tmp_v2x+tmp_v1y)'*(tmp_v2x+tmp_v1y);
%     tmp_RHS=-4*a(i)*(tmp_v1x-tmp_v2y)'-4*b(i)*(tmp_v2x+tmp_v1y)';
%     
%     % [I((i-1)*36+1:i*36),J((i-1)*36+1:i*36),S((i-1)*36+1:i*36)]=find(tmp_M*A(i));
%     M([id1 id2 id3 m+id1 m+id2 m+id3],[id1 id2 id3 m+id1 m+id2 m+id3])=...
%         M([id1 id2 id3 m+id1 m+id2 m+id3],[id1 id2 id3 m+id1 m+id2 m+id3])+tmp_M*A(i);
%     RHS([id1 id2 id3 m+id1 m+id2 m+id3])=RHS([id1 id2 id3 m+id1 m+id2 m+id3])+tmp_RHS*A(i);
%     C=C+4*a(i)^2*A(i)+4*b(i)^2*A(i);
%     
% %     result=result+tmp*A(i);
% end
% % M=sparse(I,J,S,2*m,2*m);
% M=M+M';
% RHS=-RHS;


% vectorized version
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
tmp=sign(dot(v3,outward,2));
v3=v3.*[tmp tmp tmp];
v2=cross(v3,v1);
x3=dot(p3-p1,v2,2);
%take inverse: M=[x1 0; x2 x3]^{-1}
tmp=x1.*x3;
M11=x3./tmp;
% M12=zeros(n,1);
M21=-x2./tmp;
M22=x1./tmp;
tmp1=[-M11 M11 zeros(n,1)];
tmp2=[-M21-M22 M21 M22];
tmp_v1x=[tmp1(:,1).*dot(u_vec(id1,:),v1,2) tmp1(:,2).*dot(u_vec(id2,:),v1,2) tmp1(:,3).*dot(u_vec(id3,:),v1,2) tmp1(:,1).*dot(v_vec(id1,:),v1,2) tmp1(:,2).*dot(v_vec(id2,:),v1,2) tmp1(:,3).*dot(v_vec(id3,:),v1,2)];
tmp_v1y=[tmp2(:,1).*dot(u_vec(id1,:),v1,2) tmp2(:,2).*dot(u_vec(id2,:),v1,2) tmp2(:,3).*dot(u_vec(id3,:),v1,2) tmp2(:,1).*dot(v_vec(id1,:),v1,2) tmp2(:,2).*dot(v_vec(id2,:),v1,2) tmp2(:,3).*dot(v_vec(id3,:),v1,2)];
tmp_v2x=[tmp1(:,1).*dot(u_vec(id1,:),v2,2) tmp1(:,2).*dot(u_vec(id2,:),v2,2) tmp1(:,3).*dot(u_vec(id3,:),v2,2) tmp1(:,1).*dot(v_vec(id1,:),v2,2) tmp1(:,2).*dot(v_vec(id2,:),v2,2) tmp1(:,3).*dot(v_vec(id3,:),v2,2)];
tmp_v2y=[tmp2(:,1).*dot(u_vec(id1,:),v2,2) tmp2(:,2).*dot(u_vec(id2,:),v2,2) tmp2(:,3).*dot(u_vec(id3,:),v2,2) tmp2(:,1).*dot(v_vec(id1,:),v2,2) tmp2(:,2).*dot(v_vec(id2,:),v2,2) tmp2(:,3).*dot(v_vec(id3,:),v2,2)];

tmp_re=tmp_v1x-tmp_v2y;
tmp_im=tmp_v2x+tmp_v1y;
I=[id1;id1;id1;id1;id1;id1;...
   id2;id2;id2;id2;id2;id2;...
   id3;id3;id3;id3;id3;id3;...
   m+id1;m+id1;m+id1;m+id1;m+id1;m+id1;...
   m+id2;m+id2;m+id2;m+id2;m+id2;m+id2;...
   m+id3;m+id3;m+id3;m+id3;m+id3;m+id3];
J=[id1;id2;id3;m+id1;m+id2;m+id3;...
   id1;id2;id3;m+id1;m+id2;m+id3;...
   id1;id2;id3;m+id1;m+id2;m+id3;...
   id1;id2;id3;m+id1;m+id2;m+id3;...
   id1;id2;id3;m+id1;m+id2;m+id3;...
   id1;id2;id3;m+id1;m+id2;m+id3];
S=[(tmp_re(:,1).*tmp_re(:,1)+tmp_im(:,1).*tmp_im(:,1)).*A;
   (tmp_re(:,1).*tmp_re(:,2)+tmp_im(:,1).*tmp_im(:,2)).*A;
   (tmp_re(:,1).*tmp_re(:,3)+tmp_im(:,1).*tmp_im(:,3)).*A;
   (tmp_re(:,1).*tmp_re(:,4)+tmp_im(:,1).*tmp_im(:,4)).*A;
   (tmp_re(:,1).*tmp_re(:,5)+tmp_im(:,1).*tmp_im(:,5)).*A;
   (tmp_re(:,1).*tmp_re(:,6)+tmp_im(:,1).*tmp_im(:,6)).*A;
   (tmp_re(:,2).*tmp_re(:,1)+tmp_im(:,2).*tmp_im(:,1)).*A;
   (tmp_re(:,2).*tmp_re(:,2)+tmp_im(:,2).*tmp_im(:,2)).*A;
   (tmp_re(:,2).*tmp_re(:,3)+tmp_im(:,2).*tmp_im(:,3)).*A;
   (tmp_re(:,2).*tmp_re(:,4)+tmp_im(:,2).*tmp_im(:,4)).*A;
   (tmp_re(:,2).*tmp_re(:,5)+tmp_im(:,2).*tmp_im(:,5)).*A;
   (tmp_re(:,2).*tmp_re(:,6)+tmp_im(:,2).*tmp_im(:,6)).*A;
   (tmp_re(:,3).*tmp_re(:,1)+tmp_im(:,3).*tmp_im(:,1)).*A;
   (tmp_re(:,3).*tmp_re(:,2)+tmp_im(:,3).*tmp_im(:,2)).*A;
   (tmp_re(:,3).*tmp_re(:,3)+tmp_im(:,3).*tmp_im(:,3)).*A;
   (tmp_re(:,3).*tmp_re(:,4)+tmp_im(:,3).*tmp_im(:,4)).*A;
   (tmp_re(:,3).*tmp_re(:,5)+tmp_im(:,3).*tmp_im(:,5)).*A;
   (tmp_re(:,3).*tmp_re(:,6)+tmp_im(:,3).*tmp_im(:,6)).*A;
   (tmp_re(:,4).*tmp_re(:,1)+tmp_im(:,4).*tmp_im(:,1)).*A;
   (tmp_re(:,4).*tmp_re(:,2)+tmp_im(:,4).*tmp_im(:,2)).*A;
   (tmp_re(:,4).*tmp_re(:,3)+tmp_im(:,4).*tmp_im(:,3)).*A;
   (tmp_re(:,4).*tmp_re(:,4)+tmp_im(:,4).*tmp_im(:,4)).*A;
   (tmp_re(:,4).*tmp_re(:,5)+tmp_im(:,4).*tmp_im(:,5)).*A;
   (tmp_re(:,4).*tmp_re(:,6)+tmp_im(:,4).*tmp_im(:,6)).*A;
   (tmp_re(:,5).*tmp_re(:,1)+tmp_im(:,5).*tmp_im(:,1)).*A;
   (tmp_re(:,5).*tmp_re(:,2)+tmp_im(:,5).*tmp_im(:,2)).*A;
   (tmp_re(:,5).*tmp_re(:,3)+tmp_im(:,5).*tmp_im(:,3)).*A;
   (tmp_re(:,5).*tmp_re(:,4)+tmp_im(:,5).*tmp_im(:,4)).*A;
   (tmp_re(:,5).*tmp_re(:,5)+tmp_im(:,5).*tmp_im(:,5)).*A;
   (tmp_re(:,5).*tmp_re(:,6)+tmp_im(:,5).*tmp_im(:,6)).*A;
   (tmp_re(:,6).*tmp_re(:,1)+tmp_im(:,6).*tmp_im(:,1)).*A;
   (tmp_re(:,6).*tmp_re(:,2)+tmp_im(:,6).*tmp_im(:,2)).*A;
   (tmp_re(:,6).*tmp_re(:,3)+tmp_im(:,6).*tmp_im(:,3)).*A;
   (tmp_re(:,6).*tmp_re(:,4)+tmp_im(:,6).*tmp_im(:,4)).*A;
   (tmp_re(:,6).*tmp_re(:,5)+tmp_im(:,6).*tmp_im(:,5)).*A;
   (tmp_re(:,6).*tmp_re(:,6)+tmp_im(:,6).*tmp_im(:,6)).*A];
M=sparse(I,J,S,2*m,2*m);
M=M+M'; % checked correct, same as non-vectorized version

tmp=(-4*[a a a a a a].*(tmp_v1x-tmp_v2y)-4*[b b b b b b].*(tmp_v2x+tmp_v1y)).*[A A A A A A];
I=[id1;id2;id3;m+id1;m+id2;m+id3];
J=ones(length(I),1);
S=[tmp(:,1);tmp(:,2);tmp(:,3);tmp(:,4);tmp(:,5);tmp(:,6)];
RHS=full(sparse(I,J,S,2*m,1));
RHS=-RHS; % checked correct, same as non-vectorized version

C=4*sum(a.^2.*A+b.^2.*A); % checked correct, same as non-vectorized version

end