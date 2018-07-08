function [u_vec,v_vec]=compute_sphere_uv_vec(pt,r)

% Alvin Wong, January 10, 2013
% Compute an outward/positively oriented basis of the tangent plane at each
% point of a sphere of radius r.

% assume points in pt lie on the unit sphere
% n=size(pt,1);
% for i=1:n
%     x=pt(i,1);
%     y=pt(i,2);
%     z=pt(i,3);
%     if abs(x)>=abs(y) && abs(x)>=abs(z)
%         tmp=sqrt(x^2+y^2);
%         u_vec(i,:)=[y/tmp -x/tmp 0];
%     elseif abs(y)>=abs(z) && abs(y)>=abs(x)
%         tmp=sqrt(y^2+z^2);
%         u_vec(i,:)=[0 z/tmp -y/tmp];
%     else
%         tmp=sqrt(z^2+x^2);
%         u_vec(i,:)=[-z/tmp 0 x/tmp];
%     end
%     v_vec(i,:)=cross([x y z],u_vec(i,:));
% end

if nargin==1
    r=1;
end
pt=pt/r;

x=pt(:,1);
y=pt(:,2);
z=pt(:,3);
I1=find(abs(x)>=abs(y) & abs(x)>=abs(z));
I2=find(abs(y)>=abs(z) & abs(y)>=abs(x));
I3=find(abs(z)>=abs(x) & abs(z)>=abs(y));
tmp=sqrt(x(I1).^2+y(I1).^2);
u_vec(I1,:)=[y(I1)./tmp -x(I1)./tmp zeros(length(I1),1)];
tmp=sqrt(y(I2).^2+z(I2).^2);
u_vec(I2,:)=[zeros(length(I2),1) z(I2)./tmp -y(I2)./tmp];
tmp=sqrt(z(I3).^2+x(I3).^2);
u_vec(I3,:)=[-z(I3)./tmp zeros(length(I3),1) x(I3)./tmp];
v_vec=cross(pt,u_vec);

end
