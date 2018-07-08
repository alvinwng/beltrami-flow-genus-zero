function [pt,M]=rotate_space(pt,v1,v2)

% Alvin Wong, April 25, 2014
% Rotate the whole space by an orthonormal matrix to position a direction
% vector to another direction.
% Returns the new set of points and the 3 by 3 matrix used.

if nargin==2
    v2=[0,0,1];
end

v1=v1/norm(v1);
v2=v2/norm(v2);

if isequal(v1,v2)
    M=eye(3);
    return;
end

if isequal(v1,-v2)
    if isequal(cross(v1,[0,0,1]),[0,0,0])
        u=[1,0,0];
    else
        u=cross(v1,[0,0,1]);
        u=u/norm(u);
    end
else
    u=cross(v1,v2);
    u=u/norm(u);
end

w1=cross(u,v1);
w1=w1/norm(w1);
w2=cross(u,v2);
w2=w2/norm(w2);

M=[u' v2' w2']/[u' v1' w1'];

pt=[M(1,1)*pt(:,1)+M(1,2)*pt(:,2)+M(1,3)*pt(:,3) M(2,1)*pt(:,1)+M(2,2)*pt(:,2)+M(2,3)*pt(:,3) M(3,1)*pt(:,1)+M(3,2)*pt(:,2)+M(3,3)*pt(:,3)];

end

