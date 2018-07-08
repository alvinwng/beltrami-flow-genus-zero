function v=stereo_proj_inv(x,y)

% orientation preserving (outward normal on sphere) stereographic
% projection from plane to sphere

v(:,1)=2*x./(1+x.^2+y.^2);
v(:,2)=2*y./(1+x.^2+y.^2);
v(:,3)=(1-x.^2-y.^2)./(1+x.^2+y.^2);

if ~isempty(find(x~=x & y~=y))
    v(x~=x & y~=y,:)=[0 0 -1];
end