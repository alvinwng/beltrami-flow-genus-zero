function [x,y]=stereo_proj(v)

% orientation preserving (outward normal on sphere) stereographic
% projection of an array of points on the unit sphere (n x 3 vector) onto
% 2d
%
% (0,0,1) -> (0,0)
% (x,y,0) -> (x,y)
% (0,0,-1) -> infinity

x = v(:,1) ./ (1 + v(:,3));
y = v(:,2) ./ (1 + v(:,3));