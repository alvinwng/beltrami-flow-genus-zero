function [ux,uy,vx,vy]=deriv_flow(pt,tri,r_vec,s_vec,r_dot,s_dot,outward)
% Alvin Wong, April 8, 2013
% Find the first differentials of a 3D mesh map.
if nargin==7
    [Mux,Muy,Mvx,Mvy]=deriv_flow_matrix(pt,tri,r_vec,s_vec,outward);
else
    [Mux,Muy,Mvx,Mvy]=deriv_flow_matrix(pt,tri,r_vec,s_vec);
end
ux=Mux*[r_dot;s_dot];
uy=Muy*[r_dot;s_dot];
vx=Mvx*[r_dot;s_dot];
vy=Mvy*[r_dot;s_dot];
end