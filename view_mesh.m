function p = view_mesh(ptVec,trgVec,color)

% nargin = number of input parameters
if nargin==2
    color = [1 1 0];
end

if size(ptVec,2)==2
    ptVec=[ptVec zeros(size(ptVec,1),1)];
end

% new figure
%figure

% Vertices
% Faces
p = patch('Vertices', ptVec, 'Faces', trgVec);
% FaceColor
set(p,'FaceColor',color);
% EdgeColor
set(p,'EdgeColor','k');
% Use OenGL
set(gcf,'renderer','opengl');

% set aspect ratio, [1 1 2] suppresses z
daspect([1 1 1])

% set default 3D view
view(3)

% set axis limits to range of the data
axis tight
% turns of axis
axis off

% make objects dull
material dull

camlight

lighting phong
% lighting gouraud