function [] = visualizegeometry(a,b)
% This function visualizes the geometry.
% Inputs:
% a,b   geometric parameter of rect. geometry

figure('Name','original geometry','NumberTitle','off');
% x,y are the coordinates of the 4 nodes
x = [0,a,a,0];
y = [0,0,b,b];
patch(x,y,[0.9 0.9 0.9],'EdgeColor','k','Marker','o','MarkerFaceColor','k');
axis equal 
grid on
end