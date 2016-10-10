% Assigments 3, surf plot of the solution u of the poisson eq
% By Alessandro Piccolo and Anton Sjöberg
clear all; close all;

filename = 'output.txt';
delimiterIn = ' ';
[u,delimiterOut] = importdata(filename)

figure
colormap parula
%mesh(u,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud')
%figure
%mesh(u,'FaceLighting','gouraud','LineWidth',0.3)
%colormap(parula(5))
surf(u,'FaceLighting','gouraud','LineWidth',0.3)
set(gca,'FontSize',13, 'FontWeight', 'bold');
title('Serial representation of our approximated solution to f')
xlabel('x')
ylabel('y')
zlabel('z')
