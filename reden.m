clear; close all; clc;
do_plot = true;
s = 'Mesh: Uniform mesh on Reden testcase';

pv = [  0   0;
        20  0;
        20  5;
        10  5;
        10  0.5;
        0   0.5;
        0   0];

ev = [2 6];

disp(s)
fd = { 'l_dpolygon', [], pv };
fh = @(p) heightfun(p);
[p,tri] = distmesh( fd, fh, 0.055, [-0.1,-0.1; 20.1,5.1], pv, [], 5000);

if( do_plot )
  clf
  patch( 'vertices', p, 'faces', tri, 'facecolor', [.9, .9, .9] )
  title(s)
  axis tight
  axis equal
end

%% export to vtk

vtkwrite('reden.vtk','polydata','triangle',p(:,1),p(:,2),zeros(size(p,1),1),tri);

