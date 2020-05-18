% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
 % 1) Export the required variables from pdetool and create a MATLAB script
 %    to perform operations on these.
 % 2) Define the problem completely using a MATLAB script. See
 %    http://www.mathworks.com/help/pde/examples/index.html for examples
 %    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
pdetool('snapon','on');
set(ax,'DataAspectRatio',[1 2.3999999999999999 1]);
set(ax,'PlotBoxAspectRatio',[2.5 1.6666666666666667 1]);
set(ax,'XLim',[0 5]);
set(ax,'YLim',[0 8]);
set(ax,'XTick',[ 0,...
 0.25,...
 0.5,...
 0.75,...
 1,...
 1.25,...
 1.5,...
 1.75,...
 2,...
 2.25,...
 2.5,...
 2.75,...
 3,...
 3.25,...
 3.5,...
 3.75,...
 4,...
 4.25,...
 4.5,...
 4.75,...
 5,...
]);
set(ax,'YTick',[ 0,...
 0.25,...
 0.5,...
 0.75,...
 1,...
 1.25,...
 1.5,...
 1.75,...
 2,...
 2.25,...
 2.5,...
 2.75,...
 3,...
 3.25,...
 3.5,...
 3.75,...
 4,...
 4.25,...
 4.5,...
 4.75,...
 5,...
 5.25,...
 5.5,...
 5.75,...
 6,...
 6.25,...
 6.5,...
 6.75,...
 7,...
 7.25,...
 7.5,...
 7.75,...
 8,...
]);
pdetool('gridon','on');

% Geometry description:
pdepoly([ 3.5,...
 4,...
 4,...
 3.5,...
 2,...
 2,...
 3.25,...
 3.5,...
],...
[ 0,...
 0,...
 4,...
 4.5,...
 4.5,...
 4,...
 4,...
 3.75,...
],...
 'P1');
pderect([0 4 3 7],'SQ1');
pdepoly([ 4,...
 3,...
 4,...
],...
[ 6,...
 7,...
 7,...
],...
 'P2');
pderect([0 1 4 4.5],'R1');
pderect([0 0.5 5 4.5],'SQ2');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','((P1+SQ1)-P2)+R1+SQ2')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(24,...
'dir',...
1,...
'1',...
'0')
pdesetbd(23,...
'dir',...
1,...
'1',...
'0')
pdesetbd(21,...
'dir',...
1,...
'1',...
'0')
pdesetbd(20,...
'dir',...
1,...
'1',...
'0')
pdesetbd(19,...
'dir',...
1,...
'1',...
'0')
pdesetbd(18,...
'dir',...
1,...
'1',...
'0')
pdesetbd(12,...
'dir',...
1,...
'1',...
'0')
pdesetbd(11,...
'dir',...
1,...
'1',...
'0')
pdesetbd(10,...
'dir',...
1,...
'1',...
'0')
pdesetbd(9,...
'dir',...
1,...
'1',...
'0')
pdesetbd(5,...
'dir',...
1,...
'1',...
'0')
pdesetbd(1,...
'dir',...
1,...
'1',...
'0')

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
setappdata(pde_fig,'MesherVersion','preR2013a');
pdetool('initmesh')
pdetool('jiggle')
pdetool('refine')

% PDE coefficients:
pdeseteq(1,...
'1.0',...
'0.0',...
'10',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0';...
'0.0';...
'10 ';...
'1.0'])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','1074','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');