%% Part stationary

clear
close all
clc

load('projectmesh.mat');

%--Parameters--%
%change to alter program behaviour
QUARANTINE = 'yes'; %Set to 'yes' to reduce video quality
MODE = 'stationary'; %Set to 'stationary' to solve the stationary heat problem
PLOT_ELEMENTS = 'no'; %Set to 'no' to disable plot of elements

%------preprocessor------%

nelm=length(t(1,:));
edof(:,1)=1:nelm ;
edof(:,2:4)=t(1:3,:)';
coord=p' ;
ndof=max(max(t(1:3,:)));
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);

elength = length(e);

%eldraw2(Ex,Ey,[1,4,1]);
%hold on
%eldraw2(-Ex,Ey,[1,4,1]);

%-Material parameters-%
Q = 5 * 10^(7);
T0 = 30;
Tinf = 18;
ac = 40;
ep = 0.01;

kAg = 5;
kSi = 149;
kCu = 385;

%reduce video quality if enabled%
if strcmp('yes',QUARANTINE)
    Q = 0.75*Q;
end

K = zeros(ndof);
fl = zeros(ndof, 1);
fb = zeros(ndof, 1);

%convection anonymous function
qn =@(T) ac*(T-Tinf);

ereduced = e([1 2 5],:);
conv_segments = [8 20 21]; %Boundary segments with convection

qn_segments = [1 18 19 23 24];

bc = [];
flow_elements = [];

if strcmp('stationary', MODE)
    %enforce Tinf boundary condition%
    for i  =1:elength
        if ismember(ereduced(3,i),conv_segments)
           bc = [bc
                ereduced(1,i) Tinf];
        end
    end
else
    
end


%------Solver------%
for elnbr = nelm:-1:1
    zone = t(4,elnbr);
    if zone == 3
        k = kSi;
    elseif zone == 1
        k = kAg;
    else
        k = kCu;
    end
        
    [Ke, fe] = flw2te(Ex(elnbr,:)./1000,Ey(elnbr,:)./1000,ep, k*eye(2), Q);
    index = edof(elnbr,2:end);
    if zone == 3
        fl(index) = fl(index) + fe;
    end
    K(index,index) = K(index,index) + Ke;
end
f = fl + fb;

a = solveq(K,f,bc);
ed = extract(edof, a);


%------postprocessor------%
figure(1)
axis equal

if strcmp('no', PLOT_ELEMENTS)
    patch(Ex'./1000,Ey'./1000,ed','EdgeColor','None');
    hold on
    patch(-Ex'./1000,Ey'./1000,ed','EdgeColor','None');
else
    patch(Ex'./1000,Ey'./1000,ed');
    hold on
    patch(-Ex'./1000,Ey'./1000,ed');
end
colormap hot
colorbar

%% Part convection

clear
close all
clc

load('projectmeshfine.mat');

%--Parameters--%
%change to alter program behaviour
QUARANTINE = 'no'; %Set to 'yes' to reduce video quality
PLOT_ELEMENTS = 'no'; %Set to 'no' to disable plot of elements

%------preprocessor------%

nelm=length(t(1,:));
edof(:,1)=1:nelm ;
edof(:,2:4)=t(1:3,:)';
coord=p' ;
ndof=max(max(t(1:3,:)));
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);

elength = length(e);

%eldraw2(Ex,Ey,[1,4,1]);
%hold on
%eldraw2(-Ex,Ey,[1,4,1]);

%-Material parameters-%
Q = 5 * 10^(7);
T0 = 30;
Tinf = 18;
ac = 40;
ep = 0.01;

kAg = 5;
kSi = 149;
kCu = 385;

%reduce video quality if enabled%
if strcmp('yes',QUARANTINE)
    Q = 0.75*Q;
end

K = zeros(ndof);
fl = zeros(ndof, 1);
fb = zeros(ndof, 1);

%convection anonymous function

ereduced = e([1 2 5],:);
conv_segments = [8 20 21]; %Boundary segments with convection

%qn_segments = [1 18 19 23 24];

%bc = [];
conv_nodes = [];


for i  = 1:elength
    if ismember(ereduced(3,i),conv_segments)
        conv_nodes = [conv_nodes ereduced(1:2,i)];
	end
end
    
cur_length = zeros(1,length(conv_nodes));
for i = 1:length(conv_nodes)
	p1 = conv_nodes(1,i);
    p2 = conv_nodes(2,i);
                
    xdif = p(1,p2) - p(1,p1);
    ydif = p(2,p2) - p(2,p1);
        
    cur_length(i) = sqrt(xdif.^(2) + ydif.^(2));
end

boundary_areas = cur_length*ep/1000;


for i = 1:length(conv_nodes)
   points = conv_nodes(1,i);
   fb(points) = fb(points) + ac*Tinf*boundary_areas(i);
end


%------Solver------%
for elnbr = 1:nelm
    zone = t(4,elnbr);
    Ql = 0;
    if zone == 3
        k = kSi;
        Ql = Q;
    elseif zone == 1
        k = kAg;
    else
        k = kCu;
    end
    
    
        
    [Ke, fe] = flw2te(Ex(elnbr,:)./1000,Ey(elnbr,:)./1000, ep, k*eye(2), Ql);
    index = edof(elnbr,2:end);
    K(index,index) = K(index,index)+Ke;

    fl(index) = fl(index) + fe;
end

Kc = zeros(ndof);


for i = 1:length(conv_nodes)
    node = conv_nodes(1,i);
    Kc(node, node) = Kc(node, node) + ac*boundary_areas(i);
end

Ktilde = K + Kc;

f = fl + fb;

a = solveq(Ktilde,f);

ed = extract(edof, a);


%------postprocessor------%
figure(1)
axis equal

if strcmp('no', PLOT_ELEMENTS)
    patch(Ex'./1000,Ey'./1000,ed','EdgeColor','None');
    hold on
    patch(-Ex'./1000,Ey'./1000,ed','EdgeColor','None');
else
    patch(Ex'./1000,Ey'./1000,ed');
    hold on
    patch(-Ex'./1000,Ey'./1000,ed');
end
colormap hot
colorbar


%% Iterative

clear
close all
clc

load('projectmeshfinest.mat');

nelm=length(t(1,:));
edof(:,1)=1:nelm ;
edof(:,2:4)=t(1:3,:)';
coord=p' ;
ndof=max(max(t(1:3,:)));
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);

elength = length(e);

%eldraw2(Ex,Ey,[1,4,1]);
%hold on
%eldraw2(-Ex,Ey,[1,4,1]);

%-Material parameters-%
Q = 5 * 10^(7);

%Q = 0.75*Q;

T0 = 30;
Tinf = 18;
ac = 40;
ep = 0.01;

kAg = 5;
kSi = 149;
kCu = 385;

roAg = 2500;
roSi = 2530;
roCu = 8930;

cpAg = 1000;
cpSi = 703;
cpCu = 386;

K = zeros(ndof);
C = zeros(ndof);
fl = zeros(ndof, 1);
fb = zeros(ndof, 1);

%convection anonymous function

ereduced = e([1 2 5],:);
conv_segments = [8 20 21]; %Boundary segments with convection

qn_segments = [1 18 19 23 24];

bc = [];
conv_nodes = [];


for i  = 1:elength
    if ismember(ereduced(3,i),conv_segments)
        conv_nodes = [conv_nodes ereduced(1:2,i)];
	end
end
    
cur_length = zeros(1,length(conv_nodes));
for i = 1:length(conv_nodes)
	p1 = conv_nodes(1,i);
    p2 = conv_nodes(2,i);
                
    xdif = p(1,p2) - p(1,p1);
    ydif = p(2,p2) - p(2,p1);
        
    cur_length(i) = sqrt(xdif.^(2) + ydif.^(2));
end

boundary_areas = cur_length*ep/1000;


for i = 1:length(conv_nodes)
   points = conv_nodes(1,i);
   fb(points) = fb(points) + ac*Tinf*boundary_areas(i);
end


%------Solver------%
for elnbr = 1:nelm
    zone = t(4,elnbr);
    Ql = 0;
    if zone == 3
        k = kSi;
        Ql = Q;
        x = roSi * cpSi;
    elseif zone == 1
        k = kAg;
        x = roAg * cpAg;
    else
        k = kCu;
        x = roCu * cpCu;
    end
    
    
        
    [Ke, fe] = flw2te(Ex(elnbr,:)./1000,Ey(elnbr,:)./1000, ep, k*eye(2), Ql);
    index = edof(elnbr,2:end);
    K(index,index) = K(index,index)+Ke;
    fl(index) = fl(index) + fe;
    
    Ce = plantml(Ex(elnbr,:)./1000, Ey(elnbr,:)./1000, x);
    C(index, index) = C(index, index) + Ce;
end

Kc = zeros(ndof);


for i = 1:length(conv_nodes)
    node = conv_nodes(1,i);
    Kc(node, node) = Kc(node, node) + ac*boundary_areas(i);
end

Ktilde = K + Kc;

f = fl + fb;

dt = 20;
a = T0 * ones(ndof, 1);

totaltime = 20; %minutes

for i = 0:dt:totaltime*60
    a = femImplicitEuler(C, dt, f, Ktilde, a);
end

ed = extract(edof, a);

maxTemp = max(a);
fprintf('Max temperature is %f\n', maxTemp);


%------postprocessor------%
figure(1)
axis equal

ed = extract(edof, a);
patch(Ex'./1000,Ey'./1000,ed','EdgeColor','None');
hold on
patch(-Ex'./1000,Ey'./1000,ed','EdgeColor','None');
    
colormap hot
colorbar

%% Temperature stress
clear
close all
clc

%load a;
%load('projectmesh.mat');

load afine
load('projectmesh.mat');

%-Preprocessor-%

nelm=length(t(1,:));
edof(:,1)=1:nelm ;
edof(:,2:4)=t(1:3,:)';
coord=p';
ndof=max(max(t(1:3,:)));

yDof = edof(:,2:4) + ndof;
edof = [edof(:,1) edof(:,2) yDof(:,1) edof(:,3) yDof(:,2) edof(:,4) yDof(:,3)];
     
Dof = [(1:ndof)' ((ndof+1):2*ndof)'];
%Dof = (1:ndof)';

ndof = ndof * 2 ;

[Ex,Ey]=coordxtr(edof,coord,Dof,3);
%eldraw2(Ex,Ey, [1 4 1])

% Variables
EAg = 7*10^9;
ESi = 165*10^9;
ECu = 128*10^9;

vAg = 0.3;
vSi = 0.22;
vCu = 0.36;

aAg = 4 * 10^(-5);
aSi = 2.6 * 10^(-6);
aCu = 17.6 * 10^(-6);

zDepth = 0.01;
T0 = 30;

% Boundary conditions: the feet are fixed to the ground
% the symmetry line shouldn't move in x-direction
ereduced = e([1 2 5],:);
elength = length(ereduced);

feet = 1;
symmetry_line = 14:17;

bc = [0 0];

for i = 1:elength
   if ismember(ereduced(3,i), feet)
       node1 = ereduced(1,i);
       node2 = ereduced(2,i);
       if ~ismember(bc(:,1),node1)
            bc = [bc
                  node1 0
                  node1+ndof/2 0];
       end
       if ~ismember(bc(:,1),node2)
            bc = [bc
                  node2 0
                  node2+ndof/2 0];
       end
       
   elseif ismember(ereduced(3,i), symmetry_line)
       node1 = ereduced(1,i);
       node2 = ereduced(2,i);
       if ~ismember(bc(:,1),node1)
            bc = [bc
                  node1 0];
       end
       if ~ismember(bc(:,1),node2)
            bc = [bc
            node2 0];
       end
    end
end

bc = bc(2:end,:);

%-Create D matrices-%
ptype = 2; %plane strain
DAg = hooke(ptype, EAg, vAg);
DSi = hooke(ptype, ESi, vSi);
DCu = hooke(ptype, ECu, vCu);

K = zeros(ndof);
f = zeros(ndof, 1);

ep =[ptype zDepth];

for elnbr = 1:nelm
    zone = t(4, elnbr);
    if zone == 3
        D = DSi;
        alpha = aSi;
        v = vSi;
    elseif zone == 1
        D = DAg;
        alpha = aAg;
        v = vAg;
    else
        D = DCu;
        alpha = aCu;
        v = vCu;
    end
    
    Dstress = D([1 2 4],[1 2 4]);
    
    tempNodes = edof(elnbr, [2 4 6]);
    temperatures = a(tempNodes);
    tempInElement = mean(temperatures - T0);
    
    epsilon0 = (1 + v)*alpha*tempInElement*[1 1 0]';
    es = -Dstress*epsilon0;
    
    fe = plantf(Ex(elnbr,:)./1000,Ey(elnbr,:)./1000,ep,es');
    
    Ke = plante(Ex(elnbr,:)./1000, Ey(elnbr,:)./1000, ep, D);
    index = edof(elnbr,2:end);
    K(index,index) = K(index,index)+Ke;
    f(index) = f(index) + fe;
end

u = solveq(K,f, bc);
ed = extract(edof, u);

Es = zeros(nelm,length(DSi));

for elnbr = 1:nelm
    zone = t(4, elnbr);
    if zone == 3
        D = DSi;
    elseif zone == 1
        D = DAg;
    else
        D = DCu;
    end
    
    [es, et] = plants(Ex(elnbr,:)./1000, Ey(elnbr,:)./1000, ep, D, ed(elnbr,:));
    Es(elnbr, :) = es;
end

sigmaxx = Es(:,1);
sigmayy = Es(:,2);
sigmazz = Es(:,3);
sigmaxy = Es(:,4);


vonMises = sqrt(sigmaxx.^2 + sigmayy.^2 + sigmazz.^2 -sigmaxx.*sigmayy - sigmaxx.*sigmazz - sigmayy.*sigmazz+3*sigmaxy.^2);


edSym = ed;
edSym(:,[1 3 5]) = -ed(:,[1 3 5]);

eldraw2(Ex./1000, Ey./1000, [1 4 2]);
hold on
eldraw2(-Ex./1000, Ey./1000, [1 4 2]);
hold on
eldisp2(Ex./1000, Ey./1000, ed, [2 2 1]);
hold on
eldisp2(-Ex./1000, Ey./1000, edSym, [2 2 1]);
%set(gca,'Color','k')


figure(2)

edof = edof(:,[1 2 4 6]);

Seff_nod = zeros(1,length(coord));

for i = 1:length(coord)
    [c0, c1] = find(edof(:,2:4)==i);
    Seff_nod(i) = sum(vonMises(c0))/size(c0,1);
end

seff_plot = extract(edof, Seff_nod);

patch(Ex'./1000,Ey'./1000, seff_plot','EdgeColor', 'None');
hold on
patch(-Ex'./1000,Ey'./1000, seff_plot','EdgeColor', 'None');
colormap hot
colorbar

figure(3)
eldraw2(Ex./1000, Ey./1000, [2 2 1]);
hold on
eldraw2(-Ex./1000, Ey./1000, [2 2 1]);
hold on
eliso2(Ex./1000, Ey./1000, seff_plot, 25, [1 4 2])
%hold on
%eliso2(-Ex./1000, Ey./1000, seff_plot, 15, [1 4 2])