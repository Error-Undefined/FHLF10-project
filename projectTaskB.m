%% Script for solution to task B

clear
close all
clc

%%--Parameters--%%
QUARANTINE = 'no'; %Set to 'yes' to reduce video quality by 25%
load('projectmeshfinest.mat'); %the mesh to load (consisting of variables p, e, t)
iteration_time = 20; %the time in minutes we want to iterate in our solution
dt = 40; %The time step in the iteration, in seconds

%%--Preprocessor--%%

nelm=length(t(1,:)); %number of elements
edof(:,1)=1:nelm ;
edof(:,2:4)=t(1:3,:)';
coord=p' ;
ndof=max(max(t(1:3,:)));
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);
elength = length(e); % number of boundary edges


%-Material parameters-%
Q = 5 * 10^(7); %Heat generation [W/m^3]
Tinf = 18; %Temperature outside the component [C]
T0 = 30; %Initial temperature [C]
ac = 40; %Convection constant [W/(m^2 K)]
ep = 0.01; %pertrution out of the plane [m]

%Thermal conductivities for Ag-epoxy, silicon and copper [W/(m K)]
kAg = 5;
kSi = 149;
kCu = 385;

%Densities for Ag-epoxy, silicon and copper [kg/m^3]
roAg = 2500;
roSi = 2530;
roCu = 8930;

%Specific heat constants for Ag-epoxy, silicon and copper [J/(kg K)]
cpAg = 1000;
cpSi = 703;
cpCu = 386;

%reduce video quality if enabled
if strcmp('yes',QUARANTINE)
    Q = 0.75*Q;
end

% Initialize the stiffness matrix K, the time derivative matrix C
% and force vectors fl & fb
K = zeros(ndof);
C = zeros(ndof);
fl = zeros(ndof, 1);
fb = zeros(ndof, 1);


%%--Calculate the effect of convection--%%

%Initialize the matrix that represents the contribution of convection to
%the stiffness matrix
Kc = zeros(ndof);

ereduced = e([1 2 5],:); % rows 1,2 and 5 in e are of interest
conv_segments = [8 20 21]; %Boundary segments with convection; these are
                           %8, 20 and 21 from pdetool

%Create a vector containing the nodes that have convection from ereduced
convection_nodes = [];
for i  = 1:elength
    if ismember(ereduced(3,i),conv_segments)
        convection_nodes = [convection_nodes ereduced(1:2,i)];
    end
end

%Calculate the length of the edges with convection, using pythagoras'
convection_edge_length = zeros(1,length(convection_nodes));
for i = 1:length(convection_nodes)
	p1 = convection_nodes(1,i);
    p2 = convection_nodes(2,i);
                
    xdif = p(1,p2) - p(1,p1);
    ydif = p(2,p2) - p(2,p1);
        
    convection_edge_length(i) = sqrt(xdif.^(2) + ydif.^(2));
end

% Calculate the areas of the edges with convection.
% We divide by 1000 as the edge length was calculated in millimeters and we
% want to keep everything in meters
boundary_areas = convection_edge_length*ep/1000;


%We create our force boundary; using the fact that qn = ac(T-Tinf) we can
%insert ac*Tinf*area into the force vector
for i = 1:length(convection_nodes)
   points = convection_nodes(1,i);
   fb(points) = fb(points) + ac*Tinf*boundary_areas(i);
end

%For each node that should have convection, we add a contribution of
%ac*area to the stiffness matrix.
for i = 1:length(convection_nodes)
    node = convection_nodes(1,i);
    Kc(node, node) = Kc(node, node) + ac*boundary_areas(i);
end
%%--End of convection effect calculations--%%

%------Solver------%
for elnbr = 1:nelm
    %The mesh consists of 4 zones, where zone 1 is Ag-epoxy, zone 2 & 4 is
    %copper and zone 3 is the heat generating silicon die. We can get the
    %zone in which an element lies from the 4th row in the matrix t
    zone = t(4,elnbr);
    
    %We set the material constant k according to the zone. If the element
    %lies in the silicon part, we also set the variable Ql distinct from 0
    %so the element supplies heat.
    Ql = 0;
    if zone == 3
        k = kSi;
        x = roSi * cpSi;
        Ql = Q;
    elseif zone == 1
        k = kAg;
        x = roAg * cpAg;
    else
        k = kCu;
        x = roCu * cpCu;
    end
    
    %As our Ex and Ey matrices give the coordinates in mm, we have to
    %divide these by 1000 to have SI units. We calculate the element K and
    %C matrices here.
    [Ke, fe] = flw2te(Ex(elnbr,:)./1000,Ey(elnbr,:)./1000, ep, k*eye(2), Ql);
    Ce = plantml(Ex(elnbr,:)./1000, Ey(elnbr,:)./1000, x);
   
    %Instead of using assem, we use these 4 lines with the same
    %behaviour and much better performance
    index = edof(elnbr,2:end);
    K(index,index) = K(index,index)+Ke;
    fl(index) = fl(index) + fe;
    C(index, index) = C(index, index) + Ce;
end

%Our stiffness matrix becomes the sum of the contribution from convection
%and the matrix obtained from the mesh
Ktilde = K + Kc;

%The force vector is the sum of load forces and boundary forces
f = fl + fb;

%Our nodal temperatures in a are initialized as a column matrix of T0
a = T0 * ones(ndof, 1);

%Given the parameters dt and iteration_time, we iterate the solution using
%the function femImplicitEuler; see the function for a description on how
%it works
for i = 0:dt:iteration_time*60
    a = femImplicitEuler(C, dt, f, Ktilde, a);
end

%------postprocessor------%

%We extract the max temperature from the nodal temperatures we obtained
%and use this to dynamically create the plot title
max_temp = max(a);
titleFormat = 'Temperature distribution after %d minutes;\n max. temperature is %4.2f°C';
titleString = sprintf(titleFormat, iteration_time,max_temp);

%Extract the element displacements
ed = extract(edof, a);

%Plot the results.
figure(1)
axis equal

%We can here use the symmetry to just invert the Ex-matrix and obtain a
%full plot.
patch(Ex'./1000,Ey'./1000,ed','EdgeColor','None');
hold on
patch(-Ex'./1000,Ey'./1000,ed','EdgeColor','None');

colormap hot
cbar = colorbar('southoutside');
xlabel(cbar, '°C')

title(titleString)
ylabel('Height [m]')
xlabel('Width [m]')
