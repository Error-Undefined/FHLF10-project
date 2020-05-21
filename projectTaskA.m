%% Script for solution to task a

clear
close all
clc


%--Parameters--%
%change to alter program behaviour
QUARANTINE = 'no'; %Set to 'yes' to reduce video quality by 25%
PLOT_ELEMENTS = 'no'; %Set to 'no' to remove the element outlines in the plot
load('projectmeshfine.mat'); % the mesh to load (consisting of variables p, e, t)

%------preprocessor------%

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
ac = 40; %Convection constant [W/(m^2 K)]
ep = 0.01; %pertrution out of the plane [m]

%Thermal conductivities for Ag-epoxy, silicon and copper [W/(m K)]
kAg = 5;
kSi = 149;
kCu = 385;

%reduce video quality if enabled
if strcmp('yes',QUARANTINE)
    Q = 0.75*Q;
end


%Initialize stiffnes matrix K and force vectors fl & fb
K = zeros(ndof);
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
        Ql = Q;
    elseif zone == 1
        k = kAg;
    else
        k = kCu;
    end
    
    %As our Ex and Ey matrices give the coordinates in mm, we have to
    %divide these by 1000 to have SI units. We calculate the element K
    %matrix Ke and the element force vector fe here.
    [Ke, fe] = flw2te(Ex(elnbr,:)./1000,Ey(elnbr,:)./1000, ep, k*eye(2), Ql);
    
    %Instead of using assem, we can use these three lines with the same
    %behaviour and much better performance
    index = edof(elnbr,2:end);
    K(index,index) = K(index,index)+Ke;
    fl(index) = fl(index) + fe;
end

%Our stiffness matrix becomes the sum of the contribution from convection
%and the matrix obtained from the mesh
Ktilde = K + Kc;

%The force vector is the sum of load forces and boundary forces
f = fl + fb;

%Solve and extract the element displacements; we don't have any boundary
%conditions as we have convection on the boundary
a = solveq(Ktilde,f);
ed = extract(edof, a);


%------postprocessor------%
figure(1)
axis equal

%Plot with or without element outlines, depending on what was put in the
%beginning. Because of symmetry, we can just invert the Ex matrix.
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
cbar = colorbar('southoutside');
xlabel(cbar, '°C')

title('Stationary temperature distribution in the component')
ylabel('Height [m]')
xlabel('Width [m]')

% Save the stationary temperature distribution (a), as well as the mesh 
% (p, e, t) to use in task C.
save('temperatureDistribution.mat','a', 'p', 'e', 't'); 