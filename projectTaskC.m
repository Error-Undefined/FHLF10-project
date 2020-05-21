%% Script for solution to task C

% Important! This script assumes that a file 'temperatureDistribution.mat',
% containing nodal temperatures as well as a mesh exists. This file is
% created when the script 'projectTaskA.m' is executed.

clear
close all
clc

%-Load the given temperature distribution and mesh-%
load('temperatureDistribution.mat');

%--Preprocessor--%

node_temperatures = a;
clear a

nelm=length(t(1,:));
edof(:,1)=1:nelm ;
coord=p';
ndof=max(max(t(1:3,:)));
elength = length(e); % number of boundary edges

% Because each node has 2 DoF, we let the x-DoF be the indicies given in
% the t-matrix. We then add ndof to each index; this ensures that we get an
% unique index for the y-DoF for the elements.

edof(:,[2 4 6]) = t(1:3,:)';%the x-DoF
edof(:,[3 5 7]) = ndof + t(1:3,:)';%the y-Dof

% We now also need to create a DoF matrix that specifies which DoF are
% coupled to which node. This is easily done by using the fact that
% indicies in the interval [1,ndof] are DoF in the x-direction, and
% indicies in the interval [ndof+1,2*ndof] are DoF in the y-direction. 
DoF = [(1:ndof)' (ndof+1:2*ndof)'];

% We have now created a system where each node has 2 DoF, we now need to
% update the ndof variable before extracting the element coordinates.
ndof = 2* ndof;
[Ex,Ey] = coordxtr(edof,coord,DoF,3);

%-Material parameters-%
zDepth = 0.01; %pertrution out of the plane [m]
T0 = 30; %Initial temperature [C]

% Supply to CALFEM that we have plane strain conditions and generate the
% variable element_properties for later use
ptype = 2;
element_properties = [ptype zDepth];

%Young's modulus for Ag-epoxy, silicon and copper [Pa]
EAg = 7*10^9;
ESi = 165*10^9;
ECu = 128*10^9;

%Poissons' ratio for Ag-epoxy, silicon and copper [-]
vAg = 0.3;
vSi = 0.22;
vCu = 0.36;

%Thermal expansion coefficient for Ag-epoxy, silicon and copper [1/K]
aAg = 4 * 10^(-5);
aSi = 2.6 * 10^(-6);
aCu = 17.6 * 10^(-6);

% We also create the material D matrices using hooke in advance
DAg = hooke(ptype, EAg, vAg);
DSi = hooke(ptype, ESi, vSi);
DCu = hooke(ptype, ECu, vCu);

% Initialize the stiffness matrix K and force vector fl (named f)
% The boundary force vector fb is not initialized as we have no outside
% forces.
K = zeros(ndof);
f = zeros(ndof, 1);

%-Calculate boundary conditions-%

% We generate the boundary conditions: we have that the bottom part of the
% 'feet' are fixed in x and y directions. Because of the symmetry
% conditions, any node on the symmetry line is also fixed in the x
% direction.

feet_border = 1; % Boundary at the bottom of the conductive frame, 
                     % found using pdetool
symmetry_border = [14 15 16 17]; % Boundary at the symmetry line,
                                 % found using pdetool
       
ereduced = e([1 2 5],:); % rows 1,2 and 5 in e are of interest
boundary_condition = [0 0]; %We initialize the boundary condition matrix 
                            %with [0 0] so it is not empty.

% There is no 'nice' way of finding these boundary nodes. As a node can
% appear in one or two boundaries depending on its' position,
% which is unknown, we have to iterate over the whole e-matrix, and then 
% for each element there, iterate over the current boundary condition
% matrix, and add the boundary conditions if they have not already been
% added.
for edge_number = 1:elength
    % We extract and save the two nodes as well as the boundary attribute
    % to use for later
    node1 = ereduced(1,edge_number);
    node2 = ereduced(2,edge_number);
    boundary_attribute = ereduced(3,edge_number);
    
    % Check if the edge is part of the bottom boundary
    if ismember(boundary_attribute, feet_border)
        %Check if any of the nodes are already in the boundary condition
        %matrix. If they aren't, we add in their DoF; these are given by
        %their index (x-DoF) and their index + ndof/2 (y-DoF) as we defined
        %our edof matrix that way in the beginning.
        if ~ismember(boundary_condition, node1)
            boundary_condition = [boundary_condition
                                  node1 0
                                  node1+ndof/2 0];
        end
        if ~ismember(boundary_condition, node2)
            boundary_condition = [boundary_condition
                                  node2 0
                                  node2+ndof/2 0];
        end
    end
    
    % Check if the edge is part of the symmetry boundary
    if ismember(boundary_attribute, symmetry_border)
       % The same procedure as for the feet boundary is used, exept that
       % we only set the x boundary condition.
       if ~ismember(boundary_condition, node1)
            boundary_condition = [boundary_condition
                                  node1 0];
        end
        if ~ismember(boundary_condition, node2)
            boundary_condition = [boundary_condition
                                  node2 0];
        end
    end
end

% We finally need to remove the [0 0]-pair that was inserted before we
% calculated the boundary conditions.
boundary_condition = boundary_condition(2:end,:);

%--Solver--%
for elnbr = 1:nelm
    %The mesh consists of 4 zones, where zone 1 is Ag-epoxy, zone 2 & 4 is
    %copper and zone 3 is the heat generating silicon die. We can get the
    %zone in which an element lies from the 4th row in the matrix t    
    zone = t(4, elnbr);
    
    %We now set the D matrix, expansion coefficient and Poissons' ratio
    %depending on which material the element lies in
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
    
    % To calculate our stresses, we only need the 1st, 2nd and 4th
    % rows/colums in our D matrix; see page 256, eq. 13:38 in the book
    Dstress = D([1 2 4],[1 2 4]);
    
    % To get the temperature difference in the element, we take an average
    % of the three nodes. When solving the temperature distribution problem
    % in task A, the nodal indicies are those that we use as degrees of
    % freedom in the x-direction. We can thus use columns 2, 4 and 6 in
    % our edof matrix to get the indicies of the nodes, and then compute
    % the mean temperature delta in the element
    node_indicies = edof(elnbr, [2 4 6]);
    temperatures = node_temperatures(node_indicies); %node_temperatures was loaded from the temperature distribution
    mean_temp_in_element = mean(temperatures - T0);
    
    %We now calculate the element strains and stresses caused by the
    %temperature increase, according to eq. 13:38 and 13:39 in the book
    epsilon0 = (1 + v)*alpha*mean_temp_in_element*[1 1 0]';
    element_stress = Dstress*epsilon0;
    
    % We now use the CALFEM function plantf to calculate the element forces
    fe = plantf(Ex(elnbr,:)./1000,Ey(elnbr,:)./1000,element_properties,element_stress');
    % We also use plante to calculate the element K matrix
    Ke = plante(Ex(elnbr,:)./1000, Ey(elnbr,:)./1000, element_properties, D);
    %When calling these functions, we divide the coordinates by 1000 in
    %order to convert them from millimeters to meters
    
    
    % Instead of using assem, we use the equivalent following three lines
    % of code to build our global stiffness matrix and force vector.
    index = edof(elnbr,2:end);
    K(index,index) = K(index,index)+Ke;
    f(index) = f(index) + fe;
end

% We now solve the problem using the calculated stiffness matrix, force
% vector and boundary conditions, and extract the element displacements.
dof_displacements = solveq(K, f, boundary_condition);
element_displacements = extract(edof, dof_displacements);

%We can now plot the displacement field. For the symmetry, we need to
%change the signs on the x-displacements as we mirror around the y-axis
element_displacements_symmetry = element_displacements;
element_displacements_symmetry(:,[1 3 5]) = -element_displacements(:,[1 3 5]);

figure(1)
eldraw2(Ex./1000, Ey./1000, [2 4 0]);
hold on
eldraw2(-Ex./1000, Ey./1000, [2 4 0]);
hold on
magnification_factor = eldisp2(Ex./1000, Ey./1000, element_displacements, [2 2 0]);
hold on
eldisp2(-Ex./1000, Ey./1000, element_displacements_symmetry, [2 2 0], magnification_factor);

%Dynamically generate a title containing the magnifcation factor
plot_title = sprintf('Displacement field of the component\nDisplacement shown in green\nMagnification factor: %.4f',magnification_factor);
title(plot_title);
xlabel('Width [m]');
ylabel('Height [m]');

%Using the calculated displacement field, we can calculate the
%corresponding element stresses. This has to be done element wise.
%We initialize a matrix to hold sigma_xx, sigma_yy, sigma_zz and sigma_xy 
%for all elements.
element_stresses = zeros(nelm, 4);

for elnbr = 1:nelm
    %We once again have to get the material of the element, to be able to
    %get the right D matrix.
    zone = t(4, elnbr);
    if zone == 3
        D = DSi;
    elseif zone == 1
        D = DAg;
    else
        D = DCu;
    end
    
    %We now use the CALFEM plants to calculate the element stresses for
    %each element.
    stress = plants(Ex(elnbr,:)./1000, Ey(elnbr,:)./1000, element_properties, D, element_displacements(elnbr,:));
    element_stresses(elnbr,:) = stress;
end

%Using the von Mises stress definition given in the assignment, we calculate the von Mises stress per element
sigmaxx = element_stresses(:,1);
sigmayy = element_stresses(:,2);
sigmazz = element_stresses(:,3);
sigmaxy = element_stresses(:,4);
von_mises_stress = sqrt(sigmaxx.^2 + sigmayy.^2 + sigmazz.^2 - sigmayy.*sigmaxx - sigmayy.*sigmazz - sigmazz.*sigmaxx + 3*sigmaxy.^2);

%From this, we can calculate the stress in a node as the average of the
%stresses in the surrounding elements.

node_von_mises_stress = zeros(length(coord),1);
for node_nbr = 1:length(coord)
    %We find the element indices of the elements surrounding a node. 
    [c0, c1] = find(edof(:,[2 4 6])==node_nbr);
    node_von_mises_stress(node_nbr) = sum(von_mises_stress(c0))/length(c0);
end

