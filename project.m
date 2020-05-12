clear
close all
clc

load('projectmesh.mat');

%preprocessor%

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

K = zeros(ndof);
fl = zeros(ndof, 1);
fb = zeros(ndof, 1);


ereduced = e([1 2 5],:);
conv_segments = [8 20 21]; %Boundary segments with convection

qn_segments = [1 18 19 23 24];

bc = [];
flow_elements = [];

for i  =1:elength
    if ismember(ereduced(3,i),conv_segments)
        bc = [bc
              ereduced(1,i) Tinf];
    end
end



%-Solve-%
for elnbr = 1:nelm
    zone = t(4,elnbr);
    if zone == 3
        k = kSi;
        nodes = t(1:3, elnbr);
        xy = coord(nodes,:);
        area = polyarea(xy(:,1),xy(:,2))*10^(-6);
        fl(nodes) =  Q*area*ep;
    elseif zone == 1
        k = kAg;
    else
        k = kCu;
    end
        
    Ke = flw2te(Ex(elnbr,:)./1000,Ey(elnbr,:)./1000,ep, k*eye(2));
    index = edof(elnbr,2:end);
    K(index,index) = K(index,index) + Ke;
end
f = fl + fb;

a = solveq(K,f,bc);
ed = extract(edof, a);
patch(Ex',Ey',ed');
hold on
patch(-Ex',Ey',ed');
colormap hot
colorbar
