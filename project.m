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

%eldraw2(Ex,Ey,[1,4,1]);
%hold on
%eldraw2(-Ex,Ey,[1,4,1]);


Q = 5 * 10^(7);
ep = 0.01;

kAg = 5;
kSi = 149;
kCu = 385;

K = zeros(ndof);
f = zeros(ndof, 1);

for elnbr = 1:nelm
    zone = t(4,elnbr);
    k = kCu;
    if zone == 2
        k = kSi;
    end
    if zone == 1
        k = kAg;
    end
    
    index = edof(elnbr,2:end);
        
    Ke = flw2te(Ex(elnbr,:),Ey(elnbr,:),ep, k*ones(2));
    
    fe = 0;
    
    K(index,index) = K(index,index) + Ke; 
    f(index) = f(index) + fe;
end

a = solveq(K,f);
ed = extract(edof, a);
patch(Ex',Ey',ed');
colorbar
