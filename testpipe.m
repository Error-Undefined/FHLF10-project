%clear
close all
clc

load('pipe.mat');

nelm=length(t(1,:));
edof(:,1)=1:nelm ;
edof(:,2:4)=t(1:3,:)';
coord=p' ;
ndof=max(max(t(1:3,:)));
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);

k1 = 1;
k2 = 2*k1;
k3 = 15*k1;
k4 = k1;

D1 = k1*ones(2);
D2 = k2*ones(2);
D3 = k3*ones(2);
D4 = k4*ones(2);
D = D1;

K = zeros(nelm);
f = zeros(nelm,1);
ep = 1;

for elnbr = 1:nelm
    temp = t(4,elnbr);
    if temp == 1
        D = D1;
    end
    if temp == 2
        D = D2;
    end
    if temp == 3
        D = D3;
    end
    if temp == 4
        D = D4;
    end
    Ke = flw2te(Ex(elnbr,:),Ey(elnbr,:), ep, D);
    assem(edof(elnbr,:),K,Ke);
end