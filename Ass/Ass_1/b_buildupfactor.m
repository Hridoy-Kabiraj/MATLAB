clear; clc; close all;

data= readtable('GP_fitting_parameters.xlsx');
E= data.E;
a= data.a;
b= data.b;
c= data.c;
d= data.d;
Xk=data.Xk;
a1= data.a_1;
b1= data.b_1;
c1= data.c_1;
d1= data.d_1;
Xk1= data.Xk_1;
X_values= [1,5,10,40]; %MFP


% The entries of K for the given energies never reaches exactly
% 1(calculated). Thus only the first equation is used
K= (c.*(X_values.^a)) + ((d.*(tanh((X_values./Xk)-2)-tanh(-2)))./(1-tanh(-2)));
EBF= 1 + (((b-1).*((K.^X_values)-1))./(K-1));
K1= (c1.*(X_values.^a1)) + ((d1.*(tanh((X_values./Xk1)-2)-tanh(-2)))./(1-tanh(-2)));
EABF= 1 + (((b1-1).*((K1.^X_values)-1))./(K1-1));

figure;hold on;
plot(E,EBF,'-r');
plot(E,EABF,'-g');
xlabel('Photon Energy(MeV) ');
ylabel('Build up Factor');
title('B vs E at penetration depth 1,5,10,40 MFP;Red is EBF ');
hold off;
grid on;

