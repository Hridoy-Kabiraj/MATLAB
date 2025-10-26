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
X_values= 1:40; %MFP
E_values= E([1 9 17 19]);
a_values= a([1 9 17 19]);
b_values= b([1 9 17 19]);
c_values= c([1 9 17 19]);
d_values= d([1 9 17 19]);
Xk_values=Xk([1 9 17 19]);
a1_values= a1([1 9 17 19]);
b1_values= b1([1 9 17 19]);
c1_values= c1([1 9 17 19]);
d1_values= d1([1 9 17 19]);
Xk1_values=Xk1([1 9 17 19]);

% The entries of K for the given energies never reaches exactly
% 1(calculated). Thus only the first equation is used
K= (c_values.*(X_values.^a_values)) + ((d_values.*(tanh((X_values./Xk_values)-2)-tanh(-2)))./(1-tanh(-2)));
EBF= 1 + (((b_values-1).*((K.^X_values)-1))./(K-1));
K1= (c1_values.*(X_values.^a1_values)) + ((d1_values.*(tanh((X_values./Xk1_values)-2)-tanh(-2)))./(1-tanh(-2)));
EABF= 1 + (((b1_values-1).*((K1.^X_values)-1))./(K1-1));

figure;hold on;
plot(X_values,EBF,'-r');
plot(X_values,EABF,'-g');
xlabel('Penetration Depth(MFP) ');
ylabel('Build up Factor');
title('B vs X at photon energies .015,0.15,1.5,3.0 ;Red is EBF ');
hold off;
grid on;

