clc;clear;
E_pml = load('fields_PML.mat','E_comp');
H_pml = load('fields_PML.mat','H_comp');
Z_pml = load('fields_PML.mat','Z_comp');
E_pec = load('fields_PEC.mat','E');
H_pec = load('fields_PEC.mat','H');
Z_pec = load('fields_PEC.mat','Z');


R = abs(E_pml.E_comp)./(abs(E_pec.E));
imagesc(R);
colormap
