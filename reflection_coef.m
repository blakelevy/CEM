clc;clear;
E_pml = load('fields_PML.mat','E_comp');
H_pml = load('fields_PML.mat','H_comp');
Z_pml = load('fields_PML.mat','Z_comp');
E_pec = load('fields_PEC.mat','E_comp');
H_pec = load('fields_PEC.mat','H_comp');
Z_pec = load('fields_PEC.mat','Z_comp');


R = abs(mean2(E_pml.E_comp)./mean2(E_pec.E_comp));
imagesc(R);
colorbar;
