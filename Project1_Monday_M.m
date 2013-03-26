%% Project #1
% Authors: Adedayo Lawal and Blake Levy
clc;clear;
%% Set up Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set up field characteristics%%%%%%%%%%%%%%%%
MOVIE = VideoWriter('PEC_BC_Electric_source.avi');
c = 299792458; % speed of light in free space
mu = (4*pi)*1e-7; % permiability of free space
sigma_x = 11; % conductivity for PML region X (Y)-direction
sigma_y = 11; % conductivity for PML region Y (Z)-direction
k = 1;
m = pi; % m-th degree PML grading
epsilon = 1/(mu*c^2); % permitivity of free space
e_top = epsilon; % relative permitivity of top slab (free space)
e_bottom = 4*epsilon; % relative permitivity of bottom slab
c_bottom = 1/sqrt(e_bottom*mu);% relative wave speed in bottom slab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set up gaussian pulse %%%%%%%%%%%%%%%%%%%%%%
f = 1e6; % Center Frequency  of Gaussian Pulse = 1 MHz
w = 2*pi*f; % Angular Frequency Omega
sigma = 3/w; % Time and Bandwidth of the Gaussian Pulse    
t_d = 8*sigma; % Allow for the pulse to be zero at t = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set up geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_top = c/f; % wavelength of top slab
lambda_bottom = c_bottom/f; % wavelength of bottom slab
num_of_wavelengths = 16; % propagate 16 wavelengths in X,Y-direction
b = num_of_wavelengths*lambda_top; % Width of computational domain
a = num_of_wavelengths*lambda_top; % Height of computational domain
num_of_nodes_x = num_of_wavelengths*10;
num_of_nodes_y = num_of_nodes_x;
delx = b/num_of_nodes_x; % space discretization
dely = a/num_of_nodes_y;
delta = dely;
delt = .999*delx/(sqrt(2)*c); % time discretization
pml_offset_x = 0; % additional thickness of boundary in X-direction
pml_offset_y = 0; % additional thickness of boundary in Y-direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set up E,H,T, matrices%%%%%%%%%%%%%%%%%%%%%%
H_xz = zeros(2,num_of_nodes_x + 2*pml_offset_x,num_of_nodes_y + 2*pml_offset_y);
H_xy = zeros(2,num_of_nodes_x + 2*pml_offset_x,num_of_nodes_y + 2*pml_offset_y); % E-field - row one: L+1, row two: L
E_y = zeros(2,num_of_nodes_x + 2*pml_offset_x,num_of_nodes_y + 2*pml_offset_y); % H_y-field - row one: L+1/2, row two: L-1/2
E_z = zeros(2,num_of_nodes_x + 2*pml_offset_x,num_of_nodes_y + 2*pml_offset_y); % H_z-field - row one: L+1/2, row two: L-1/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set up PML Exy and Exz matrices%%%%%%%%%%%%%

Time = 3*num_of_nodes_x; % total time steps
max_color_past = 1; % color map adjust
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set up Source %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source_x = (2*pml_offset_x +num_of_nodes_x)/2; % x-position of source centered on X-axis
source_y = floor((2/3)*(2*pml_offset_y +num_of_nodes_y)); % y-position of source on top slab
E_store = zeros(Time,num_of_nodes_x,num_of_nodes_y);
M = zeros(1,Time); % create source in time-domain
R = zeros(Time,num_of_nodes_y);
f1 = figure(1);
f2 = figure(2);
for L = 2:Time
    M(L) = exp(-(((L-1)*delt-t_d)^2)/(2*sigma^2));    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
open(MOVIE);
for L = 1:floor(.9*Time) % Time March

    tic; % start timer
%     for j = 1:num_of_nodes_y + 2*pml_offset_y % Z-direction (up/down)
%         for i = 1:num_of_nodes_x + 2*pml_offset_x % Y- direction (left/right)
    for i = 1:num_of_nodes_y + 2*pml_offset_y % Z-direction (up/down)
        for j = 1:num_of_nodes_x + 2*pml_offset_x % Y- direction (left/right)
%%%%%%%%%%%%% Boundary region of computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ((i >= num_of_nodes_x+2*pml_offset_x) ||...
                    (j>=num_of_nodes_y+2*pml_offset_y) || (i == 1)...
                    || (j == 1))
                location = 'boundary';
%%%%%%%%%%%%% PML region of computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PML Left interface %%%%%%%%%%%%%%%%%%%%%%%                
            elseif((i <= pml_offset_x) && ...
                    ((j > pml_offset_y) && (j < pml_offset_y + num_of_nodes_y/2)))              
                location = 'PML_Left_Bottom';
            elseif((i <= pml_offset_x) && ...
                    ((j < pml_offset_y+num_of_nodes_y) && (j >= pml_offset_y + num_of_nodes_y/2)))              
                location = 'PML_Left_Top';                
            elseif((i <= pml_offset_x) && ...
                    (j <= pml_offset_y))              
                location = 'PML_Left_Bottom_Corner';   
            elseif((i <= pml_offset_x) && ...
                    (j >= pml_offset_y+num_of_nodes_y))              
                location = 'PML_Left_Top_Corner';   
%%%%%%%%%%%%%%%%%%%%%%%%% PML Right interface %%%%%%%%%%%%%%%%%%%%%%%                
            elseif((i >= pml_offset_x+num_of_nodes_x) && ...
                    ((j > pml_offset_y) && (j < pml_offset_y + num_of_nodes_y/2)))              
                location = 'PML_Right_Bottom';
            elseif((i >= pml_offset_x+num_of_nodes_x) && ...
                    ((j < pml_offset_y+num_of_nodes_y) && (j >= pml_offset_y + num_of_nodes_y/2)))              
                location = 'PML_Right_Top';                
            elseif((i >= pml_offset_x+num_of_nodes_x) && ...
                    (j <= pml_offset_y))              
                location = 'PML_Right_Bottom_Corner';   
            elseif((i >= pml_offset_x+num_of_nodes_x) && ...
                    (j >= pml_offset_y+num_of_nodes_y))              
                location = 'PML_Right_Top_Corner';     
%%%%%%%%%%%%%%%%%%%%%%%%% PML Bottom interface %%%%%%%%%%%%%%%%%%%%%%%                
            elseif((j <= pml_offset_y) && ...
                    ((i > pml_offset_x) && (i < pml_offset_x + num_of_nodes_x)))              
                location = 'PML_Bottom';
%%%%%%%%%%%%%%%%%%%%%%%%% PML Top interface %%%%%%%%%%%%%%%%%%%%%%%                                
            elseif((j >= pml_offset_y+num_of_nodes_y) && ...
                    ((i > pml_offset_x) && (i < pml_offset_x + num_of_nodes_x)))
                location = 'PML_Top';                              
%%%%%%%%%%%%% interface region of computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            elseif (j == num_of_nodes_y/2 + pml_offset_y)
                location = 'interface';
%%%%%%%%%%%%% source region of computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ((i == any(source_x)) && (j == any(source_y)))
                location = 'source';
%%%%%%%%%%%%% lower region of computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            elseif j < pml_offset_y + num_of_nodes_y/2
                location = 'lower';
%%%%%%%%%%%%% upper region of computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            else
                location = 'upper';
            end
            switch location
                case 'boundary'
                    % PEC condition for boundaries, tangential E-fields are continuous
                    H_xz(:,i,j) = 0;
                    H_xy(:,i,j) = 0;
                    E_y(:,i,j) = 0;
                    E_z(:,i,j) = 0; 
                case 'PML_Left_Bottom'
                    pml_e = e_bottom;                              
                    sigma_mx = k*sigma_x;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = 0; sigma_ey = pml_e*sigma_my/mu; 
                    sigma_my = sigma_my*((pml_offset_x - i)/pml_offset_x)^m;
                    sigma_mx = sigma_mx*((pml_offset_x - i)/pml_offset_x)^m;
                    sigma_ey = sigma_ey*((pml_offset_x - i)/pml_offset_x)^m;
                    sigma_ex = sigma_ex*((pml_offset_x - i)/pml_offset_x)^m;                    
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ...
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);
                case 'PML_Right_Bottom'
                    pml_e = e_bottom;                              
                    sigma_mx = k*sigma_x;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = 0; sigma_ey = pml_e*sigma_my/mu;
                    sigma_my = sigma_my*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;
                    sigma_mx = sigma_mx*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;
                    sigma_ey = sigma_ey*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;
                    sigma_ex = sigma_ex*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;                    
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ...
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);                
                case 'PML_Left_Top'
                    pml_e = e_top;                              
                    sigma_mx = sigma_x;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = 0; sigma_ey = pml_e*sigma_my/mu;                                                        
                    sigma_my = sigma_my*((pml_offset_x - i)/pml_offset_x)^m;
                    sigma_mx = sigma_mx*((pml_offset_x - i)/pml_offset_x)^m;
                    sigma_ey = sigma_ey*((pml_offset_x - i)/pml_offset_x)^m;
                    sigma_ex = sigma_ex*((pml_offset_x - i)/pml_offset_x)^m;                    
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ...
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);
                case 'PML_Right_Top'
                    pml_e = e_top;                              
                    sigma_mx = sigma_x;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = 0; sigma_ey = pml_e*sigma_my/mu;
                    sigma_my = sigma_my*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;
                    sigma_mx = sigma_mx*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;
                    sigma_ey = sigma_ey*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;
                    sigma_ex = sigma_ex*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;                    
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ...
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);                
                case 'PML_Bottom'
                    pml_e = e_bottom;                              
                    sigma_mx = 0; %sigma_x/4;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = sigma_y; sigma_ey = pml_e*sigma_my/mu;
                    sigma_my = sigma_my*((pml_offset_y-j)/pml_offset_y)^m;
                    sigma_mx = sigma_mx*((pml_offset_y-j)/pml_offset_y)^m;
                    sigma_ey = sigma_ey*((pml_offset_y-j)/pml_offset_y)^m;
                    sigma_ex = sigma_ex*((pml_offset_y-j)/pml_offset_y)^m;                    
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ...
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);                      
                case 'PML_Top'
                    pml_e = e_top;                              
                    sigma_mx = 0; %sigma_x;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = sigma_y; sigma_ey = pml_e*sigma_my/mu;
                    sigma_my = sigma_my*((j-(pml_offset_y+num_of_nodes_y))/pml_offset_y)^m;
                    sigma_mx = sigma_mx*((j-(pml_offset_y+num_of_nodes_y))/pml_offset_y)^m;
                    sigma_ey = sigma_ey*((j-(pml_offset_y+num_of_nodes_y))/pml_offset_y)^m;
                    sigma_ex = sigma_ex*((j-(pml_offset_y+num_of_nodes_y))/pml_offset_y)^m;                    
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ...
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);                      
                case 'PML_Left_Bottom_Corner'
                    pml_e = e_bottom;                              
                    sigma_mx = sigma_x;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = sigma_y; sigma_ey = pml_e*sigma_my/mu; 
                    sigma_my = sigma_my*((pml_offset_x - i)/pml_offset_x)^m;
                    sigma_mx = sigma_mx*((pml_offset_x - i)/pml_offset_x)^m;
                    sigma_ey = sigma_ey*((pml_offset_x - i)/pml_offset_x)^m;
                    sigma_ex = sigma_ex*((pml_offset_x - i)/pml_offset_x)^m;                    
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ...
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);                            
                case 'PML_Right_Bottom_Corner'
                    pml_e = e_bottom;                              
                    sigma_mx = sigma_x;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = sigma_y; sigma_ey = pml_e*sigma_my/mu; 
                    sigma_my = sigma_my*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;
                    sigma_mx = sigma_mx*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;
                    sigma_ey = sigma_ey*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;
                    sigma_ex = sigma_ex*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;                    
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ...
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);                
                case 'PML_Left_Top_Corner'
                    pml_e = e_top;                              
                    sigma_mx = sigma_x;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = sigma_y; sigma_ey = pml_e*sigma_my/mu;
                    sigma_my = sigma_my*((pml_offset_x - i)/pml_offset_x)^m;
                    sigma_mx = sigma_mx*((pml_offset_x - i)/pml_offset_x)^m;
                    sigma_ey = sigma_ey*((pml_offset_x - i)/pml_offset_x)^m;
                    sigma_ex = sigma_ex*((pml_offset_x - i)/pml_offset_x)^m;                    
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ...
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);
                case 'PML_Right_Top_Corner'
                    pml_e = e_top;                              
                    sigma_mx = sigma_x;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = sigma_y; sigma_ey = pml_e*sigma_my/mu; 
                    sigma_my = sigma_my*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;
                    sigma_mx = sigma_mx*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;
                    sigma_ey = sigma_ey*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;
                    sigma_ex = sigma_ex*((i-(pml_offset_x + num_of_nodes_x))/pml_offset_x)^m;                     
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ...
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);                 
                case 'interface'
                    pml_e = (e_top+e_bottom)/2;                              
                    sigma_mx = 0;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = 0; sigma_ey = pml_e*sigma_my/mu;                                                        
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ...
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);                        
                case 'source'       
                      E_xy(1,source_x,source_y) = -0.5*(delt/e_top)*M(L);
                      E_xz(1,source_x,source_y) = -0.5*(delt/e_top)*M(L);        
%                       H_y(1,source_x,source_y) = -0.5*(delt/mu)*My(L);
%                       H_z(1,source_x,source_y) = -0.5*(delt/mu)*Mz(L);
                case 'lower'
                    pml_e = e_bottom;                              
                    sigma_mx = 0;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = 0; sigma_ey = pml_e*sigma_my/mu;                                                        
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ... 
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);                        
                case 'upper'
                    pml_e = e_top;                              
                    sigma_mx = 0;
                    sigma_ex = pml_e*sigma_mx/mu;
                    sigma_my = 0; sigma_ey = pml_e*sigma_my/mu;                                                        
                    % Finite Difference Equation (4) from our notes
                    H_z(1,i,j) = ((2*delt)/(delta*(2*mu+sigma_mx*delt)))*...
                      (E_xz(2,i+1,j)-E_xz(2,i,j)+E_xy(2,i+1,j)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_mx*delt))*H_z(2,i,j) - ...
                      ((delt*sigma_mx)/(2*mu+sigma_mx*delt))*H_z(2,i,j);
                    % Finite Difference Equation (3) from our notes
                    H_y(1,i,j) = -1*((2*delt)/(delta*(2*mu+sigma_my*delt)))*...
                      (E_xz(2,i,j+1)-E_xz(2,i,j)+E_xy(2,i,j+1)-E_xy(2,i,j)) + ...
                      (2*mu/(2*mu+sigma_my*delt))*H_y(2,i,j) - ...
                      ((delt*sigma_my)/(2*mu+sigma_my*delt))*H_y(2,i,j);    
                    % Finite Difference Equation (2) from our notes
                    E_xy(1,i,j) = -1*((2*delt)/(delta*(2*pml_e+sigma_ey*delt)))*...
                    (H_y(1,i,j)-H_y(1,i,j-1)) + ...
                    (2*pml_e/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j) - ...
                    ((delt*sigma_ey)/(2*pml_e+sigma_ey*delt))*E_xy(2,i,j);                                   
                    % Finite Difference Equation (1) from our notes
                    E_xz(1,i,j) = ((2*delt)/(delta*(2*pml_e+sigma_ex*delt)))*...
                    (H_z(1,i,j)-H_z(1,i-1,j)) + ...
                    (2*pml_e/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j) - ...
                    ((delt*sigma_ex)/(2*pml_e+sigma_ex*delt))*E_xz(2,i,j);                      
                otherwise

            end 
    % test comment  
        end
    % Hard source condition
%     Hard source condition including PML region
    E_xy(1,source_x,source_y) = -0.5*(delt/e_top)*((delt/e_top)^-1)*M(L);
%     Hard source condition including PML region
    E_xz(1,source_x,source_y) = -0.5*(delt/e_top)*((delt/e_top)^-1)*M(L);     
    % Hard source condition
%     H_z(1,source_x,source_y) = -0.5*(delt/mu)*Mz(L);
%     H_y(1,source_x,source_y) = -0.5*(delt/mu)*My(L);
    end
    % Update the row vectors
    H_y(2,:,:) = H_y(1,:,:);   
    H_z(2,:,:) = H_z(1,:,:);
    E_xy(2,:,:) = E_xy(1,:,:);    
    E_xz(2,:,:) = E_xz(1,:,:);
%     display(toc); % Stop timer
    %E = reshape(E_x(1,:,:),[(num_of_nodes_x + 2*pml_offset_x) (num_of_nodes_y + 2*pml_offset_y)]);
    % Reshape E-field matrix for image output
    E = reshape(E_xy(1,:,:)+E_xz(1,:,:),[(num_of_nodes_x + 2*pml_offset_x) (num_of_nodes_y + 2*pml_offset_y)]);    
    H = reshape(sqrt((H_y(1,:,:).^2)+ (H_z(1,:,:).^2)),[(num_of_nodes_x + 2*pml_offset_x) (num_of_nodes_y + 2*pml_offset_y)]);
    Z = abs(E)./abs(H);
   
    % Reshape E-field matrix for color bar ouput
    color = reshape(E_xy(1,:,:)+E_xz(1,:,:),1,(num_of_nodes_x + 2*pml_offset_x)*(num_of_nodes_y + 2*pml_offset_y));
    max_color_present = max(abs(color));
%     if max_color_present > max_color_past
%         c_max = max_color_present;
%         max_color_past = max_color_present;        
%     else
%         c_max = max_color_past;
%     end
%     max_color_present = max(abs(E(:)));
%     if max_color_present > max_color_past
%         c_max = max_color_present;
%         max_color_past = max_color_present;
%     else
%         c_max = max_color_past;
%     end
    
    
    
    % E-field computational domain
    E_comp = E((pml_offset_x+1:(num_of_nodes_x + pml_offset_y)),(pml_offset_y+1:(num_of_nodes_y + pml_offset_y)));
    H_comp = H((pml_offset_x+1:(num_of_nodes_x + pml_offset_y)),(pml_offset_y+1:(num_of_nodes_y + pml_offset_y)));
    Z_comp = Z((pml_offset_x+1:(num_of_nodes_x + pml_offset_y)),(pml_offset_y+1:(num_of_nodes_y + pml_offset_y)));
%     if L == floor(2*Time/3)
%         save('fields_PEC.mat', 'E_comp', 'H_comp', 'Z_comp');
%     end
%     R(L,:) = E(pml_offset_x,((pml_offset_x+1):(pml_offset_y+num_of_nodes_y)))...
%         ./E(pml_offset_x+1,((pml_offset_x+1):(pml_offset_y+num_of_nodes_y)));
    H_y_latest = reshape(H_y(1,:,:),[(num_of_nodes_x + 2*pml_offset_x) (num_of_nodes_y + 2*pml_offset_y)]);
    E_store(L,:,:) = E_comp;
    set(0, 'CurrentFigure', f1)
    imagesc(abs(E)/max(abs(E(:))))
    frame = getframe;
    writeVideo(MOVIE,frame);    
%     imagesc(Z);
%     imagesc(abs(E_comp)/max_color_present)
%     imagesc(abs(E), [0 .25])
%     colorbar;
    
    set(0, 'CurrentFigure', f2)    
    plot(1:(num_of_nodes_y+2*pml_offset_y), E(floor(pml_offset_x+num_of_nodes_x/2),:))%,...
%         1:(num_of_nodes_y+2*pml_offset_y), E(floor(pml_offset_x/2),:));
    pause(.1)

end
close(MOVIE);