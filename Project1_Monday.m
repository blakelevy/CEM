%% Project #1
% Authors: Adedayo Lawal and Blake Levy
clc;clear;
%% Set up Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set up field characteristics%%%%%%%%%%%%%%%%
c = 299792458; % speed of light in free space
mu = (4*pi)*1e-7; % permiability of free space
epsilon = 1/(mu*c^2); % permitivity of free space
e_top = epsilon; % relative permitivity of top slab (free space)
e_bottom = 3*epsilon; % relative permitivity of bottom slab
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
num_of_wavelengths = 50; % propagate 120 wavelengths in X,Y-direction
b = num_of_wavelengths*lambda_top; % Width of computational domain
a = num_of_wavelengths*lambda_top; % Height of computational domain
num_of_nodes_x = 1000;
num_of_nodes_y = num_of_nodes_x;
delx = b/num_of_nodes_x; % space discretization
dely = a/num_of_nodes_y;
delta = dely;
delt = delx/(sqrt(2)*c_bottom); % time discretization
pml_offset_x = 0; % additional thickness of boundary in X-direction
pml_offset_y = 0; % additional thickness of boundary in Y-direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set up E,H,T, matrices%%%%%%%%%%%%%%%%%%%%%%
E = zeros(3,num_of_nodes_x,num_of_nodes_y); % E-field - row one: L+1, row two: L, row 3: L-1
H_y = zeros(3,num_of_nodes_x,num_of_nodes_y); % H_y-field - row one: L+1, row two: L, row 3: L-1
H_z = zeros(3,num_of_nodes_x,num_of_nodes_y); % H_z-field - row one: L+1, row two: L, row 3: L-1
Time = 2*num_of_nodes_x; % total time steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for L = 1:Time % Time March
    for j = 1:num_of_nodes_x + 2*pml_offset_y % Z-direction (up/down)
        for i = 1:num_of_nodes_y + 2*pml_offset_x % Y- direction (left/right)
%%%%%%%%%%%%% PML region of computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if((i <= pml_offset_x) || (j <= pml_offset_y) ||...
                    (i >= num_of_nodes_x) || (j >= num_of_nodes_y))              
                location = 'PML';
%%%%%%%%%%%%% interface region of computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            elseif (j == num_of_nodes_y/2 + pml_offset_y)
                location = 'interface';
%%%%%%%%%%%%% source region of computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ((i == source_x) && (j == source_y))
                location = 'source';
%%%%%%%%%%%%% lower region of computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            elseif j < pml_offset_y + num_of_nodes_y/2
                location = 'lower';
%%%%%%%%%%%%% upper region of computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            else
                location = 'upper';
            end
switch location
    case 'pml'
        
        
    case 'interface'
        
        
        
    case 'source'
        
        
        
    case 'lower'
        % Finite Difference Equation (3) from our notes
        H_z(1,i,j) = (delt/mu)*(((E(2,i+1,j) - E(2,i,j))/delta) - (s_bottom...
            *((H_z(1,i+1,j)+H_z(3,i+1,j))/2))) + H_z(3,i-1,j);
        % Finite Difference Equation (2) from our notes
        H_y(1,i,j+1) = -(delt/mu)*(((-E(2,i,j+1)+E(2,i,j))/delta)+(s_bottom...
            *(H_y(1,i,j+1) - H_y(3,i,j+1))/2))+H_y(3,i,j+1);        
        % Finite Difference Equation (1) from our notes
        E(1,i,j) = E(2,i,j) + (delt/e_bottom)*((1/delta)*...
            (H_z(L+1,i+1,j)-H_z(L+1,i-1,j) - H_y(L+1,i,j+1) + H_y(L+1,i,j-1)));                    
    case 'upper'
        
    otherwise
        
end
       
                    
            
            
            
            
            
            
            
        end
    end
    % Update the row vectors
    H_y(3,:,:) = H_y(2,:,:);
    H_y(2,:,:) = H_y(1,:,:);
    H_z(3,:,:) = H_z(2,:,:);
    H_z(2,:,:) = H_z(1,:,:);    
    E(3,:,:) = E(2,:,:);
    E(2,:,:) = E(1,:,:);
end


































