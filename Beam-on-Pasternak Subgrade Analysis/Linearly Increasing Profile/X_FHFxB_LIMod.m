% Script to determine the calibration factor chi for free-head fixed-base
% pile embeded in soil with linearly increasing modulus and eta = 0. 
%
% The script solves the BVP of beam in Kerr-equivalent Pasternak model with
% an inital assumption of chi (X_init). The pile head displacment is
% compared to FE analysis results with the same parameters. The percentage
% difference in the pile head displacment is then compared against a
% tolerance value (solDiff_Tol). If the the value is exceeded, chi is
% incremented (by X_inc) and the anlysis redone. The analysis is iterated
% untill the diviation is lower than the tolerance value. The final chi
% value is saved in the X_out vector. 

clc;
clear;
close all;

% Control Parameters
X_init = 0.01;
X_inc = 0.001;
solDiff_Tol = 1.0;
opts = bvpset('RelTol', 1.0e-6);

% Read analysis input parameters
pData_file = '../Parameters_LIModulus_FHFxB_eta_zero.csv';
pData = readmatrix(pData_file);
baseName = '../Abaqus_Output/LIModulus_FHFxB_eta_zero_Model_'; % FE output file

% Initalize the output
output = zeros(size(pData, 1), 8);

% Conduct the anlysis for each entry in the analysis data file
for j = 1:size(pData, 1)
    % Get parameters
    r = pData(j, 1);
    Lr = pData(j, 2);
    eta = pData(j, 3);
    nu_s = pData(j, 4);
    Ep_mr = pData(j, 5);
    nu_p = pData(j, 6);
    Ep = pData(j, 7);
    H = pData(j, 8);
    M = -pData(j, 9);

    % Read the FEM output file
    solFEMfileName = baseName + string(j-1); % the index is j-1, the abaqus automation script uses python and python uses 0 as the first index
    solFEM = readmatrix(solFEMfileName);

    % Inital guess for BVP and X
    L = r*Lr;
    zmesh = 0:0.05:L;
    uinit = [solFEM(1,2); 0; 0; 0];
    init_guess = bvpinit(zmesh, uinit);
    X = X_init;

    % Conduct initial analysis to establish solDiff
    p = [r Lr eta nu_s Ep_mr nu_p Ep H M X];
    solBVP = bvp5c(@(z, u)Lin_govDiffEq(z, u, p), @(u0, uL) Lin_FHFxB_BC(u0, uL, p), init_guess, opts);
    solDiff = abs(solFEM(1, 2) - solBVP.y(1,1))/abs(solFEM(1, 2))*100; %diviation in pile head displacment

    % Iterate untill solDiff < = solDiff_Tol
    while solDiff > solDiff_Tol
        X = X + X_inc;
        p = [r Lr eta nu_s Ep_mr nu_p Ep H M X];
        solBVP = bvp5c(@(z, u) Lin_govDiffEq(z, u, p), @(u0, uL) Lin_FHFxB_BC(u0, uL, p), init_guess, opts);
        solDiff = abs((solFEM(1, 2) - solBVP.y(1,1))/solFEM(1,2))*100;
        disp(X);
        disp(solDiff);
        if X > 5
            solDiff = 0;
            continue
        end
    end

    % output
    output(j, 1) = r;
    output(j, 2) = Lr;
    output(j, 3) = nu_s;
    output(j, 4) = Ep_mr;
    output(j, 5) = nu_p;
    output(j, 6) = Ep;
    output(j, 7) = X;
    output(j, 8) = solDiff;
end

% Write output to an Excel file
writematrix(output, 'X_LIMod_FHFxB.xlsx');





