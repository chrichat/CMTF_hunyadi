% Demo of coupled polynomial code
% -------------------------------
%
% Run the code to get EEG and fMRI first.

% Author(s): Nico Vervliet
%
% Version History:
% - 2015/11/11   NV      Initial version


%% Set some parameters
E = permute(EEG, [2 3 1]); % set first mode to last (not necessary)
F = fMRI;
R = 3;
d = 3;

%% A two stage initialization is computed
% 1. Compute CPD of E
Uinit = cpd(E, R);
Uinit{1} = Uinit{1}*(-1); 
Uinit{3} = Uinit{3}*(-1); % I like to work with positive coordinates
                          
% 2. Compute restricted CPD of F
% struct_poly2 uses different points for each column
clear model sol
model.variables.c = randn(R,d);
model.variables.e = randn(size(F,2),R);
model.factors.D = {'c', @(z,t) struct_poly2(z,t,Uinit{3})}; 
model.factors.E = {'e'};
model.factorizations.fmri.data = F;
model.factorizations.fmri.cpd = {'D', 'E'};

options.MaxIter = 1000;
[sol, output] = sdf_nls(model,options);

%% Compute the solution for the coupled problem, using the good init.
clear model;
model.variables.a = Uinit{1};
model.variables.b = Uinit{2};
model.variables.ct = {sol.variables.c, Uinit{3}}; % merge coefs and coords
model.variables.e = sol.variables.e;

model.factors.A = 1;
model.factors.B = 2;
model.factors.C = {3, @(z,t) struct_select(z,t,2)}; % extract coords
model.factors.D = {3, @struct_bori}; % add poly constraint
model.factors.E = 4;

model.factorizations.eeg.data = E;
model.factorizations.eeg.cpd = {1,2,3};
model.factorizations.fmri.data = F;
model.factorizations.fmri.cpd = {4,5};

clear options;
options.MaxIter = 5000;
options.Display = 100;
options.TolFun = 1e-10;
options.TolX = 1e-8;
options.RelWeights = [1 1];

[sol, output] = sdf_nls(model,options);

%% Check result
% Compute the errors for best uncoupled solutions
frob(cpdres(E,Uinit))/frob(E)
frob(cpdres(F,cpd(F,R)))/frob(F)

% Compute the errors for the coupled solution
frob(cpdres(E, {sol.factors.A, sol.factors.B, sol.factors.C}))/frob(E)
frob(cpdres(F, {sol.factors.D, sol.factors.E}))/frob(F)

% Plot factor matrix D in function of the coords to see if it is a polynomial
plot(sol.factors.C, sol.factors.D,'.')
