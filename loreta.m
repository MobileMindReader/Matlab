% Loreta implementation

% Data loading
load('../EEG data/16 EEG/A01G.mat');
A = load('model data/forwardmodel_spheres_reduced.txt');
L = load('model data/spatialCoherenceSmooth0-2_reduced.txt');
LInv = load('model data/spatialCoherenceSmooth0-2_reduced_inverse.txt');

Y = data{1,1}.X;
channels = data{1,1}.channels; % Mismatch between constants and test data

% Setup
sourceDim = 1028;
sources = zeros(length(Y)-2, sourceDim);
% sourceCovar = zeros(sourceDim, sourceDim);
% sensorCovar = zeros(length(channels)-2, length(channels)-2);

% Initial parameter values
alphaInv = 0.0100;
betaInv = 0.3781;



%%

% TODO: Subtract mean from signal 

sensorCovarInv = alphaInv * (A*L')*(L*A') + betaInv*eye(length(channels)-2);
sensorCovar = inv(sensorCovarInv);

for sample = 1 : length(Y)-2
    
    sources(sample,:) = alphaInv * (A' * sensorCovar) * Y(sample,1:14)';
    
end


%% Parameter updating with EM-algorithm

sourceCovar = alphaInv*eye(sourceDim) - alphaInv * A' * sensorCovar * A * alphaInv;

tempAlphaInv = sum(sum(LInv * sourceCovar'));