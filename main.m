clear;
close all;
clc; 

%% ================
% Orignal author: Rizwan Ahmad (ahmad.46@osu.edu)
% Updated by: Mihir Joshi (joshi.373@osu.edu)

%% ================
% Add folders to path
restoredefaultpath
addpath(genpath('./functions/'));

%% ================
% Modify Flag to Select a Sampling Method
fg = "pr4d";
% Options: "vista, "gro", "cava", "opra", "pr4d"

%% ================
% VISTA Parameters 
param_vista.PE    = 160; % Size of of phase encoding (PE) grid
param_vista.FR    = 64; % Number of frames
param_vista.n     = 12;  % Number of samples (readouts) per frame
param_vista.M     = param_vista.FR*param_vista.n; % Total number of samples
param_vista.s     = 1.6; % 1<=s<=10 controls sampling density; 1: uniform density, 10: maximally non-uniform density
param_vista.sig   = param_vista.PE/6; % Std of the Gaussian envelope for sampling density
param_vista.w     = max((param_vista.PE/param_vista.n)/10 + 0.25, 1); % Scaling of time dimension; frames are "w" units apart
param_vista.beta = 1.4; % Exponent of the potenital energy term.
% Probably you don't need to change the following VISTA paramters
% If unsure, leave them empty and the default value will be employed.
param_vista.sd   = 10; % Seed to generate random numbers; a fixed seed should reproduce VISTA
param_vista.nIter= []; % Number of iterations for VISTA (defualt: 120)
param_vista.ss   = []; % Step-size for gradient descent. Default value: 0.25; 
param_vista.tf   = []; % Step-size in time direction wrt to phase-encoding direction; use zero for constant temporal resolution. Default value: 0.0
param_vista.g    = []; % Every gth iteration is relocated on a Cartesian grid. Default value: floor(param.nIter/6)
param_vista.uni  = []; % At param.uni iteration, reset to equivalent uniform sampling. Default value: floor(param.nIter/2)
param_vista.sz   = []; % Display size of samples. Default value: 3.5
param_vista.dsp  = 5;  % Display frequency (verbosity), every dsp-th iteration will be displayed. Default value: 5
param_vista.fs   = [];  % Does time average has to be fully sampled, 0 for no, 1 for yes. Default value: 1
param_vista.fc   = []; % What fraction of the time-averaged should be fully sampled. Default value: 1
param_vista.fl   = []; % Start checking fully sampledness of time-averge at fl^th iteration. Default value: floor(param.nIter*5/6)

%% ================
% GRO Parameters
param_gro.PE   = 160;  % Size of of phase encoding (PE) grid
param_gro.FR   = 64;   % Numbe of frames
param_gro.n    = 12;   % Number of samples (readouts) per frame
param_gro.M    = param_gro.FR*param_gro.n; % Total number of samples
param_gro.E    = 1;    % Number of encoding, E=1 for cine, E=2 for flow (phase-contrast MRI)
param_gro.tau  = 1;    % Extent of shift between frames, tau = 1 or 2: golden ratio shift, tau>2: tiny golden ratio shift
param_gro.s    = 2.2;  % s>=1. larger values means higher sampling density in the middle (default: 2.2)
param_gro.alph = 3;    % alph>1. larger alpha means sharper transition from high-density to low-density regions (default: 3)
param_gro.PF   = 0;    % for partial fourier; discards PF samples from one side (default: 0)
param_gro.dsp  = 1;    % Display figures: 0 no, 1 yes

%% ================
% CAVA Parameters
param_cava.PE   = 120;  % Size of of phase encoding (PE) grid
param_cava.FR   = 48;  % Nominal number of frames (for display only)
param_cava.n    = 6;   % Nominal number of samples per frame per encoding (for display one)
param_cava.M    = param_cava.FR*param_cava.n; % Total number of samples
param_cava.E    = 2;   % Number of encoding, E=1 for cine, E=2 for flow (phase-contrast MRI)
param_cava.tau  = 1;   % Size of jumps, tau = 1 or 2: golden ratio jump, tau>2: tiny golden ratio jump
param_cava.s    = 2.2;  % s>=1. larger values means higher sampling density in the middle (default: 2.2)
param_cava.alph = 3;  % alph>1. larger alpha means sharper transition from high-density to low-density regions (default: 3)
param_cava.dsp  = 1;   % Display figures: 0 no, 1 yes

%% ================
% OPRA Parameters
param_opra.PE   = [96, 60]; % Size of of phase encoding ky-kz grid (larger number first)
param_opra.FR   = 80;  % Number of frames (for display only)
param_opra.n    = 30;  % Nominal number of samples per frame (for display one)
param_opra.L    = 10;  % Number of samples in each "L" leaflet (keep it even)
param_opra.M    = param_opra.FR*param_opra.n; % Total number of samples
param_opra.phi  = pi/12;  % Angular jump (in radian) from the end of one leaflet to the start of the other
param_opra.s    = 3; % s>=0; s=1 uniform density, s>1 means more samples in the middle
param_opra.ar   = 0; % controls the aspect ratio/orientation of the high density ellipsoid region in the middle
param_opra.gs   = 6^(1/2); % a second irrational number for radial advancement
param_opra.cg   = 1/3; % cg>=0, less than 1/2 for no gap in the middle, more than 1/2 for a hole at kx=ky=0 of size round(cg). 
                      % This hole can accomodate self-gating samples at L/2, 3L/2, 5L/2...
param_opra.dsp  = 1;   % Display figures: 0 no, 1 yes

%% ================
% PR4D Parameters
param_pr4d.PE   = [96, 60]; % Size of the phase encoding ky-kz grid (larger number first)
param_pr4d.FR   = 80; % Nominal number of frames (for display only)
param_pr4d.n    = 30; % Nominal number of samples per frame per encoding (for display one)
param_pr4d.M    = param_pr4d.FR*param_pr4d.n; % Total number of samples
param_pr4d.s    = 3; % s>=0; s=1 uniform density, s>1 means more samples in the middle
param_pr4d.ar   = 0; % controls the aspect ratio/orientation of the high density ellipsoid region in the middle
param_pr4d.gs   = 35^(1/3); % a second irrational number for angular advancement
param_pr4d.E    = 4; % Number of encodings E>=1
param_pr4d.cg   = 1/3; % cg>=0, less than 1/2 for no gap in the middle, more than 1/2 for a hole at kx=ky=0 of size round(cg)
                     % This hole can accomodate self-gating samples (not included)
param_pr4d.dsp  = 1;   % Display figures: 0 no, 1 yes


%% ================
% Run Sampling Method Based on Flag
if fg == "vista" % VISTA
    param_vista = check_param(param_vista); % Check parameters
    [PEInd, FRInd, samp] = vista_fun(param_vista); 
    % ky(i) = PEInd(i), 'i' is the order of acquisition
    % t(i) = FRInd(i) 
    % samp = binary mask in ky-t domain
elseif fg == "gro" % GRO
    [PEInd, FRInd, samp] = gro_fun(param_gro);   
    % ky(i,e) = PEInd(i,e), where 'e' is encoding and 'i' is the order of acquisition
    % t(i) = FRInd(i)
    % samp = binary mask in ky-t-encoding domain
elseif fg == "cava" % CAVA
    [PEInd, FRInd, samp] = cava_fun(param_cava); 
    % ky(i,e) = PEInd(i,e), where 'e' is encoding and 'i' is the order of acquisition
    % t(i) = FRInd(i)
    % samp = binary mask in ky-t-encoding domain
elseif fg == "opra" % OPRA
    [PEInd] = opra_fun(param_opra); 
    % ky(i) = PEInd(i,1); kz(i) = PEInd(i,2) where 'i' is the order of acquisition
elseif fg == "pr4d" % PR4D
    [PEInd] = pr4d_fun(param_pr4d); 
    % ky(i,e) = PEInd(i,1,e); kz(i,e) = PEInd(i,2,e) where 'e' is encoding and
    % 'i' is the order of acquisition
end