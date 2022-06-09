function [param] = check_param(param)
% Check if the parameters are assigned reasonable values. If no value is
% assigned by the user, the default value will be used for that parameter.
% Author: Rizwan Ahmad (ahmad.46@osu.edu)


% param.p
if ~isfield(param, 'PE') || isempty(param.PE)
    error('No value has been assigned to param.PE');
elseif param.PE < 2 || rem(param.PE,1) ~= 0
    error('The value assigned to param.PE must be an integer greater than 1');
end

% param.FR
if ~isfield(param, 'FR') || isempty(param.FR)
    error('No value has been assigned to param.FR');
elseif param.FR < 2 || rem(param.FR,1) ~= 0
    error('The value assigned to param.FR must be an integer greater than 1');
elseif param.FR > 32
    fprintf(['----------------------------------- Tip ------------------------------------' '\n' ...
             'For faster processing, reduce the number of frames (to lets say 32) and then' '\n' ...
             'cyclically reuse these frames to achieve any arbitrarilty number of frames.' '\n'...
             '----------------------------------------------------------------------------' '\n']);
end

% param.R
param.R = param.PE/param.n; % acceleration rate
if param.R < 1 || param.R > param.PE
    error('The value assigned to param.R must be an integer between 1 and param.PE');
end

% param.s
if ~isfield(param, 's') || isempty(param.s)
    error('No value has been assigned to param.s');
elseif param.s < 1 || param.s > 10
    error('The value assigned to param.s must be between 1 and 10');
end

% param.sig
if ~isfield(param, 'sig') || isempty(param.sig)
    error('No value has been assigned to param.sig');
elseif param.sig <= 0
    error('The value assigned to param.sig must be greater than zero');
end


%% Parameters with preset default values
% param.nIter
if ~isfield(param, 'nIter') || isempty(param.nIter)
    param.nIter = 120;
elseif param.nIter < 20 || param.nIter > 1024 || rem(param.nIter,1) ~= 0
    error('The value assigned to param.nIter must be an integer between 20 and 1024');
end

% param.ss
if ~isfield(param, 'ss') || isempty(param.ss)
    param.ss = 0.25; % default
elseif param.ss <= 0
    error('The value assigned to param.ss must be greater than zero');
end

% param.tf
if ~isfield(param, 'tf') || isempty(param.tf)
    param.tf = 0.0; % default
elseif param.tf < 0
    error('The value assigned to param.tf must be greater than or equql to zero');
elseif param.tf > 0
    warning('param.tf > 0 will destroy the constant temporal resolution.');
end

% param.beta
if ~isfield(param, 'beta') || isempty(param.beta)
    param.beta = 1.4; % default
elseif param.beta <= 0 || param.beta > 10
    error('The value assigned to param.beta must be between 0 and 10');
end

% param.g
if ~isfield(param, 'g') || isempty(param.g)
    param.g = floor(param.nIter/6); % default
elseif param.g < 5 || param.g > round(param.nIter/4) || rem(param.g,1) ~= 0
    error('The value assigned to param.g must be between 1 and param.nIter/4');
end

% param.uni
if ~isfield(param, 'uni') || isempty(param.uni)
    param.uni = round(param.nIter/2); % default
elseif param.uni < round(param.nIter/4) || param.uni > round(param.nIter/2) || rem(param.uni,1) ~= 0
    error('The value assigned to param.uni must be between param.nIter/4 and param.nIter/2');
end

% param.fl
if ~isfield(param, 'fl') || isempty(param.fl)
    param.fl= round(5/6*param.nIter); % default
elseif param.fl <= param.nIter/2 || param.fl > param.nIter || rem(param.fl,1) ~= 0
    error('The value assigned to param.fl must be between param.nIter/2 and param.nIter');
end

% param.W
if ~isfield(param, 'W') || isempty(param.W)
%     param.W = max(param.R/8, 1); % default
    param.W = param.R/10 + 0.25;
elseif param.W <= 0
    error('The value assigned to param.W must be greater than zero');
end

% param.sz
if ~isfield(param, 'sz') || isempty(param.sz)
    param.sz = 3.5; % default
elseif param.sz <= 0
    error('The value assigned to param.sz must be greater than zero');
end

% param.dsp
if ~isfield(param, 'dsp') || isempty(param.dsp)
    param.dsp = 5; % default
elseif param.dsp < 0 || param.dsp > param.nIter
    error('The value assigned to param.dsp must be between 0 and nIter');
end

% param.fs
if ~isfield(param, 'fs') || isempty(param.fs)
    param.fs = 1; % default
elseif (param.fs ~= 0 && param.fs ~= 1 && param.fs < param.R) || rem(param.fs,1)
    error('The value assigned to param.fs must be a non-negative integer');
end

% param.fc
if ~isfield(param, 'fc') || isempty(param.fc)
    param.fc = 1; % default
elseif (param.fc <= 0 || param.fc > 1)
    error('The value assigned to param.fc must be zero or one');
end
