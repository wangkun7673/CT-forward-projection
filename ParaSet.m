% the fan-beam or cone-beam system parameters initialization
% scanning mode, the default value is cone-beam, the other is fan-beam
param.bSysType = 1; % 1 means that the system is cone-beam CT system and 0 means fan-beam CT system

% the distance from x-ray source to the object in cm
param.fDso = 115;
% the distance from x-ray source to the detector in cm
param.fDsd = 190;

% the number of rotation angle
param.nAngle = 360;

% the resolution of the detector
param.nDetU = 256; % u direction
param.nDetV = 256; % v direction
param.fDetUnit = 0.01;

% the resolution of the reconstructed volume
param.nVoxelX = 256;
param.nVoxelY = 256;
param.nVoxelZ = 256;
param.fVoxel = param.fDetUnit * param.fDso / param.fDsd;

% weather to use Nvidia GPU
param.bGPU = 1;% the default value is 1, which means that gpu is supported
