function [len] = Ax( param, angleIndex, phan )
% Ax: this function is used to generate a single projection image
% param: the system setting
% Dense_uranium: the data for uranium, the first column is the energy bin, 
% the second column is the mass attenuation coefficients and the third
% column is the mass energy-absorption coefficients
% uranium_pos: the positions where the volume contains the relevant
% material

% the param.nVoxelX, param.nVoxelY, param.nVoxelZ should be equal to the volume size
assert(isequal(param.nVoxelX, length(phan(1,1,:))), 'bad geometry setting', 'the param should be equal to the volume');

% the input volume should three-dimensional
assert(size(phan, 3)>1, 'Invalid Input', 'image should be 3D');

% sin_angle = sind(angleIndex);
% cos_angle = cosd(angleIndex);

% the return value is three-dimensional consisting of the sum of the
% intersecting length for the three materials
% proj = Ax_mex(param, angleIndex, fTable, uranium_pos, iron_pos, LiH_pos);

len = zeros(param.nDetV, param.nDetU, 'single');

pt0.x = param.fDso;
pt0.y = 0;
pt0.z = 0;

% pt_det = zeros(param.nDetV, param.nDetU, 'single');
x_temp = -(param.fDsd - param.fDso);
y_temp = (-(param.nDetU-1)/2:1:(param.nDetU-1)/2) * param.fDetUnit;
z_temp = (-(param.nDetV-1)/2:1:(param.nDetV-1)/2) * param.fDetUnit;

for ii=1:param.nDetV
    for jj=1:param.nDetU
       pt1.x = pt0.x * cosd(angleIndex) - pt0.y * sind(angleIndex);
       pt1.y = pt0.x * sind(angleIndex) + pt0.y * cosd(angleIndex);
       pt1.z = 0;
       
       pt2.x = x_temp * cosd(angleIndex) - y_temp(jj) * sind(angleIndex);
       pt2.y = x_temp * sind(angleIndex) + y_temp(jj) * cosd(angleIndex);
       pt2.z = z_temp(ii);
       
       len(jj, ii) = Nf_M3D(param, pt1, pt2, phan);
       
    end
end