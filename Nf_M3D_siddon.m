function [leng] = Nf_M3D_siddon(param, pt1, pt2, phan)
% Nf_M3D function is used to calculate the inetersection length of a single
% ray and the phantom volume
% param is the struct array of system parameters
% pt1 is the coordinate of x-ray source point
% pt2 is the coordinate of detector unit
% phan1 is the phantom of uranium
% phan2 is the phantom of iron
% phan3 is the phantom of LiH

leng = 0.0;
% lambda = 1e-12;

Nx = param.nVoxelX + 1;% the number of planes in x direction
Ny = param.nVoxelY + 1;% the number of planes in y direction
Nz = param.nVoxelZ + 1;% the number of planes in z direction

X_planes = (-param.nVoxelX/2:param.nVoxelX/2) * param.fVoxel;
Y_planes = (-param.nVoxelY/2:param.nVoxelY/2) * param.fVoxel;
Z_planes = (-param.nVoxelZ/2:param.nVoxelZ/2) * param.fVoxel;

ray.x = pt2.x - pt1.x;% the x coordinate of the ray vector
ray.y = pt2.y - pt1.y;% the y coordinate of the ray vector
ray.z = pt2.z - pt1.z;% the z coordinate of the ray vector

axm = min((X_planes(1) - pt1.x) / ray.x, (X_planes(Nx) - pt1.x) / ray.x);
aym = min((Y_planes(1) - pt1.y) / ray.y, (Y_planes(Ny) - pt1.y) / ray.y);
azm = min((Z_planes(1) - pt1.z) / ray.z, (Z_planes(Nz) - pt1.z) / ray.z);

axM = max((X_planes(1) - pt1.x) / ray.x, (X_planes(Nx) - pt1.x) / ray.x);
ayM = max((Y_planes(1) - pt1.y) / ray.y, (Y_planes(Ny) - pt1.y) / ray.y);
azM = max((Z_planes(1) - pt1.z) / ray.z, (Z_planes(Nz) - pt1.z) / ray.z);

am = max(max(axm, aym), azm);
aM = min(min(axM, ayM), azM);

ax = (X_planes - pt1.x) / (ray.x );
ay = (Y_planes - pt1.y) / (ray.y );
az = (Z_planes - pt1.z) / (ray.z );

if am>=aM
    return;% the ray does't intersect with the volume
end

% for x
if pt1.x <= pt2.x
    i_min = ceil(Nx - (X_planes(end) - am * (pt2.x - pt1.x) - pt1.x) / param.fVoxel);
    i_max = floor(1 + (pt1.x + aM * (pt2.x - pt1.x) - X_planes(1)) / param.fVoxel);
    ax = ax(i_min:i_max);
else
    i_min = ceil(Nx - (X_planes(end) - aM * (pt2.x - pt1.x) - pt1.x) / param.fVoxel);
    i_max = floor(1 + (pt1.x + am * (pt2.x - pt1.x) - X_planes(1)) / param.fVoxel);
    ax = ax(i_max:-1:i_min);
end

% for y
if pt1.y <= pt2.y
    j_min = ceil(Ny - (Y_planes(end) - am * (pt2.y - pt1.y) - pt1.y) / param.fVoxel);
    j_max = floor(1 + (pt1.y + aM * (pt2.y - pt1.y) - Y_planes(1)) / param.fVoxel);
    ay = ay(j_min:j_max);
else
    j_min = ceil(Ny - (Y_planes(end) - aM * (pt2.y - pt1.y) - pt1.y) / param.fVoxel);
    j_max = floor(1 + (pt1.y + am * (pt2.y - pt1.y) - Y_planes(1)) / param.fVoxel);
    ay = ay(j_max:-1:j_min);
end

% for z
if pt1.z <= pt2.z
    k_min = ceil(Nz - (Z_planes(end) - am * (pt2.z - pt1.z) - pt1.z) / param.fVoxel);
    k_max = floor(1 + (pt1.z + aM * (pt2.z - pt1.z) - Z_planes(1)) / param.fVoxel);
    az = az(k_min:k_max);
else
    k_min = ceil(Nz - (Z_planes(end) - aM * (pt2.z - pt1.z) - pt1.z) / param.fVoxel);
    k_max = floor(1 + (pt1.z + am * (pt2.z - pt1.z) - Z_planes(1)) / param.fVoxel);
    az = az(k_max:-1:k_min);
end

a = unique(sort([am, ax, ay, az, aM]));
a = round(a*1.0e8)/1.0e8;

total_len = sqrt(ray.x * ray.x + ray.y * ray.y + ray.z * ray.z); % total length of the ray
len = zeros(length(a)-1, 1);

% the intersection length
for i=1:length(len)
    len(i) = total_len * (a(i+1) - a(i));
end

index = zeros(length(len), 3);

for i=1:length(len)
    a_mid = (a(i+1) + a(i)) / 2;
    xx = (pt1.x + a_mid * (pt2.x - pt1.x) - X_planes(1)) / param.fVoxel;
    yy = (pt1.y + a_mid * (pt2.y - pt1.y) - Y_planes(1)) / param.fVoxel;
    zz = (pt1.z + a_mid * (pt2.z - pt1.z) - Z_planes(1)) / param.fVoxel;
    if abs(xx)<=1e-6
        xx = 0;
    end
    if abs(yy)<=1e-6
        yy = 0;
    end
    if abs(zz)<=1e-6
        zz = 0;
    end
    index(i,1) = floor(xx+1);
    index(i,2) = floor(yy+1);
    index(i,3) = floor(zz+1);
end

for ii=1:length(len)
    if index(ii,1)>=1 && index(ii,1)<=param.nVoxelY && ...
            index(ii,2)>=1 && index(ii,2)<=param.nVoxelX &&...
            index(ii,3)>=1 && index(ii,3)<=param.nVoxelZ
        leng = leng + len(ii) * phan(index(ii,2), index(ii,1), index(ii,3));
    end
end

end

