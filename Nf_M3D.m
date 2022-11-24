function [len] = Nf_M3D(param, pt1, pt2, phan)
% Nf_M3D function is used to calculate the inetersection length of a single
% ray and the phantom volume
% param is the struct array of system parameters
% pt1 is the coordinate of x-ray source point
% pt2 is the coordinate of detector unit
% phan is the phantom

len = 0.0;

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

if am>=aM
    return;% the ray does't intersect with the volume
end

% for x
if pt1.x < pt2.x
    if am == axm
        i_min = 2;% the first plane that the ray intersects with after the ray enters the voxel space
    else
        i_min = ceil((pt1.x + am * ray.x - X_planes(1)) / param.fVoxel);
    end
    
    if aM == axM
        i_max = Nx;
    else
        i_max = floor((pt1.x + aM * ray.x - X_planes(1)) / param.fVoxel);
    end
else
    if am == axm
        i_max = Nx - 1;
    else
        i_max = floor((pt1.x + am * ray.x - X_planes(1)) / param.fVoxel);
    end
    
    if aM == axM
        i_min = 1;
    else
        i_min = ceil((pt1.x + aM * ray.x - X_planes(1)) / param.fVoxel);
    end
end

% for y
if pt1.y < pt2.y
    if am == aym
        j_min = 2;% the first plane that the ray intersects with after the ray enters the voxel space
    else
        j_min = ceil((pt1.y + am * ray.y - Y_planes(1)) / param.fVoxel);
    end
    
    if aM == ayM
        j_max = Ny;
    else
        j_max = floor((pt1.y + aM * ray.y - Y_planes(1)) / param.fVoxel);
    end
else
    if am == aym
        j_max = Ny - 1;
    else
        j_max = floor((pt1.y + am * ray.y - Y_planes(1)) / param.fVoxel);
    end
    
    if aM == ayM
        j_min = 1;
    else
        j_min = ceil((pt1.y + aM * ray.y - Y_planes(1)) / param.fVoxel);
    end
end

% for z
if pt1.z < pt2.z
    if am == azm
        k_min = 2;% the first plane that the ray intersects with after the ray enters the voxel space
    else
        k_min = ceil((pt1.z + am * ray.z - Z_planes(1)) / param.fVoxel);
    end
    
    if aM == azM
        k_max = Nz;
    else
        k_max = floor((pt1.z + aM * ray.z - Z_planes(1)) / param.fVoxel);
    end
else
    if am == azm
        k_max = Nz - 1;
    else
        k_max = floor((pt1.z + am * ray.z - Z_planes(1)) / param.fVoxel);
    end
    
    if aM == azM
        k_min = 1;
    else
        k_min = ceil((pt1.z + aM * ray.z - Z_planes(1)) / param.fVoxel);
    end
end

lambda = 1e-12;
if pt1.x < pt2.x
    ax = (i_min * param.fVoxel + X_planes(1) - pt1.x) / (ray.x + lambda);
else
    ax = (i_max * param.fVoxel + X_planes(1) - pt1.x) / (ray.x + lambda);
end
if pt1.y < pt2.y
    ay = (j_min * param.fVoxel + Y_planes(1) - pt1.y) / (ray.y + lambda);
else
    ay = (j_max * param.fVoxel + Y_planes(1) - pt1.y) / (ray.y + lambda);
end
if pt1.z < pt2.z
    az = (k_min * param.fVoxel + Z_planes(1) - pt1.z) / (ray.z + lambda);
else
    az = (k_max * param.fVoxel + Z_planes(1) - pt1.z) / (ray.z + lambda);
end

aminc = min(min(ax, ay), az);
i = int32(floor((pt1.x + (aminc + am)/2 * ray.x - X_planes(1)) / param.fVoxel));
j = int32(floor((pt1.y + (aminc + am)/2 * ray.y - Y_planes(1)) / param.fVoxel));
k = int32(floor((pt1.z + (aminc + am)/2 * ray.z - Z_planes(1)) / param.fVoxel));

ac = am;
axu = param.fVoxel / abs(ray.x);
ayu = param.fVoxel / abs(ray.y);
azu = param.fVoxel / abs(ray.z);

if pt1.x < pt2.x
    iu = 1;
else
    iu = -1;
end

if pt1.y < pt2.y
    ju = 1;
else
    ju = -1;
end

if pt1.z < pt2.z
    ku = 1;
else
    ku = -1;
end

maxlength = sqrt(ray.x * ray.x + ray.y * ray.y + ray.z * ray.z);
sum_1 = 0.0;

Np = (i_max - i_min + 1) + (j_max - j_min + 1) + (k_max - k_min + 1);
for ii=0:Np
    if ax == aminc
        if i>=0 && i<param.nVoxelX && j>=0 && j<param.nVoxelY && k>=0 && k<param.nVoxelZ
            sum_1 = sum_1 + (ax - ac) * phan(j+1, i+1, k+1);
        end
        i = i + iu;
        ac = ax;
        ax = ax + axu;
    elseif ay == aminc
        if i>=0 && i<param.nVoxelX && j>=0 && j<param.nVoxelY && k>=0 && k<param.nVoxelZ
            sum_1 = sum_1 + (ay - ac) * phan(j+1, i+1, k+1);
        end
        j = j + ju;
        ac = ay;
        ay = ay + ayu;
    elseif az == aminc
        if i>=0 && i<param.nVoxelX && j>=0 && j<param.nVoxelY && k>=0 && k<param.nVoxelZ
            sum_1 = sum_1 + (az - ac) * phan(j+1, i+1, k+1);
        end
        k = k + ku;
        ac = az;
        az = az + azu;
    end
    aminc = min(min(ax, ay), az);
end

len = sum_1 * maxlength;% the sum of intersection length with the uranium phantom

end

