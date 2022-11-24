% this program is used to test jacob_ray_projection.cu file
clear
clc
% system parameters setting
ParaSet;
% phantom generation
phan = single(phantom3dAniso('modified shepp-logan', param.nVoxelX));
% figure, imshow(phan(:,:,128), [])

proj = Ax(param, 0, phan);
proj_siddon = Ax_siddon(param, 0, phan);
err = proj - proj_siddon;

figure, hold on 
subplot(131), imshow(proj, []), title('Jacobs algorithm')
subplot(132), imshow(proj_siddon, []), title('Siddon algorithm')
subplot(133), imshow(err, []), colorbar, title('error between Jacobs and Siddon')