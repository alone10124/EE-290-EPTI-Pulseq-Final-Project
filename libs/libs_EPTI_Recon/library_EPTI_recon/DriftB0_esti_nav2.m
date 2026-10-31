
function b0_avg = DriftB0_esti_nav2(nav1,nav2,dt)
% Estimate the global B0 offset / field drift using 1D phasecor navigators

nav1 = ifftc(nav1,2);
nav2 = ifftc(nav2,2);

temp = nav2.*conj( nav1 ) ./( abs(nav2).*abs(nav1) + eps);
temp = temp(3,:,:,:,:,:).*conj(temp(1,:,:,:,:,:)).* abs(nav1(1,:,:,:,:,:));
temp = temp(:,:,:,1);
[cleanData,out_ind] = rmoutliers(angle(temp));
msk = sum(abs(temp),4);
ind = msk>(max(msk)/2);
ind = ind.*~out_ind;
[npe,nx,nz,~] = size(temp);
msk = reshape(ind,npe,nx,nz);
temp = sum(temp,4);
temp = sum(temp(:).*msk(:))./sum(msk(:));
b0_avg = angle(temp)/2/dt/2/pi;

