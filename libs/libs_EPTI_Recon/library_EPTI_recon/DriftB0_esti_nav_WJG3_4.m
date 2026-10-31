function b0_delta = DriftB0_esti_nav_WJG3_4(lowres_refI,lowresI_rep,dt,ceter_index,msk2,motion_param_file,BART_savepath)
% Estimate the  B0 offset
% msk1 is the mask for lowres_refI
% msk2 is the mask for lowresI_rep,by applying the motion parameters to the
% msk
warning("off")
[nro,npe,nslice,~] = size(lowres_refI);
b0_delta = zeros(nro,npe,nslice);
time_gap = (ceter_index(3)+ceter_index(4))/2-(ceter_index(1)+ceter_index(2))/2;
rmse = zeros(1,nslice);
unstable_slice = [];
% fitting
navI11 = lowres_refI(:,:,:,1);
navI21 = lowresI_rep(:,:,:,1);
if (~isempty(motion_param_file))
    navI21_r = warp_wjg(real(navI21),motion_param_file,BART_savepath);
    navI21_i = warp_wjg(imag(navI21),motion_param_file,BART_savepath);
    navI21 = navI21_r+1i*navI21_i;
end
navI12 = lowres_refI(:,:,:,2);
navI22 = lowresI_rep(:,:,:,2);
if (~isempty(motion_param_file))
    navI22_r = warp_wjg(real(navI22),motion_param_file,BART_savepath);
    navI22_i = warp_wjg(imag(navI22),motion_param_file,BART_savepath);
    navI22 = navI22_r+1i*navI22_i;
end


temp1 = navI21.*conj( navI11 ) ./( abs(navI21).*abs(navI11) + eps);
temp2 = navI22.*conj( navI12 ) ./( abs(navI22).*abs(navI12) + eps);
PH1 = temp2.*conj(temp1);
b0_avg = angle(PH1)/time_gap/dt/2/pi;

for loopi = 1:nslice
    msk = (abs(msk2(:,:,loopi))>0.9);
    ind = find(msk>0);
    b0_temp = b0_avg(:,:,loopi);
    %remove the outlier
    outlier_idx = isoutlier(b0_temp(ind));
    outlier_idx = find(outlier_idx==0);
    ind = ind(outlier_idx);
    %fitting
    if (~isempty(ind))&&(length(ind)>1/40*nro*npe)
        [row, col] = ind2sub([nro,npe], ind);
        z = b0_temp(ind);
        [sf, gof] = fit([row, col],z,'poly22');
        [X,Y] = meshgrid(1:npe,1:nro);
        b0_avg2 = sf([Y(:),X(:)]);
        b0_avg2 = reshape(b0_avg2,[nro,npe]);
        % subplot(122);imshow(b0_avg2.*msk(:,:,slice_group(select_slice)),[-2,2]);colormap jet;
        b0_delta(:,:,loopi) = b0_avg2;
        rmse(loopi) = gof.rmse;
    else
        b0_delta(:,:,loopi) = zeros(nro,npe);
        unstable_slice = [unstable_slice,loopi];
    end
end

% if the fitting is unstable then using the global value
global_B0 = sum(b0_delta(:,:,setdiff(1:nslice,unstable_slice)),[1,2,3]);
totle_pixel = sum(b0_delta(:,:,setdiff(1:nslice,unstable_slice))~=0,[1,2,3]);
global_B0 = global_B0/(totle_pixel+eps);
for loopi = unstable_slice
    b0_delta(:,:,loopi)=global_B0;
end

% % imshow3((b0_delta).*msk2,[-5,5]);colormap jet
% % imshow3(b0_delta.*msk,[-5,5]);colormap jet
