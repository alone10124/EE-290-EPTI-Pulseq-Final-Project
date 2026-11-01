function lowresI1= recon_lowresI2(kdata_rep1,ky_acq,sens_map_all,ceter_index,ceter_size, parameters,linearphase_PF)
% reconstruct the low resolution images with the densely sampled data
warning("off")
sms = parameters.SMS;
ny = parameters.Rseg;
[necho,nro,~,ncc,nslice] = size(kdata_rep1);
sms_phase = mod([1:length(ky_acq)]+1,2)+1;
Matrix_encode = Generate_SMS_encode_Matrix2(nro,ny,sms,sms);
Matrix_encode = Matrix_encode(:,:,:,1,1,1,end:-1:1);
lowresI1 = zeros(nro,ny,nslice*sms,2);         % different echos

% reconstruct the low reslution images
for dif_slice = 1:nslice
    slice_group = dif_slice:nslice:nslice*sms;
    slice_group = slice_group(end:-1:1);
    smap = sens_map_all(:,:,slice_group,:);
    %combine the data1
    ind = [ceter_index(1):ceter_index(1)+ceter_size,ceter_index(2):ceter_index(2)+ceter_size];
    odd_ind = ind(find(sms_phase(ind)==1));
    even_ind = ind(find(sms_phase(ind)==2));
    combineK1_1 = zeros(ny,nro,1,ncc,1,1,sms);
    combineK1_1(ky_acq(even_ind),:,1,:,1,1,1,1) = kdata_rep1(even_ind,:,1,:,dif_slice);
    combineK1_1(ky_acq(odd_ind),:,1,:,1,1,2,1) = kdata_rep1(odd_ind,:,1,:,dif_slice);
    combineK1_1 = permute(combineK1_1,[2,1,3,4,5,6,7,8]);

    %combine the data2
    ind = [ceter_index(3):ceter_index(3)+ceter_size,ceter_index(4):ceter_index(4)+ceter_size];
    odd_ind = ind(find(sms_phase(ind)==1));
    even_ind = ind(find(sms_phase(ind)==2));
    combineK2_1 = zeros(ny,nro,1,ncc,1,1,sms);
    combineK2_1(ky_acq(even_ind),:,1,:,1,1,1,1) = kdata_rep1(even_ind,:,1,:,dif_slice);
    combineK2_1(ky_acq(odd_ind),:,1,:,1,1,2,1) = kdata_rep1(odd_ind,:,1,:,dif_slice);
    combineK2_1 = permute(combineK2_1,[2,1,3,4,5,6,7,8]);
    %recon
    im_combo1=sum(sum(ifft2c(combineK1_1).*conj(smap),4) .*conj(Matrix_encode(:,:,:,1,1,1,:)),7).*linearphase_PF;
    im_combo2=sum(sum(ifft2c(combineK2_1).*conj(smap),4) .*conj(Matrix_encode(:,:,:,1,1,1,:)),7).*linearphase_PF;

    lowresI1(:,:,slice_group,1) = squeeze(im_combo1);
    lowresI1(:,:,slice_group,2) = squeeze(im_combo2);
end

% reconstruct with iteration