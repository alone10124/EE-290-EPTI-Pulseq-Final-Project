function linear_fit_coeff = EPTI_SMS_Estimate_Linear_Coeff_WJG(parameters,data_phasor1d)
% EPTI estimate the coefficient for deghosting
% By Jian Wu MGH 2025
% Fuyixue Wang & Zijing Dong, MGH, 2023
nslice = parameters.nslice;
SMS = parameters.SMS;
slice_sep = nslice/SMS;
DestNx = parameters.DestNx;
echo_num = size(data_phasor1d,2);
pe_lines = size(data_phasor1d,7);
nchannel = size(data_phasor1d,3);
dim = [DestNx,echo_num,nchannel,1,1,1,pe_lines,2,1,slice_sep];
linear_fit_coeff = zeros(2,nchannel,slice_sep);
for slice_index = 1:slice_sep
    data_phasor1d_slice = mean(data_phasor1d(:,:,:,:,:,:,:,:,:,slice_index:slice_sep:end),10);
    linear_fit_coeff_tmp = mrir_artifact_ghost_compute(data_phasor1d_slice);
    linear_fit_coeff_tmp = reshape( median( reshape(linear_fit_coeff_tmp, [2 1 dim(3:7) 1 1 1]), 7) , 2,[]) ; % for saving to use for tailor ghost correction   
    linear_fit_coeff(:,:,slice_index) = linear_fit_coeff_tmp;
end