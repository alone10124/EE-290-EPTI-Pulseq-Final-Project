function kdata = EPTI_SMS_Deghost_WJG(linear_fit_coeff,data)
% EPTI read deghost
% By Jian Wu MGH 2025
% Fuyixue Wang & Zijing Dong, MGH, 2023
[DestNx,echo_num,nchannel,~,~,~,pe_lines,~,~,nslice] = size(data);
kdata = zeros([DestNx,echo_num,nchannel,pe_lines,nslice]);
for slice_index = 1:nslice
    temp_coeff = linear_fit_coeff(:,:,slice_index);
    tempdata = data(:,:,:,:,:,:,:,:,:,slice_index);
    data_hybrid = mrir_iDFT_freqencode(tempdata); 
    data_hybrid = mrir_artifact_ghost_correct(data_hybrid, temp_coeff);
    data_hybrid = sum(data_hybrid,8);
    data_hybrid = mrir_fDFT_freqencode(data_hybrid); 
    kdata(:,:,:,:,slice_index) = squeeze(data_hybrid);
end