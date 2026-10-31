function [twix_obj,data_temp,data_phasor1d_temp] = ReadRaw_EPTI_calib_XA30_WJG(filename_calib)
% EPTI read calibraiton data
% By Jian Wu MGH 2025
% Fuyixue Wang & Zijing Dong, MGH, 2023

% read the data
twix_obj = mapVBVD(filename_calib);

if length(twix_obj)>1
    twix_obj = twix_obj{2};
end

data = single(twix_obj.image.unsorted()) ;
nslice = twix_obj.hdr.Config.NSlc;
dim = twix_obj.image.sqzSize;
Nx = dim(1);
nchannel = dim(2);
echo_num = floor((cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(38))*twix_obj.hdr.Config.PhaseEncodingLines)/100);
echo_num = echo_num-mod(echo_num,2);
Rseg = ceil(cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(43)));
PAT_factor = ceil(cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(44)));
pe_lines = 18 + ceil((Rseg-18)/PAT_factor);
dim = [Nx,nchannel,echo_num,nslice,pe_lines];
data = reshape(data,dim);
data = permute(data,[1,3,2,5,4]); 
dim = [Nx,echo_num,nchannel,1,1,1,pe_lines,2,1,nslice];
data_temp = zeros(dim);
data_temp(:,1:2:end,:,:,:,:,:,1,:,:) = squeeze(data(:,1:2:end,:,:,:));
data_temp(:,2:2:end,:,:,:,:,:,2,:,:) = data(:,2:2:end,:,:,:);

% read the nav data
data_phasor1d = single(twix_obj.phasecor.unsorted()) ;
dim = [Nx,nchannel,3,nslice,pe_lines];
data_phasor1d = reshape(data_phasor1d,dim);
data_phasor1d = permute(data_phasor1d,[1,3,2,5,4]); 
dim = [Nx,3,nchannel,1,1,1,pe_lines,2,1,nslice];
data_phasor1d_temp = zeros(dim);
data_phasor1d_temp(:,1:2:end,:,:,:,:,:,1,:,:) = squeeze(data_phasor1d(:,1:2:3,:,:,:));
data_phasor1d_temp(:,2:2:end,:,:,:,:,:,2,:,:) = data_phasor1d(:,2,:,:,:);