function [data_temp,data_phasor1d_temp,parameters] = ReadRaw_EPTI_calib_pulseq_Siemens(filename_calib,parameters)


twix_obj = mapVBVD(filename_calib);
data = single(twix_obj{end}.image.unsorted()) ;
dim = twix_obj{2}.image.sqzSize;
Nx = dim(1);
nchannel = dim(2);
echo_num = length(unique(twix_obj{1,2}.image.Lin));
data = reshape(data,[Nx,nchannel,echo_num,parameters.N_slices,parameters.pe_lines]);
data_phasor1d = data(:,:,1:3,:,:);
data_phasor1d = squeeze(data_phasor1d);
data_phasor1d = permute(data_phasor1d,[1,3,2,5,4]); 
dims = [Nx,3,nchannel,1,1,1,parameters.pe_lines,2,1,parameters.N_slices];
data_phasor1d_temp = zeros(dims);
data_phasor1d_temp(:,1:2:end,:,:,:,:,:,1,:,:) = squeeze(data_phasor1d(:,1:2:3,:,:,:));
data_phasor1d_temp(:,2:2:end,:,:,:,:,:,2,:,:) = data_phasor1d(end:-1:1,2,:,:,:);

data = data(:,:,4:end,:,:);
data = squeeze(data);
data = squeeze(permute(data,[1,3,2,5,4])); 
dims = [Nx,echo_num-3,nchannel,1,1,1,parameters.pe_lines,2,1,parameters.N_slices];
data_temp = zeros(dims);
data_temp(:,1:2:end,:,:,:,:,:,1,:,:) = squeeze(data(:,1:2:end,:,:,:));
data_temp(:,2:2:end,:,:,:,:,:,2,:,:) = squeeze(data(end:-1:1,2:2:end,:,:,:));

%update the parameters
parameters.readout = Nx;
parameters.nchannel = nchannel;
if(isempty(parameters.Rseg ))
    parameters.Rseg = parameters.DestNx ;
    parameters.Nseg = 1;
end

parameters.iEffectiveEpiEchoSpacing =  (parameters.rampUpTime+parameters.rampDownTime+parameters.flatTopTime)*1e6;
% nii header
% data_size = [parameters.reconNx_,parameters.Rseg*parameters.Nseg, parameters.N_slices];
% temp = zeros(data_size);
% [nifti_header]=MyNIFTIWrite(temp,twix_obj{end},'','');
% nifti_header.FrequencyDimension=1;
% nifti_header.PhaseDimension=2;
% nifti_header.SpatialDimension=3;
% parameters.nifti_header = nifti_header;


opt.SMSRefscan = 0; 
opt.ReturnStruct=1;    % Zijing settings to read data w/o phasecor navigator
opt.ReadMultipleRepetitions = 0;
opt.SqueezeChannels = 0;
opt.ReturnCellArray = 1; 
opt.SkipRTFeedback = 0;
opt.OfflineSave = 0;
meas= read_meas_dat_nonav_XA60(filename_calib,opt);
meas.evp.NImageCols = parameters.reconNx_;
meas.evp.NImageLins = parameters.N_phase_encode;
%readin first rep and save prot,evp,(patrefscan,patrefscan_phascor)
parameters.nifti_header = meas;