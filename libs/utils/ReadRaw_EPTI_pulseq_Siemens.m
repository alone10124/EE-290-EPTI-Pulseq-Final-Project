function [twix_obj,EPTIdata,data_phase,parameters] = ReadRaw_EPTI_pulseq_Siemens(filename_EPTI,parameters)

twix_obj = mapVBVD(filename_EPTI);
EPTIdata = single(twix_obj{end}.image.unsorted()) ;
dim = twix_obj{2}.image.sqzSize;
Nx = dim(1);
nchannel = dim(2);
echo_num = parameters.echo_num;
slice = parameters.nslice/parameters.SMS;
rep_num = size(EPTIdata,3)/slice/echo_num;

EPTIdata = reshape(EPTIdata,[Nx,nchannel,echo_num,slice,rep_num]);

data_phase = [];

%update the parameters
parameters.readout = Nx;
parameters.nchannel = nchannel;
parameters.lPhaseEncodingLines = echo_num;
parameters.nechoGE = echo_num;
parameters.iEffectiveEpiEchoSpacing =  (parameters.rampUpTime+parameters.rampDownTime+parameters.flatTopTime)*1e6;
parameters.NslicesEX = parameters.SMS;
% nii header
% data_size = [parameters.Rseg*parameters.Nseg,parameters.reconNx_,parameters.N_slices];
% temp = zeros(data_size);
% [nifti_header]=MyNIFTIWrite(temp,twix_obj{end},'','');
% nifti_header.FrequencyDimension=2;
% nifti_header.PhaseDimension=1;
% nifti_header.SpatialDimension=3;
% parameters.nifti_header = nifti_header;

% opt.SMSRefscan = 0; 
% opt.ReturnStruct=1;    % Zijing settings to read data w/o phasecor navigator
% opt.ReadMultipleRepetitions = 0;
% opt.SqueezeChannels = 0;
% opt.ReturnCellArray = 1; 
% opt.SkipRTFeedback = 0;
% opt.OfflineSave = 0;
% 
% %readin first rep and save prot,evp,(patrefscan,patrefscan_phascor)
% parameters.nifti_header = read_meas_dat_nonav_XA60(filename_EPTI,opt);