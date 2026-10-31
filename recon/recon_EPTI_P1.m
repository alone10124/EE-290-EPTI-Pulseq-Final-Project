% This package provides open-source tools for Echo Planar Time-resolved Imaging (EPTI), including sequence, reconstruction, and data processing tools. 
% It is developed using the vendor-neutral Pulseq platform, enabling rapid prototyping and flexible imaging across different scanners.

% !!Access Request: To use this resource, please complete the request form here https://martinos.org/~fw089/sequences.html. 
% This is required for usage tracking in compliance with the requirements of our BRAIN Initiative U24 grant.

% This specific code is the reconstruction code for:
% Sequence: 2D single-shot EPTI
% Contrast: gradient-echo contrast
% Field strength: 3T
% Resolution: 3 mm iso

% See README file for reconstruction installation and output details.  

% It is included with pre-configured parameters for the example protocol.
% Users can modify parameters to tailor protocols to their specific needs.

% ------------------------------------------------------------------------------------------------------------
% Version date: 09/2025;
% Developed by: Jian Wu, Zijing Dong, Fuyixue Wang;
% Advanced MRI Acquisition Lab (https://mriacq-lab.martinos.org);
% Updates and additional resources see EPTI dissemination website: https://martinos.org/~fw089/. 
% -----------------------------------------------------------------------------------------------------
%%
clear;
clc
addpath(genpath('../libs'))
%% bart path
bart_path = '';                        % path of BART toolbox
setenv('TOOLBOX_PATH', bart_path);
addpath(strcat(getenv('TOOLBOX_PATH'), '/matlab'));
%% write data path and filename
directory = '';           % please put the working directory, data will be saved here
filename_calib = fullfile(directory,'');   % path for rawdata file of the calibration scan
filename_EPTI = fullfile(directory,'');   % path for rawdata file of the imaging scan
save_filename = 'P1_pulseq';                                                        % filename for output data
sequence_name1 = 'seq/P1_Calib.seq';                                                % filename for Pulseq
sequence_name2 = 'seq/P1_EPTI_sms2.seq';                                            % filename for Pulseq
[sub_dir,subdir1_acq,subdir2_recon,subdir3_nii] = EPTI_makeDir(directory,save_filename);
%% set parameters
ISMRMRD_data = 0;	    % '1': ismrmrd data,'0':siemens raw data
start_rep = 4;          % first dynamic to be reconstructed; dynamics before 'start_rep' will be discarded as dummy scans to avoid non-steady-state.
B0_update = 1;          % '0': global b0 drift update using standard 3-line navigator;  '1': (new feature still under development) high-order b0 update using the self-navigators.
moco_recon = 0;         % '0': regular recon; '1': (new feature still under development) motion-aware recon that considers motion during B0 update.
GPU_flag = 1;           % '1': GPU; '0': CPU. GPU is recommended for faster reconstruction.
StartEcho_forPhasecor = 1;                  % start from which echo in calibration data to calculate odd-even echo phase difference for all echoes (for 7T data, recommend skipping early echoes that may have strong eddy-current variation)

% reconstruction parameters
recon_algorithm.regularization = 1;         % regularization for subspace reconstruction. '0': linear, no regularization; '1': locally low-rank
if recon_algorithm.regularization == 1
    recon_algorithm.rho = 1e-2;             % rho for ADMM
    recon_algorithm.lambda = 1e-3;          % lambda for regularization
    recon_algorithm.block_size = 6;         % block size for LLR
end
recon_algorithm.PF_recontype = 1;           % If there is partial fourier in the data, what type of partial fourier reconstruction to use. '0': zero-filling; '1': POCS
recon_algorithm.apodization = 1;            % perform apodization to reduce Gibbs artifact. '1': yes; '0': no
recon_algorithm.apodize_factor = 0.1;       % apodization factor
%% Preparation and preprocessing
parameters = [];
parameters = EPTI_SMS_Ext_parameters_pulseq(sequence_name2,parameters);
% read calibration data
parameters_calib =[];
parameters_calib = EPTI_SMS_Ext_parameters_pulseq(sequence_name1,parameters_calib);
parameters_calib.Rseg = parameters.Rseg;
parameters_calib.Nseg = parameters.Nseg;
parameters_calib.N_phase_encode = parameters.Rseg*parameters.Nseg; 
[data,data_phasor1d,parameters_calib] = ReadRaw_EPTI_calib_pulseq(filename_calib,parameters_calib,ISMRMRD_data);
data = data(:,:,:,:,:,:,:,:,:,parameters_calib.slice_order);
data_phasor1d = data_phasor1d(:,:,:,:,:,:,:,:,:,parameters_calib.slice_order);
% regridding
data_phasor1d = EPI_trapezoid_regridding_WJG(parameters_calib,data_phasor1d);
data = EPI_trapezoid_regridding_WJG(parameters_calib,data);
% read imaging data
[EPTIdata,~,parameters] = ReadRaw_EPTI_imaging_pulseq(filename_EPTI,parameters,ISMRMRD_data);
parameters.nifti_header  =parameters_calib.nifti_header;
[Nx,nchannel,echo_num,~,rep_num] = size(EPTIdata);
parameters_calib.PF_shift = parameters.PF_shift;
parameters_calib.Rseg = parameters.Rseg;
parameters_calib.Nseg = parameters.Nseg;
% estimate the coefficient
linear_fit_coeff = EPTI_SMS_Estimate_Linear_Coeff_WJG(parameters,data_phasor1d);  % estimate the coefficient, using the navigators from the calibration
% odd-even phase correction
kdata = EPTI_SMS_Deghost_WJG(repmat(linear_fit_coeff,[1 1 parameters.SMS]),data); 
kdata = permute(kdata,[2,1,4,3,5]);
% Calculate sensitivity and B0
[sens_map_all,P_dB_all,Phase0_odd_all,Phase0_even_all,gccmtx_aligned,im_GRE_comb] = EPTI_recon_calib_WJG(parameters_calib,kdata);
save(fullfile(subdir1_acq,['PreRecon_SensB0_',save_filename,'.mat']),'sens_map_all','P_dB_all','Phase0_odd_all','Phase0_even_all','gccmtx_aligned','im_GRE_comb','data_phasor1d','-v7.3');
load(fullfile(subdir1_acq,['PreRecon_SensB0_',save_filename,'.mat']))
disp('---------------- EPTI Image Reconstruction start ----------------');
%% Image reconstruction
first_cal = 1; % using to control the first time to calculate something
msk = EPTI_genmask_WJG(im_GRE_comb);
% subspace parameters
BART_savepath = [directory,'/tmp/',save_filename,'_BART_tmp/'];
mkdir(BART_savepath);
K = 4;
parameters.pf_echo = 0;
[Phi,TEs_all,S] = EPTI_gen_Subspace_Basis_GE(parameters,K);      % generate subspace bases: GE
scaling_for_Phi = (S(1:K).^(1/3))';
Phi = Phi.*scaling_for_Phi;

% image matrix size parameters
dt = parameters.iEffectiveEpiEchoSpacing*1e-6;
nx = parameters.DestNx;
ny = parameters.Nseg *parameters.Rseg;
echo_cut1 = 2;
echo_cut2 = 0;
if (parameters.PF_shift>0)  % k-t partial fourier phase shift
    mid=ceil((ny)/2);
    shift = ceil((ny+parameters.PF_shift)/2)-mid;
    linearphase_PF = permute(exp(-1i*2*pi* shift/ny*((0:1:ny-1)-mid)'),[2 1]);
else
    linearphase_PF = 1;
end
[ky_acq,ind] = calc_Pattern_WJG4(parameters.Nseg, parameters.Rseg,echo_num,parameters.Rpe,2,parameters.PF_shift);
scaling_all = calculate_imagescaling_EPTI_useCalib(im_GRE_comb);
scaling_all = scaling_all/mean(scaling_all(:));
[nro,npe,~] = size(msk);
b0_drift= zeros(nro,npe,parameters.nslice,rep_num);

% start image reconstruction
dyn_to_recon = start_rep:rep_num;
for rep = dyn_to_recon
    disp(['-------- reconstructing dynamic #',num2str(rep),'-------- ']);
    EPTIdata_rep = EPTIdata(:,:,:,:,rep);
    EPTIdata_rep = permute(EPTIdata_rep,[1,3,2,5,4]);
    dims = [Nx,echo_num,nchannel,1,1,1,1,2,1,parameters.nslice/parameters.SMS];
    data_temp = zeros(dims);
    data_temp(:,1:2:end,:,:,:,:,:,1,:,:) = squeeze(EPTIdata_rep(:,1:2:end,:,:,:));
    data_temp(:,2:2:end,:,:,:,:,:,2,:,:) = squeeze(EPTIdata_rep(end:-1:1,2:2:end,:,:,:));
    EPTIdata_rep = data_temp(:,:,:,:,:,:,:,:,:,parameters.slice_order);
    % regridding
    EPTIdata_rep = EPI_trapezoid_regridding_WJG(parameters,EPTIdata_rep);
    % odd-even phase correction
    kdata = EPTI_SMS_Deghost_WJG(linear_fit_coeff,EPTIdata_rep); 
    kdata = permute(kdata,[2,1,4,3,5]);
    % CC
    dim = size(kdata);
    dim(4) = size(gccmtx_aligned,2);
    kdata_cc=zeros(dim);
    for t=1:size(kdata,1)
        for echo = 1:size(kdata,6)
            CCDATA_aligned = CC(permute(kdata(t,:,:,:,:,echo),[2 3 5 4 1 6]),gccmtx_aligned, 1);
            kdata_cc(t,:,:,:,:,echo)=permute(CCDATA_aligned,[1 2 4 3]);
        end
    end
    schannel = size(kdata_cc,4);
    % putRawtoPattern
    kdata = putRawtoPattern_SMS_zigzag_PF_selfnav_WJ4(kdata_cc,parameters,3,1);
  
    if(first_cal==1)
        if (B0_update==1)
            warning("!!!!!DO NOT IGNORE THIS MESSAGE: high-order B0 update is only applicable when using sampling pattern 1 during acquisition!!!");
            k_nav_first = kdata_cc;
        else
            k_nav_first = kdata_cc(1:3,:,:,:,:);
        end
        save(fullfile(subdir2_recon,['TEs_all.mat']),'TEs_all','-v7.3');
    end    
    if (B0_update==1)
        k_nav_rep = kdata_cc;
    else
        k_nav_rep = kdata_cc(1:3,:,:,:,:);
    end
    if (B0_update==1)
        assert(parameters.SMS<=2);
        if (first_cal==1)
            lowres_refI= recon_lowresI2(k_nav_first,ky_acq,sens_map_all,ind,4, parameters,linearphase_PF);
        end
        lowresI_rep= recon_lowresI2(k_nav_rep,ky_acq,sens_map_all,ind,4, parameters,linearphase_PF);

        % b0 update that also considers motion
        if (first_cal==1)||(moco_recon==1)
            [motion_param_file,~] = motion_estimation_wjg(lowresI_rep(:,:,:,2),lowres_refI(:,:,:,2),BART_savepath);
        end
        msk2 = warp_wjg(msk,motion_param_file,BART_savepath);
        B0_drift_Hz = DriftB0_esti_nav_WJG3_4(lowresI_rep,lowres_refI,dt,ind,msk2,motion_param_file,BART_savepath); % field drift/global B0 change estimation

        b0_drift(:,:,:,rep) = -B0_drift_Hz;
        P_dB_all2 = warp_wjg(P_dB_all,motion_param_file,BART_savepath); % register to TRN
        P_dB_all2 = P_dB_all2 - B0_drift_Hz;
    else
        B0_drift_Hz = DriftB0_esti_nav(k_nav_first,k_nav_rep,dt); % global field drift/B0 change estimation
        B0_drift_Hz = repmat(B0_drift_Hz,[parameters.SMS 1]);
        P_dB_all2 = P_dB_all+permute(B0_drift_Hz,[2 3 1]);
    end   
    slice_sep = parameters.nslice/parameters.SMS;
    scaling_dyn = calculate_imagescaling_EPTIv2(kdata).*scaling_all*2/parameters.SMS; % adjust  the scale
    im_epti = zeros(nx,ny+parameters.PF_shift,parameters.nslice,echo_num-echo_cut1-echo_cut2);
    fprintf('slice: ');
    for dif_slice=1:slice_sep          % slice loop
        fprintf('%d ',dif_slice);
        slice_group = dif_slice:slice_sep:slice_sep*parameters.SMS;
        slice_group = slice_group(end:-1:1);
        kdata_tmp = kdata(:,:,:,:,:,dif_slice);
        scaling = max(scaling_dyn(slice_group));
        sens_map = sens_map_all(:,:,slice_group,:);
        P_dB = P_dB_all2(:,:,slice_group);
        Phase0_odd = Phase0_odd_all(:,:,:,slice_group);
        Phase0_even = Phase0_even_all(:,:,:,slice_group);
        Phase0_odd=exp(1i*Phase0_odd).*linearphase_PF;
        Phase0_even=exp(1i*Phase0_even).*linearphase_PF;
        Phase0_use_t = zeros(size(Phase0_odd,1),size(Phase0_odd,2),size(Phase0_odd,4),size(kdata_tmp,1));
        Phase0_use_t(:,:,:,1:2:end) = permute(repmat(exp(1i*angle(mean(Phase0_odd(:,:,min(StartEcho_forPhasecor,size(Phase0_odd,3)):end,:),3))),[1 1 length(1:2:size(Phase0_use_t,4)) 1]),[1 2 4 3]);
        Phase0_use_t(:,:,:,2:2:end) = permute(repmat(exp(1i*angle(mean(Phase0_even(:,:,min(StartEcho_forPhasecor,size(Phase0_even,3)):end,:),3))),[1 1 length(2:2:size(Phase0_use_t,4)) 1]),[1 2 4 3]);
        [im_recon,a] = SubspaceRecon_Process_v5(kdata_tmp,sens_map,P_dB,...
            Phase0_use_t,TEs_all,Phi,parameters,scaling,BART_savepath,recon_algorithm,GPU_flag); 
        % figure; imshow3(sos(im_recon(:,:,:,11:end-10),4),[]);
        im_epti(:,:,dif_slice:slice_sep:end,:)=im_recon(end:-1:1,end:-1:1,end:-1:1,echo_cut1+1:end-echo_cut2);
    end
    save(fullfile(subdir2_recon,['Recon_EPTI_',save_filename,'_Dyn_',num2str(rep),'.mat']),'im_epti','-v7.3');
    first_cal = 0;
end
disp('---------------- image reconstruction finished ----------------')
%% EPTI data preprocessing 
disp('---------------- image preprocessing start ----------------')
nDyn = length(dyn_to_recon);
for loopi = 1:nDyn
    load(fullfile(subdir2_recon,['Recon_EPTI_',save_filename,'_Dyn_',num2str(dyn_to_recon(loopi)),'.mat']));
    if (loopi==1)
        [dim1,dim2,dim3,necho] = size(im_epti);
        im_epti_final = zeros(dim1,dim2,dim3,nDyn,necho);
    end
    im_epti_final(:,:,:,loopi,:) = im_epti;
end
im_epti_final = im_epti_final*10000;
im_epti_sos = squeeze(sos(im_epti_final,5));
saveniftidata_EPTI_pulseq(im_epti_sos,parameters,fullfile(subdir3_nii,[save_filename,'_im_sos.nii']));

im_epti_Tave = squeeze(mean(abs(im_epti_final),4));
saveniftidata_EPTI_pulseq(im_epti_Tave,parameters,fullfile(subdir3_nii,[save_filename,'_im_Tave_echoes.nii']));

for dif = 1:necho
    disp(['---------------------saving echo',num2str(dif)]);
    saveniftidata_EPTI_pulseq(squeeze(abs(im_epti_final(:,:,:,:,dif))),parameters,fullfile(subdir3_nii,'echoes',[save_filename,'_im_echo',num2str(dif),'.nii']));
end
disp('---------------- T2 mapping ----------------')
load([subdir2_recon,'/TEs_all.mat'])
TEs_all_save = TEs_all(echo_cut1+1:end-echo_cut2);
writematrix(TEs_all_save,fullfile(subdir3_nii,[save_filename,'__TEs.txt']));
im_epti_final = permute(im_epti_final,[1,2,3,5,4]);
% combined image
[img_comb,T2s_map] = optimalCombineEchoes(im_epti_final,TEs_all_save);
saveniftidata_EPTI_pulseq(T2s_map,parameters,fullfile(subdir3_nii,[save_filename,'__T2s.nii']));
saveniftidata_EPTI_pulseq(img_comb,parameters,fullfile(subdir3_nii,[save_filename,'_im_comb.nii'])); 
disp('done!');
