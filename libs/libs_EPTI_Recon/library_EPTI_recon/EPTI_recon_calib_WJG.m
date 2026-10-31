function [sens_map_all,P_dB_all,Phase0_odd_all,Phase0_even_all,gccmtx_aligned,im_GRE_comb] = EPTI_recon_calib_WJG(parameters,kdata)
% EPTI reconstruct the calibration data
% By Jian Wu MGH 2025
% Fuyixue Wang & Zijing Dong, MGH, 2023

[echo_num,DestNx,~,nchannel,nslice] = size(kdata);
espirit_num_acs = 24;                       % ESPIRiT parameters for 3D coil sensitivity estimation, calbiration size, should be larger than # of slices
espirit_c = 0.1;                            % ESPIRiT parameters for 3D coil sensitivity estimation, threshold for eigenvalue
parameters_calib = [];
meas = [];
% parameters_calib.iEffectiveEpiEchoSpacing = (parameters.rampUpTime+parameters.rampDownTime+parameters.flatTopTime)*1e6;
parameters_calib.iEffectiveEpiEchoSpacing = parameters.iEffectiveEpiEchoSpacing;
parameters_calib.alTE = parameters.alTE; 
parameters_calib.nechoGE = echo_num;
parameters_calib.sWipMemBlock_alFree(35) = parameters.PAT_factor;
parameters_calib.sWipMemBlock_alFree(34)=parameters.CNseg;
parameters_calib.sSliceArray_lSize = nslice;
meas.prot.iMaxNoOfRxChannels = nchannel;
meas.prot.lBaseResolution = DestNx;
meas.prot.Rseg = parameters.Rseg;
meas.prot.Nseg = parameters.Nseg;
[sens_map_all,P_dB_all,Phase0_odd_all,Phase0_even_all,gccmtx_aligned,im_GRE_comb] = EPTI_Preprocessing_B0_sens_BARTv5(kdata,meas.prot,parameters_calib,espirit_num_acs,espirit_c);
