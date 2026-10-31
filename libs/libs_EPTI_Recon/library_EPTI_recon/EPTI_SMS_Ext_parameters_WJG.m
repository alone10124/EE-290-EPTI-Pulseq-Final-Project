function parameters = EPTI_SMS_Ext_parameters_WJG(twix_obj)
% EPTI extract parameters
% By Jian Wu MGH 2025
% Fuyixue Wang & Zijing Dong, MGH, 2023
dim = twix_obj.image.sqzSize;
DestNx = twix_obj.hdr.Meas.lBaseResolution;
nslice = twix_obj.hdr.Config.NSlc;
Nx = dim(1);

parameters=[];
parameters.nslice = twix_obj.hdr.Config.NSlc;
parameters.SMS = twix_obj.hdr.Meas.MultiBandFactor;
parameters.DestNx = twix_obj.hdr.Meas.lBaseResolution;
parameters.dim = twix_obj.image.sqzSize;
parameters.Nx =parameters.dim(1);
parameters.readout = Nx;
parameters.rampUpTime = twix_obj.hdr.Meas.alRegridRampupTime(1)*1e-6;
parameters.rampDownTime = twix_obj.hdr.Meas.alRegridRampdownTime(1)*1e-6;
parameters.acqDelayTime = twix_obj.hdr.Meas.alRegridDelaySamplesTime(1)*1e-6;
parameters.flatTopTime = twix_obj.hdr.Meas.alRegridFlattopTime(1)*1e-6;
parameters.dwellTime = twix_obj.hdr.Meas.alDwellTime(1)*1e-9;
parameters.N_phase_encode = DestNx;
parameters.reconNx_ = DestNx;
parameters.FOV_1 = twix_obj.hdr.Config.RoFOV;
parameters.position = [0,0,0];
parameters.read_dir = [0,0,0];
parameters.N_slices = nslice;
parameters.DestNx = DestNx;
parameters.alTE = twix_obj.hdr.Meas.alTE;
parameters.PAT_factor = cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(44));
parameters.CNseg = cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(43));

parameters.Rseg = cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(33));
parameters.Nseg = cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(34));
parameters.Rpe = cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(35));
parameters.PF_shift = cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(39));
parameters.rep_num = twix_obj.image.NRep;
parameters.NslicesEX = twix_obj.hdr.Meas.MultiBandFactor;
parameters.SMS_shift = twix_obj.hdr.Meas.MultiBandFactor;
% parameters.iEffectiveEpiEchoSpacing = (parameters.rampUpTime+parameters.rampDownTime+parameters.flatTopTime)*1e6;
parameters.iEffectiveEpiEchoSpacing = twix_obj.hdr.Meas.iEffectiveEpiEchoSpacing;

parameters.alTE(1) = parameters.alTE(1)-1000-parameters.iEffectiveEpiEchoSpacing;  
if (isempty(parameters.PF_shift ))
    parameters.PF_shift  = 0;
    parameters.PF_shift = 0;
end

if(isempty(cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(39))))
    PF_shift = 0;
else
    PF_shift = cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(39));
end

PhaseEncodingLines = twix_obj.hdr.Config.PhaseEncodingLines;
GRE_line = PhaseEncodingLines*cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(38))/100;
GRE_line = GRE_line-mod(GRE_line,2);
if (twix_obj.image.NEco==2)
    SE_line = PhaseEncodingLines*cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(36))/100;
    SE_line = SE_line-mod(SE_line,2);
else
    SE_line = 0;
end
parameters.GRE_line = GRE_line;
parameters.SE_line = SE_line;
%acquisiton pattern
sampling_pattern = twix_obj.image.Lin;
try 
    sampling_pattern = reshape(sampling_pattern,[GRE_line+SE_line,parameters.nslice/parameters.SMS,parameters.Nseg*parameters.rep_num]);
    sampling_pattern_GE = squeeze(sampling_pattern(1:GRE_line,1,1:parameters.Nseg));
    parameters.sampling_pattern_GE = sampling_pattern_GE';
    sampling_pattern_SE = squeeze(sampling_pattern(GRE_line+1:GRE_line+SE_line,1,1:parameters.Nseg));
    parameters.sampling_pattern_SE = sampling_pattern_SE';
catch ME % for calibration
    parameters.sampling_pattern_GE = [];
    parameters.sampling_pattern_SE = [];
end
