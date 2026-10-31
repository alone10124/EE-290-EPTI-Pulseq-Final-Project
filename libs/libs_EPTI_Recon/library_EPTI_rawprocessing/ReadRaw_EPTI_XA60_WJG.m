function [twix_obj,EPTIdata,data_phase] = ReadRaw_EPTI_XA60_WJG(filename_EPTI)
% EPTI read calibraiton data
% By Jian Wu MGH 2025
% Fuyixue Wang & Zijing Dong, MGH, 2023

% read the data

twix_obj = mapVBVD(filename_EPTI);
if length(twix_obj)>1
    twix_obj = twix_obj{2};
end

nslice = twix_obj.hdr.Config.NSlc;
SMS = twix_obj.hdr.Meas.MultiBandFactor;
dim = twix_obj.image.sqzSize;
Nx = dim(1);
nchannel = dim(2);
EPTIdata = single(twix_obj.image.unsorted()) ;
Nseg = cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(34));

if(isempty(cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(39))))
    PF_shift = 0;
else
    PF_shift = cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(39));
end
PhaseEncodingLines = twix_obj.hdr.Config.PhaseEncodingLines;
GRE_line = PhaseEncodingLines*cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(39))/100;
if ((isempty(GRE_line)))
    GRE_line=0;
end
if (GRE_line<2)
    GRE_line = 1;
else
    GRE_line = GRE_line-mod(GRE_line,2);
end
if (twix_obj.image.NEco==2)
    SE_line = PhaseEncodingLines*cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(37))/100;
    SE_line = SE_line-mod(SE_line,2);
else
    SE_line = 0;
end
echo_num = GRE_line+SE_line;

rep_num = twix_obj.image.NRep;
EPTIdata = reshape(EPTIdata,[Nx,nchannel,echo_num,nslice/SMS,Nseg*rep_num]);

%read phase
if (isfield(twix_obj, 'phasecor'))
    data_phase = single(twix_obj.phasecor.unsorted()) ;
    dim = [Nx,nchannel,3,nslice/SMS,Nseg*rep_num];
    data_phase = reshape(data_phase,dim);
else
    data_phase=[];
end
