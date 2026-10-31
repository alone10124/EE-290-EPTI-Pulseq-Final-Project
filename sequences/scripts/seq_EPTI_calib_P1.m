% This package provides open-source tools for Echo Planar Time-resolved Imaging (EPTI), including sequence, reconstruction, and data processing tools. 
% It is developed using the vendor-neutral Pulseq platform, enabling rapid prototyping and flexible imaging across different scanners.

% !!Access Request: To use this resource, please complete the request form here https://martinos.org/~fw089/sequences.html. 
% This is required for usage tracking in compliance with the requirements of our BRAIN Initiative U24 grant.

% This specific code is the source code for generating the .seq files for:
% Sequence: 2D single-shot EPTI (calibration scan)
% Contrast: gradient-echo contrast
% Field strength: 3T
% Resolution: 3 mm iso

% It is included with pre-configured parameters for the example protocol.
% Users can modify parameters to tailor protocols to their specific needs.

% ------------------------------------------------------------------------------------------------------------
% Version date: 09/2025;
% Developed by: Jian Wu, Zijing Dong, Fuyixue Wang;
% Advanced MRI Acquisition Lab (https://mriacq-lab.martinos.org);
% Updates and additional resources see EPTI dissemination website: https://martinos.org/~fw089/. 
% -----------------------------------------------------------------------------------------------------
%%
%B0 estimation of EPTI
%%!!!!!!! ESP is important!!!!!!!
% This code will check the forbiden ESP, but we still need to make sure that the
% ESP is available in siemens protocol. 
close all
clear 
clc
% warning('off')
addpath(genpath('../libs'))
% parameters
blip_dur = 0.12e-3; % from the EPTI_singleshot, the same as the sequence in IDEA
fov = 216e-3; 
Nx = 72; Ny = 72;            % Define FOV and resolution 
thickness = 3e-3;            % slice thinckness 
sliceGap = 0e-3;
PAT_factor = 3;
ACS_npe = 18; 
Nseg = 49;
Nslices = 42;
minTE = 1;                  % set TE to minimum
TE = 10.1e-3;
minTR = 0;
TR = 1200e-3;
ro_os = 2;                  % oversampling factor
NAV_ETL = 3;                % navigator for deghost
ETL = 14;
ESP = 0.49e-3;
rfDur = 2.05e-3;
RF_angle = 50/180*pi;
B0 = 2.89; 
rfSpoilingInc = 0;                % RF spoiling increment
folder = 'seq/';
mkdir(folder)
sequence_name = 'P1_Calib.seq'; 
% Set system limits
sys = mr.opts('MaxGrad',64,'GradUnit','mT/m',...
    'MaxSlew',200,'SlewUnit','T/m/s',...        
    'rfRingdownTime', 60e-6, 'rfDeadtime', 100e-6,...
    'adcDeadTime', 20e-6, 'adcRasterTime', 2.5e-6, ...
    'gradRasterTime', 10e-6, 'blockDurationRaster', 10e-6, ...
    'B0', B0,'rfRasterTime',1e-6);
sys2 = mr.opts('MaxGrad',64,'GradUnit','mT/m',...
    'MaxSlew',160,'SlewUnit','T/m/s',...        
    'rfRingdownTime', 60e-6, 'rfDeadtime', 100e-6,...
    'adcDeadTime', 20e-6, 'adcRasterTime', 2.5e-6, ...
    'gradRasterTime', 10e-6, 'blockDurationRaster', 10e-6, ...
    'B0', B0,'rfRasterTime',1e-6);
assert((round(blip_dur/sys.gradRasterTime)*sys.gradRasterTime-blip_dur) < 1e-9 )  
% Create fat-sat pulse
sat_ppm=-3.31;
sat_freq=sat_ppm*1e-6*sys.B0*sys.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys2,'Duration',5.12e-3,'dwell',10e-6,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase    
gz_fs = mr.makeTrapezoid('z',sys2,'delay',mr.calcDuration(rf_fs),'duration',2e-3,'Area',-0.1/1e-4); % spoil up to 0.1mm
gz_fs.delay = 0;
% Create slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(RF_angle,'system',sys2,'Duration',rfDur,...
    'SliceThickness',thickness,'apodization',0.3,'timeBwProduct',4,'use','excitation');
% Define gradients
deltak = 1/fov;
kWidth = Nx*deltak;
% Gyblip
blip_dur1 = ceil(2*sqrt(deltak/sys.maxSlew)/sys.gradRasterTime/2)*sys.gradRasterTime*2; % we round-up the duration to 2x the gradient raster time
assert(blip_dur1<=blip_dur+eps);
%Gx
extra_area=blip_dur/2*blip_dur/2*sys.maxSlew; 
gx = mr.makeTrapezoid('x',sys,'Area',kWidth+extra_area,'duration',ESP);
actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx = mr.scaleGrad(gx, kWidth/actual_area);
% calculate ADC
readoutTime = ESP - blip_dur;
adcDwellNyquist = deltak/gx.amplitude/ro_os;
adcDwell = round(adcDwellNyquist/sys.adcRasterTime)*sys.adcRasterTime;
adcSamples = round(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);
time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
adc.delay = round((gx.riseTime+gx.flatTime/2-time_to_center)/sys.rfRasterTime)*sys.rfRasterTime;  %5/9/2024
% Pre-phasing gradients
gxPre = mr.makeTrapezoid('x',sys2,'Area',-gx.area/2);
PEindex = cal_PEindex(ACS_npe,Nseg,PAT_factor);
pe_lines = length(PEindex);
gyPre = mr.makeTrapezoid('y',sys2,'Area',max(abs(PEindex))*deltak,'duration',ESP);
duration_Pre = mr.calcDuration(gxPre,gyPre,gzReph)*2; %extend the duration to lower the PNS
gxPre = mr.makeTrapezoid('x',sys2,'Area',gxPre.area,'duration',duration_Pre);
gzReph = mr.makeTrapezoid('z',sys2,'Area',gzReph.area,'duration',duration_Pre);
% spoil gradients
gxSpoil=mr.makeTrapezoid('x','Area',-2*Nx*deltak,'system',sys2);
gySpoil=mr.makeTrapezoid('y','Area',-2*Ny*deltak,'system',sys2);
gzSpoil=mr.makeTrapezoid('z','Area',-3/thickness,'system',sys2);
duration_spoil = mr.calcDuration(gxSpoil,gySpoil,gzSpoil)*3;%extend the duration to lower the PNS
gxSpoil = mr.makeTrapezoid('x',sys2,'Area',gxSpoil.area,'duration',duration_spoil);
gySpoil = mr.makeTrapezoid('y',sys2,'Area',gySpoil.area,'duration',duration_spoil);
gzSpoil = mr.makeTrapezoid('z',sys2,'Area',gzSpoil.area,'duration',duration_spoil);
gz_fs = mr.makeTrapezoid('z',sys2,'Area',gz_fs.area,'duration',duration_spoil);
% Calculate delay times
ETL_pre = round(ETL/2-1); % PE steps prior to ky=0, excluding the central line
durationToCenter = (ETL_pre+0.5)*mr.calcDuration(gx);
durationofEPI = ETL*mr.calcDuration(gx);
nav_duration = NAV_ETL*mr.calcDuration(gx);
rfCenterInclDelay = rf.delay + mr.calcRfCenter(rf);
if (minTE == 1)
    delayTE1 = sys.blockDurationRaster*2;
    TE = ceil((mr.calcDuration(rf,gz) - rfCenterInclDelay  + duration_Pre + nav_duration + ESP+ durationToCenter)/sys.gradRasterTime)*sys.gradRasterTime;
else
    delayTE1 = ceil((TE - mr.calcDuration(rf,gz) + rfCenterInclDelay - nav_duration - durationToCenter - duration_Pre - duration_Pre)/sys.gradRasterTime)*sys.gradRasterTime;
    TE = ceil(( mr.calcDuration(rf,gz) - rfCenterInclDelay  + duration_Pre + nav_duration + durationToCenter + ESP)/sys.gradRasterTime)*sys.gradRasterTime+delayTE1;
end
assert(delayTE1>=0);
durationSeq = mr.calcDuration(rf_fs)+duration_spoil+mr.calcDuration(rf,gz)+duration_Pre+delayTE1+nav_duration+ESP+durationofEPI+duration_spoil;
if (minTR == 1)
    delayTE2 =  sys.blockDurationRaster*2;
else
    delayTE2 = ceil((TR-Nslices*durationSeq)/Nslices/sys.gradRasterTime)*sys.gradRasterTime;
end
assert(delayTE2>=0);
% calculate slice order
SliceOrder = tdr_sliceorder(Nslices,2);
[~,SliceOrder] = sort(SliceOrder);
% slice positions
excitePositions=(thickness+sliceGap)*((0:(Nslices-1)) - (Nslices-1)/2);
excitePositions = -excitePositions;
excitePositions=excitePositions(SliceOrder); % reorder slices for an interleaved acquisition (optional)
% Define sequence blocks
sys.adcRasterTime = adcDwell;
seq=mr.Sequence(sys);         % Create a new sequence object
seq.addBlock(mr.makeDelay(sys.blockDurationRaster*2));
seq.addBlock(mr.makeLabel('SET','REP', 0));
seq.addBlock(mr.makeLabel('SET','AVG', 29));
rf_phase = 0;
rf_inc = 0;
for nseg = 1:length(PEindex)
    seq.addBlock(mr.makeLabel('SET','SLC', 0));
    actual_area = PEindex(nseg)*deltak;
    scale = actual_area/gyPre.area;
    for s = 1:Nslices
        seq.addBlock(rf_fs);
        seq.addBlock(gz_fs,gxSpoil,gySpoil);
        rf.freqOffset=gz.amplitude*excitePositions(s);
      % RF spoiling
        rf.phaseOffset = rf_phase/180*pi;
        adc.phaseOffset = rf_phase/180*pi;
        rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase = mod(rf_phase+rf_inc, 360.0);

        seq.addBlock(rf,gz);
        seq.addBlock(gxPre,gzReph);
        seq.addBlock(mr.makeDelay(delayTE1));
        seq.addBlock(mr.makeLabel('SET','LIN', 0));
        for i = 1:NAV_ETL
            mr.makeLabel('SET','NAV',1);
            mr.makeLabel('SET','SEG', (-1)^mod(i,2)==1) ;
            mr.makeLabel('SET','REV', (-1)^mod(i+1,2)==1);
            seq.addBlock(mr.scaleGrad(gx,(-1)^mod(i+1,2)),adc);  
            seq.addBlock(mr.makeLabel('INC','LIN', 1));
        end
        seq.addBlock(mr.scaleGrad(gyPre, scale),mr.scaleGrad(gx,(-1)^mod(i+2,2)));
        for i = 1:ETL
            mr.makeLabel('SET','NAV',0);
            mr.makeLabel('SET','SEG', (-1)^mod(i,2)==1) ;
            mr.makeLabel('SET','REV', (-1)^mod(i,2)==1);
            seq.addBlock(mr.scaleGrad(gx,(-1)^mod(i+1,2)),adc);  
            seq.addBlock(mr.makeLabel('INC','LIN', 1));
        end
        seq.addBlock(mr.makeLabel('INC','SLC', 1));  
        seq.addBlock(mr.scaleGrad(gxSpoil,(-1)^mod(i+2,2)),gySpoil,gzSpoil);
        seq.addBlock(mr.makeDelay(delayTE2));
    end
    seq.addBlock(mr.makeLabel('INC','REP', 1));
end
% time checking
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
seq.plot('timeRange',[(durationSeq+delayTE2)*0,(durationSeq+delayTE2)*10]);             % Plot sequence waveforms
seq.setDefinition('FOV', [fov fov max(excitePositions)-min(excitePositions)+thickness]);
seq.setDefinition('Resolution', [fov/Nx fov/Ny thickness]);
seq.setDefinition('Name', sequence_name);
seq.setDefinition('TargetGriddedSamples',Nx); % number of samples after gridding (with oversamping)
seq.setDefinition('TE',TE);  
seq.setDefinition('ETL',ETL+NAV_ETL);  
seq.setDefinition('CNseg',Nseg); 
seq.setDefinition('pe_lines',pe_lines); 
seq.setDefinition('PAT_factor',PAT_factor); 
seq.setDefinition('ReadoutOversamplingFactor',ro_os);
seq.setDefinition('ReceiverGainHigh',1);
seq.setDefinition('TrapezoidGriddingParameters', [gx.riseTime gx.flatTime gx.fallTime adc.delay-gx.delay]); % rise,flat,fall,adc_delay,adc_dur
seq.setDefinition('SlicePositions', excitePositions);      %for 2D
seq.setDefinition('actual_SlicePositions', excitePositions); 
seq.setDefinition('SliceThickness', thickness);           %for 2D
seq.setDefinition('SliceGap', sliceGap);                  %for 2D
seq.write([folder,sequence_name]);

sys_type = 'prisma';  %prisma,Connectome2,skyra  
% check gradient spectrum
esp = mr.calcDuration(gx);
gradSpectrum_bb_v1(seq,sys,sys_type, 3e3, esp);
% check pns
check_pns(seq,sys_type);


