% This package provides open-source tools for Echo Planar Time-resolved Imaging (EPTI), including sequence, reconstruction, and data processing tools. 
% It is developed using the vendor-neutral Pulseq platform, enabling rapid prototyping and flexible imaging across different scanners.

% !!Access Request: To use this resource, please complete the request form here https://martinos.org/~fw089/sequences.html. 
% This is required for usage tracking in compliance with the requirements of our BRAIN Initiative U24 grant.

% This specific code is the source code for generating the .seq files for:
% Sequence: 2D single-shot EPTI (imaging scan)
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
close all
clear 
clc
% warning('off')
addpath(genpath('../libs'))
% parameters
blip_dur = 0.12e-3;
fovx = 216e-3; fovy = 216e-3; 
Nx = 72; Ny = 72;           % Define FOV and resolution  
thickness = 3e-3;           % slice thinckness  
MB = 2;                     % MB factors, only support the MB <= 4
Nslices = 42;
sliceGap = 0;
sliceSep = (thickness+sliceGap)*(Nslices/MB);
minTE = 1;                  % set TE to minnum
minTR = 0;
TE1 = 4e-3;                 % TE of the first pe line
TR = 1800e-3;
ro_os = 2;                  % oversampling factor
NAV_ETL = 3;                % navigator for deghost
ETL = 128-NAV_ETL;
Rseg = 72;
Rpe = 6;
Nseg = 1;
PFshift = 0;
ESP = 0.49e-3;
RF_angle = 60;
rfDur = 4e-3;
rfTB = 4;                   %TBP of the RF pulse
B0 = 2.89; 
folder = 'seq/';
mkdir(folder)
sequence_name = ['P1_EPTI_sms',num2str(MB),'.seq'];
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
    'B0', B0,'rfRasterTime',1e-6);% different module use the different restriction
assert((round(blip_dur/sys.gradRasterTime)*sys.gradRasterTime-blip_dur) < 1e-9 )  
assert(mod(Nslices,MB) == 0)   
% Create fat-sat pulse 
sat_ppm=-3.31;
sat_freq=sat_ppm*1e-6*sys.B0*sys.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys2,'Duration',5.12e-3,'dwell',10e-6,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase    
gz_fs = mr.makeTrapezoid('z',sys2,'delay',mr.calcDuration(rf_fs),'duration',2e-3,'Area',-0.1/1e-4); % spoil up to 0.1mm
gz_fs.delay = 0;
% SMS RF Create slice selection pulse and gradient
sysGE = toppe.systemspecs('maxGrad', 0.5*sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', 0.4*sys.maxSlew/sys.gamma/10, ...           % G/cm/ms
    'maxRF', 0.25);
[rf, gzRF, freq] = getsmspulse(RF_angle, thickness, rfTB, rfDur, ...
MB, sliceSep, sysGE, sys2, 'doSim', 0, 'type', 'st', 'ftype', 'ls'); 
rf.use='excitation';
%seperate the gzRF
[~,ind1] = max(gzRF.waveform);
[~,ind2] = min(gzRF.waveform);
temp_signal = abs(gzRF.waveform(ind1:ind2));
[~,ind3] = min(temp_signal);
ind4 = ind1+ind3-1;
amplitude = gzRF.waveform(ind1);
duration = ind4*sys.gradRasterTime;
FlatTime = (ind3-ind1)*sys.gradRasterTime;
gz = mr.makeTrapezoid('z', sys, ...  % readout trapezoid, temporary object
    'Amplitude', amplitude, ...
    'FlatTime', FlatTime, ...
    'riseTime',(ind1-1)*sys.gradRasterTime, ...
    'fallTime',(ind1-1)*sys.gradRasterTime);
gzReph = mr.makeTrapezoid('z', sys2, ...  % readout trapezoid, temporary object
    'Area', -gz.area/2, ...
    'Duration',1.2e-3);

% Define gradients
deltakx = 1/fovx;
deltaky = 1/fovy;
kWidth = Nx*deltakx;
% Gyblip + Gzblip duration
[ky_acq,~] = calc_Pattern_WJG3(Nseg,Rseg,ETL+NAV_ETL,Rpe,2,PFshift);
ky_step = diff(ky_acq);
ky_step = ky_step(3:end);
GyblipAreaMax = deltaky*max(abs(ky_step));
blip_dur1 = ceil(2*sqrt(GyblipAreaMax/sys.maxSlew)/sys.gradRasterTime/2)*sys.gradRasterTime*2; % we round-up the duration to 2x the gradient raster time
[kz_step,deltakz] = calc_Pattern_SMS_blip(MB,ETL,sliceSep);
gzarea = deltakz*max(abs(kz_step)); 
blip_dur2 = ceil(2*sqrt(gzarea/sys.maxSlew)/sys.gradRasterTime/2)*sys.gradRasterTime*2; 
blip_dur3 = max(blip_dur1,blip_dur2);
assert(blip_dur3<=blip_dur+eps);
% Gx
extra_area=blip_dur/2*blip_dur/2*sys.maxSlew; 
gx = mr.makeTrapezoid('x',sys,'Area',kWidth+extra_area,'duration',ESP);
actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx = mr.scaleGrad(gx, kWidth/actual_area);
% ADC
readoutTime = ESP - blip_dur;
adcDwellNyquist = deltakx/gx.amplitude/ro_os;
adcDwell = round(adcDwellNyquist/sys.adcRasterTime)*sys.adcRasterTime;
adcSamples = round(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);
time_to_center = adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
adc.delay = round((gx.riseTime+gx.flatTime/2-time_to_center)/sys.rfRasterTime)*sys.rfRasterTime;  %5/9/2024

% Pre-phasing gradients
gxPre = mr.makeTrapezoid('x',sys2,'Area',-gx.area/2);
duration_Pre = mr.calcDuration(gzReph);
gzReph = mr.makeTrapezoid('z',sys2,'Area',gzReph.area,'duration',duration_Pre);
gyPre = mr.makeTrapezoid('y',sys2,'Area',0*deltaky);
gxPre = mr.makeTrapezoid('x',sys2,'Area',gxPre.area,'duration',duration_Pre);
gyPre = mr.makeTrapezoid('y',sys2,'Area',gyPre.area,'duration',duration_Pre);
% spoil gradients
gxSpoil=mr.makeTrapezoid('x','Area',-2*Nx*deltakx,'system',sys2);
gySpoil=mr.makeTrapezoid('y','Area',-2*Ny*deltaky,'system',sys2);
gzSpoil=mr.makeTrapezoid('z','Area',-4/thickness,'system',sys2);
duration_spoil = mr.calcDuration(gxSpoil,gySpoil,gzSpoil)*2;%extend the duration tp lower the PNS
gxSpoil = mr.makeTrapezoid('x',sys2,'Area',gxSpoil.area,'duration',duration_spoil);
gySpoil = mr.makeTrapezoid('y',sys2,'Area',gySpoil.area,'duration',duration_spoil);
gzSpoil = mr.makeTrapezoid('z',sys2,'Area',gzSpoil.area,'duration',duration_spoil);
gz_fs = mr.makeTrapezoid('z',sys2,'Area',gz_fs.area,'duration',duration_spoil);
% Calculate delay times 
durationToCenter = 0.5*mr.calcDuration(gx);
durationofEPI = (ETL+NAV_ETL)*mr.calcDuration(gx);
rfCenterInclDelay = rf.delay + mr.calcRfCenter(rf);
if minTE == 1
    delayTE1 = sys.blockDurationRaster*2;
    TE1 = mr.calcDuration(rf,gz)-rfCenterInclDelay+mr.calcDuration(gxPre,gyPre,gzReph)+ ...
    sys.blockDurationRaster*2+((ETL+NAV_ETL)/2-0.5)*ESP;
else
    delayTE1 = ceil((TE1 - mr.calcDuration(rf,gz) + rfCenterInclDelay - ((ETL+NAV_ETL)/2-0.5)*ESP-duration_Pre)/sys.gradRasterTime)*sys.gradRasterTime;
end
assert(delayTE1>=0);
trig_out=mr.makeDigitalOutputPulse('ext1','duration', 100e-6,'delay',500e-6); % possible channels: 'osc0','osc1','ext1'
durationSeq = mr.calcDuration(rf_fs)+duration_spoil+mr.calcDuration(rf,gz)+duration_Pre+delayTE1+durationofEPI+duration_spoil;
if (minTR == 1)
    delayTE2 =  sys.blockDurationRaster*2;
else
    actual_slice = Nslices/MB;
    delayTE2 = ceil((TR-mr.calcDuration(trig_out)-actual_slice*durationSeq)/actual_slice/sys.gradRasterTime)*sys.gradRasterTime;
end
assert(delayTE2>=0);
% calculate slice order
SliceOrder = tdr_sliceorder(Nslices/MB,2);
[~,SliceOrder] = sort(SliceOrder);
% slice positions
excitePositions=(thickness+sliceGap)*((0:(Nslices/MB-1)) - (Nslices/MB-1)/2);
excitePositions = -excitePositions;
excitePositions=excitePositions(SliceOrder); % reorder slices for an interleaved acquisition (optional)
% Define sequence blocks
sys.adcRasterTime = adcDwell;
seq=mr.Sequence(sys);         % Create a new sequence object
seq.addBlock(mr.makeDelay(sys.blockDurationRaster*2));
seq.addBlock(mr.makeLabel('SET','SLC', 0));
seq.addBlock(trig_out);
for s=1:Nslices/MB
    seq.addBlock(rf_fs);
    seq.addBlock(gz_fs,gxSpoil,gySpoil);
    rf.freqOffset = gz.amplitude*excitePositions(s);    
    seq.addBlock(rf,gz);
    seq.addBlock(gxPre,gyPre,gzReph);  
    seq.addBlock(mr.makeDelay(delayTE1));
    seq.addBlock(mr.makeLabel('SET','LIN', 0));
    phi = 0;
    adc.phaseOffset = 0;
    for i = 1:NAV_ETL
        mr.makeLabel('SET','NAV',1);
        mr.makeLabel('SET','REV', (-1)^mod(i+1,2)==1);
        mr.makeLabel('SET','SEG', (-1)^mod(i,2)==1) ;
        if i == NAV_ETL
            Gy_blip_now = calc_EPTIGyBlip_WJG(ky_step(1),ky_step(1),deltaky,blip_dur,gx,sys,1);
            Gz_blip_now = calc_EPTIGzBlip_WJG(kz_step(1),kz_step(1),deltakz,blip_dur,gx,sys,1);
            seq.addBlock(mr.scaleGrad(gx,(-1)^mod(i+1,2)),Gy_blip_now,Gz_blip_now,adc);  
            seq.addBlock(mr.makeLabel('INC','LIN', 1));  
        else
            seq.addBlock(mr.scaleGrad(gx,(-1)^mod(i+1,2)) ,adc);  
            seq.addBlock(mr.makeLabel('INC','LIN', 1));
        end
    end
    for i = 1:ETL
        temp_gzarea = deltakz*kz_step(i);
        phi = phi+2*pi*temp_gzarea*(excitePositions(s)+sliceSep*(MB-1)/2); 
        adc.phaseOffset = phi;
        mr.makeLabel('SET','NAV',0);
        mr.makeLabel('SET','REV', (-1)^mod(i,2)==1);
        mr.makeLabel('SET','SEG', (-1)^mod(i,2)==1) ;
        if i == ETL
             Gy_blip_now = calc_EPTIGyBlip_WJG(ky_step(i),ky_step(i),deltaky,blip_dur,gx,sys,3);
             Gz_blip_now = calc_EPTIGzBlip_WJG(kz_step(i),kz_step(i),deltakz,blip_dur,gx,sys,3);
             seq.addBlock(mr.scaleGrad(gx,(-1)^mod(i,2)),Gy_blip_now,Gz_blip_now,adc);
        else
            Gy_blip_now = calc_EPTIGyBlip_WJG(ky_step(i),ky_step(i+1),deltaky,blip_dur,gx,sys,2);
            Gz_blip_now = calc_EPTIGzBlip_WJG(kz_step(i),kz_step(i+1),deltakz,blip_dur,gx,sys,2);
            seq.addBlock(mr.scaleGrad(gx,(-1)^mod(i,2)),Gy_blip_now,Gz_blip_now,adc);
        end
        seq.addBlock(mr.makeLabel('INC','LIN', 1));
    end
    seq.addBlock(mr.makeLabel('INC','SLC', 1));
    seq.addBlock(gxSpoil,gySpoil,gzSpoil);
    seq.addBlock(mr.makeDelay(delayTE2));
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
% in order to shift the FOV
SliceOrder = tdr_sliceorder(Nslices,2);
[~,SliceOrder] = sort(SliceOrder);
slicePositions=(thickness+sliceGap)*((0:(Nslices-1)) - (Nslices-1)/2);
slicePositions=-slicePositions(SliceOrder); % reorder slices for an interleaved acquisition (optional)

seq.plot('timeRange',[(durationSeq+delayTE2)*0,(durationSeq+delayTE2)*3]);             % Plot sequence waveforms
seq.setDefinition('FOV', [fovx fovy max(slicePositions)-min(slicePositions)+thickness]);
seq.setDefinition('Resolution', [fovx/Nx fovy/Ny thickness]);
seq.setDefinition('Name', sequence_name);
seq.setDefinition('TargetGriddedSamples',Nx); % number of samples after gridding (with oversamping)
seq.setDefinition('TE',TE1);  
seq.setDefinition('GRE_Lines',ETL+NAV_ETL);  
seq.setDefinition('Nseg',Nseg);
seq.setDefinition('Rpe',Rpe);  
seq.setDefinition('Rseg',Rseg);  
seq.setDefinition('PFshift',PFshift);  
seq.setDefinition('SMS',MB);  
seq.setDefinition('pattern_GE',ky_acq);  
seq.setDefinition('ReadoutOversamplingFactor',ro_os);
seq.setDefinition('ReceiverGainHigh',1);
seq.setDefinition('TrapezoidGriddingParameters', [gx.riseTime gx.flatTime gx.fallTime adc.delay-gx.delay]); % rise,flat,fall,adc_delay,adc_dur
seq.setDefinition('SlicePositions', slicePositions);      %for 2D
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
