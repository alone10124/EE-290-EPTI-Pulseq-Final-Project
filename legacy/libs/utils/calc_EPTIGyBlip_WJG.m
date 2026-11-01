% calc_EPTIGyBlip_WJG calculates the required phase-encoding gradient "blip" waveform for EPTI 
% acquisitions, depending on the specified state (start, middle, or end).
%
% Inputs:
%   StartStep:    Current phase encoding step at the start of the blip [scalar]
%   EndStep:      Current phase encoding step at the end of the blip [scalar]
%   deltak:       K-space step size in the phase direction [scalar]
%   blip_dur:     Duration of the phase-encoding blip [scalar]
%   gx:           Readout gradient object (used for alignment) [object]
%   sys:          Scanner system structure for gradient shapes [struct]
%   state:        Blip state (1 = start, 2 = middle, 3 = end) [scalar]
%
%
% Outputs:
%   Gy_blip:      Final phase-encoding gradient waveform for the respective state [object]
% 
% Usage:
%   Used within EPTI sequence generation to create phase-encoding blips synchronized to readout,

function Gy_blip = calc_EPTIGyBlip_WJG(StartStep,EndStep,deltak,blip_dur,gx,sys,state)
%state 1 start
%state 2 middle
%state 3 end
if (state ==2)
    gyblip1 = mr.makeTrapezoid('y',sys,'amplitude',StartStep*deltak/blip_dur*2,'riseTime',blip_dur/2, 'flatTime', 0);
    gyblip2 = mr.makeTrapezoid('y',sys,'amplitude',EndStep*deltak/blip_dur*2,'riseTime',blip_dur/2, 'flatTime', 0);
    % split the blip into two halves and produnce a combined synthetic gradient
    gblip1_parts = mr.splitGradientAt(gyblip1, blip_dur/2, sys);
    [gyblip1up, gyblip1down,~]=mr.align('right',gblip1_parts(1),'left',gblip1_parts(2),gx);
    gblip2_parts = mr.splitGradientAt(gyblip2, blip_dur/2, sys);
    [gyblip2up, gyblip2down,~]=mr.align('right',gblip2_parts(1),'left',gblip2_parts(2),gx);
    Gy_blip=mr.addGradients({gyblip1down,gyblip2up }, sys);
elseif (state ==1)
        gyblip1 = mr.makeTrapezoid('y',sys,'amplitude',StartStep*deltak/blip_dur*2,'riseTime',blip_dur/2, 'flatTime', 0);
        gblip1_parts = mr.splitGradientAt(gyblip1, blip_dur/2, sys);
        [Gy_blip, ~,~]=mr.align('right',gblip1_parts(1),'left',gblip1_parts(2),gx);
elseif(state ==3)
        gyblip1 = mr.makeTrapezoid('y',sys,'amplitude',StartStep*deltak/blip_dur*2,'riseTime',blip_dur/2, 'flatTime', 0);
        gblip1_parts = mr.splitGradientAt(gyblip1, blip_dur/2, sys);
        [~, Gy_blip,~]=mr.align('right',gblip1_parts(1),'left',gblip1_parts(2),gx);
end


