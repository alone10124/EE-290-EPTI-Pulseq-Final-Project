% cal_PEindex  Calculate phase-encoding line indices
%
% Inputs:
%   ACS_npe     - Number of ACS (auto-calibration signal) phase-encoding lines
%   Nseg        - Total number of phase-encoding lines (full k-space lines)
%   PAT_factor  - Parallel imaging acceleration factor (e.g., GRAPPA)
%
% Output:
%   index       - Array of phase-encoding line indices (centered around 0)

function index = cal_PEindex(ACS_npe,Nseg,PAT_factor)
n_pe = ACS_npe + ceil((Nseg-ACS_npe)/PAT_factor);
index = zeros(1,n_pe);
for i = 1:n_pe
    if i<(ceil((Nseg-ACS_npe)/PAT_factor)/2)+1
        index(i) = -(Nseg-1)/2 + PAT_factor*(i-1);
    elseif i >= ((ceil((Nseg-ACS_npe)/PAT_factor)/2) + ACS_npe + 1)
        index(i) = -(Nseg-1)/2 + PAT_factor*(i-ACS_npe+round(ACS_npe/PAT_factor)-1);
    else
        index(i) = -(ACS_npe)/2 + (i-floor(ceil((Nseg-ACS_npe)/PAT_factor)/2)-1);
    end
end