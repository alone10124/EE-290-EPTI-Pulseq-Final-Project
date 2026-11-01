function make_reference_fixtures()
% MAKE_REFERENCE_FIXTURES  Generate MATLAB reference outputs for the Python port's tests. Run once.
%
%   >> cd legacy; make_reference_fixtures
%
% It writes three files into tests/ref/:
%   P1_Calib.seq, P1_EPTI_sms2.seq   from the unchanged MGH scripts in legacy/sequences/ (Siemens Prisma settings)
%   fixtures.mat                     helper outputs + black-box input/output pairs for the encrypted P-code
%
% Requirements: MATLAB with the Signal Processing and Curve Fitting toolboxes (used by getsmspulse/makeslr).
% You do NOT need scanner data, BART, or AFNI. The script opens a few sequence-diagram/PNS figures; just close them.
%
% Written for the EE 290 EPTI Python port (see docs/EPTI_PORT_PLAN.md). The Python tests in tests/ load fixtures.mat
% with scipy.io.loadmat, so everything is saved with -v7.

here   = fileparts(mfilename('fullpath'));
refdir = fullfile(here, '..', 'tests', 'ref');
if ~exist(refdir, 'dir'), mkdir(refdir); end
addpath(genpath(fullfile(here, 'libs')));

%% 1) Reference .seq files from the original scripts
% The scripts write into ./seq/ relative to the current folder, and call addpath(genpath('../libs')).
seqdir = fullfile(here, 'sequences');
run_script_isolated(seqdir, 'seq_EPTI_calib_P1.m');
run_script_isolated(seqdir, 'seq_EPTI_imaging_P1.m');
copyfile(fullfile(seqdir, 'seq', 'P1_Calib.seq'),     refdir);
copyfile(fullfile(seqdir, 'seq', 'P1_EPTI_sms2.seq'), refdir);

F = struct();
F.matlab_version = version;
F.created = char(datetime('now'));

%% 2) Sequence helpers (arguments as used in seq_EPTI_imaging_P1.m / seq_EPTI_calib_P1.m / recon_EPTI_P1.m)
[F.ky_WJG3, F.ind_WJG3] = calc_Pattern_WJG3(1, 72, 128, 6, 2, 0);   % imaging L109
[F.ky_WJG4, F.ind_WJG4] = calc_Pattern_WJG4(1, 72, 128, 6, 2, 0);   % recon L115
for mb = 1:4
    [kz, dkz] = calc_Pattern_SMS_blip(mb, 125, 0.063);              % imaging L114 (sliceSep = 3 mm * 42/2)
    F.(sprintf('kz_mb%d', mb))  = kz;
    F.(sprintf('dkz_mb%d', mb)) = dkz;
end
F.PEindex_18_49_3 = cal_PEindex(18, 49, 3);                          % calib L103
% tdr_sliceorder is encrypted P-code: record enough cases to pin down its rule.
for n = [1:12, 20, 21, 22, 41, 42, 43]
    F.(sprintf('sliceorder_mode2_N%d', n)) = tdr_sliceorder(n, 2);
end

%% 3) RF pulse design (TOPPE / John Pauly SLR), as in seq_EPTI_imaging_P1.m L57-84
F.dzrf_200_4_st_ls = toppe.utils.rf.jpauly.dzrf(200, 4, 'st', 'ls');
B0 = 2.89;
sys  = mr.opts('MaxGrad',64,'GradUnit','mT/m','MaxSlew',200,'SlewUnit','T/m/s', ...
    'rfRingdownTime',60e-6,'rfDeadtime',100e-6,'adcDeadTime',20e-6,'adcRasterTime',2.5e-6, ...
    'gradRasterTime',10e-6,'blockDurationRaster',10e-6,'B0',B0,'rfRasterTime',1e-6);
sys2 = mr.opts('MaxGrad',64,'GradUnit','mT/m','MaxSlew',160,'SlewUnit','T/m/s', ...
    'rfRingdownTime',60e-6,'rfDeadtime',100e-6,'adcDeadTime',20e-6,'adcRasterTime',2.5e-6, ...
    'gradRasterTime',10e-6,'blockDurationRaster',10e-6,'B0',B0,'rfRasterTime',1e-6);
sysGE = toppe.systemspecs('maxGrad', 0.5*sys.maxGrad/sys.gamma*100, ...
    'maxSlew', 0.4*sys.maxSlew/sys.gamma/10, 'maxRF', 0.25);
F.sysGE = sysGE;
% base (single-band) SLR pulse on the 4 us GE raster: 60 deg, 0.3 cm, TBW 4, 4 ms
[F.makeslr_rf_G, F.makeslr_gz_Gcm] = toppe.utils.rf.makeslr(60, 0.3, 4, 4, 1e-6, sysGE, ...
    'type', 'st', 'ftype', 'ls', 'writeModFile', false);
% multiband pulse converted to Pulseq units (MB = 2, slice separation 63 mm)
[rf, gzRF, freq] = getsmspulse(60, 3e-3, 4, 4e-3, 2, 3e-3*21, sysGE, sys2, ...
    'doSim', 0, 'type', 'st', 'ftype', 'ls');
F.sms_rf_signal  = rf.signal(:);
F.sms_rf_t       = rf.t(:);
F.sms_rf_delay   = rf.delay;
F.sms_rf_center  = mr.calcRfCenter(rf);
F.sms_gz_waveform = gzRF.waveform(:);
F.sms_gz_delay   = gzRF.delay;
F.sms_freq       = freq;

%% 4) P-code black boxes: EPI ghost correction, readout DFTs, apodization
% Synthetic navigator with a KNOWN odd/even phase error phi0 + phi1*x applied to the reversed (SEG=2) line.
% Layout follows data_phasor1d in recon_EPTI_P1.m: [Nx, 3 navs, coil, 1, 1, 1, pe_shots, SEG(2)].
rng(0);
Nx = 144; nch = 4; npe = 5;
x = ((0:Nx-1)' - Nx/2) / Nx;                          % readout position, normalized to [-0.5, 0.5)
obj = exp(-x.^2/0.02) .* (1 + 0.3*cos(6*pi*x));       % smooth 1D object
F.ghost_phi0 = 0.4;                                   % rad
F.ghost_phi1 = 3.0;                                   % rad per unit normalized x
ghost = exp(1i*(F.ghost_phi0 + F.ghost_phi1*x));
nav  = zeros(Nx, 3, nch, 1, 1, 1, npe, 2);
data = zeros(Nx, 6, nch, 1, 1, 1, npe, 2);            % a short echo train: odd echoes SEG 1, even echoes SEG 2
for c = 1:nch
    coil = exp(1i*0.5*c) * (1 + 0.2*c*x);
    for p = 1:npe
        s = obj .* coil * (1 + 0.05*p);
        k_fwd = fftc(s, 1);
        k_rev = fftc(s .* ghost, 1);
        nav(:,1,c,1,1,1,p,1) = k_fwd;
        nav(:,3,c,1,1,1,p,1) = k_fwd;
        nav(:,2,c,1,1,1,p,2) = k_rev;
        data(:,1:2:end,c,1,1,1,p,1) = repmat(k_fwd, [1 3]);
        data(:,2:2:end,c,1,1,1,p,2) = repmat(k_rev, [1 3]);
    end
end
F.ghost_nav_in  = nav;
F.ghost_data_in = data;
F.ghost_coeff_raw = mrir_artifact_ghost_compute(nav);                     % per coil and PE shot
params = struct('nslice', 1, 'SMS', 1, 'DestNx', Nx);
F.ghost_linear_fit_coeff = EPTI_SMS_Estimate_Linear_Coeff_WJG(params, nav);   % [2, coil, slice], recon L81
F.ghost_deghost_out = EPTI_SMS_Deghost_WJG(F.ghost_linear_fit_coeff, data);  % [Nx, echo, coil, pe, slice], recon L83

F.dft_in = randn(Nx, 3, 2) + 1i*randn(Nx, 3, 2);
F.idft_freqencode_out = mrir_iDFT_freqencode(F.dft_in);
F.fdft_freqencode_out = mrir_fDFT_freqencode(F.dft_in);
F.apod_ones72_dim1 = mrir_filter_raw_apodize_1d(ones(72, 72), 1, 0.1);   % SubspaceRecon_Process_v5.m L98
F.apod_ones72_dim2 = mrir_filter_raw_apodize_1d(ones(72, 72), 2, 0.1);

%% 5) Parameters read back from the .seq files, and the subspace basis (recon L60-63, L98)
p_img = EPTI_SMS_Ext_parameters_pulseq(fullfile(refdir, 'P1_EPTI_sms2.seq'), []);
p_cal = EPTI_SMS_Ext_parameters_pulseq(fullfile(refdir, 'P1_Calib.seq'), []);
F.params_imaging = p_img;
F.params_calib   = p_cal;
p_img.nechoGE = p_img.echo_num;   % set by the raw reader in the original pipeline
p_img.pf_echo = 0;                % recon L97
[F.Phi, F.TEs_all, F.S] = EPTI_gen_Subspace_Basis_GE(p_img, 4);

save(fullfile(refdir, 'fixtures.mat'), '-struct', 'F', '-v7');
fprintf('Wrote reference outputs to %s\n', refdir);
end


function run_script_isolated(folder, name)
% Run one of the MGH scripts from its own folder. The scripts begin with `clear`, which only clears this
% helper's workspace, so we can't keep local variables across the call; recompute the return folder afterwards.
cd(folder);
run(fullfile(folder, name));
close all;
cd(fileparts(mfilename('fullpath')));
end
