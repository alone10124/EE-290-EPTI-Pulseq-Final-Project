function ge_raw_export(rawfile, outfile, varargin)
% ge_raw_export - Turn GE raw data into a plain HDF5 file in ADC order, for the Python recon.
%
%   tv6 (P-file) -- needs the _tv6.tar that was actually run, to map P-file slots back to ADC order:
%     matlab -batch "cd tools/ge; ge_raw_export('P12345.7', '../../out/scan/epti_3mm_ge_calib.h5', 'tar', '../../sequences/epti_3mm_ge_calib_tv6.tar')"
%   pge2 (ScanArchive) -- needs GE's Orchestra SDK for MATLAB on the path (ask Hannah for it):
%     matlab -batch "addpath('C:/orchestra-sdk/matlab'); cd tools/ge; ge_raw_export('ScanArchive_xxx.h5', '../../out/scan/epti_3mm_ge_calib.h5')"
%
% Output: HDF5 with datasets /adc_real and /adc_imag, MATLAB size [nSamples, nCoils, nADC] (single). In Python
% (row-major) that reads as [nADC, nCoils, nSamples], the layout epti.io.ismrmrd_io.read_acquisitions returns, in
% the order the ADC events occur in the .seq file. epti/io/ge_raw.py reads it.
%
% Why MATLAB: the P-file header layout changes with every GE software release, and toppe.utils.loadpfile (maintained
% by the TOPPE authors) handles that; ScanArchive files can only be read with GE's Orchestra SDK.
%
% Options:
%   'tar'         tv6 .tar file used for this scan (required for P-files)
%   'pfileData'   [nSamples nCoils nSlices nEchoes nViews] array to use instead of reading rawfile (for tests)

arg.tar = '';
arg.pfileData = [];
arg = vararg_pair_local(arg, varargin);

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'toppe'));                                          % +toppe (loadpfile, readloop)
addpath(fullfile(here, 'PulCeq', 'matlab'));                               % +pge2 (loaddata), textprogressbar

[~, ~, ext] = fileparts(rawfile);
isScanArchive = strcmpi(ext, '.h5') && isempty(arg.pfileData);

if isScanArchive
    % pge2: frames come out in acquisition order, one per ADC event, each [nSamples nCoils].
    if exist('GERecon', 'file') == 0
        error('ScanArchive needs GE''s Orchestra SDK (GERecon) on the MATLAB path');
    end
    frames = pge2.utils.loaddata(rawfile);
    d = cat(3, frames{:});                                                  % [nSamples nCoils nADC]
    source = 'pge2 ScanArchive';
else
    % tv6: the P-file stores each ADC in a (dabslice, dabecho, dabview) slot set by scanloop.txt
    if isempty(arg.tar)
        error('For a P-file, pass the tv6 .tar that was run: ge_raw_export(pfile, out, ''tar'', ''x_tv6.tar'')');
    end
    if isempty(arg.pfileData)
        % slices 2..N: loadpfile skips dabslice 0 (baseline slot), so output slice index k == dabslice k
        dat = toppe.utils.loadpfile(rawfile, [], [], [], 'quiet', true, 'returnAsDouble', false);
    else
        dat = arg.pfileData;
    end
    [hasADC, loop] = read_tar(arg.tar);
    m = loop(:, 1);                                     % module index per row, 0 = pure delay
    rows = loop(m > 0, :);
    rows = rows(hasADC(rows(:, 1)) == 1, :);
    sl = rows(:, 7); ec = rows(:, 8) + 1; vw = rows(:, 9);   % dabslice (1-based), dabecho (0-based), dabview (1-based)
    sz = size(dat, 1:5);
    if any(sl > sz(3)) || any(ec > sz(4)) || any(vw > sz(5))
        error(['P-file is smaller than scanloop.txt expects (P-file [slices echoes views] = [%d %d %d], ' ...
            'scanloop max [%d %d %d]). Is this the right .tar?'], sz(3:5), max(sl), max(ec), max(vw));
    end
    d = zeros(sz(1), sz(2), numel(sl), 'like', single(1i));
    for n = 1:numel(sl)
        d(:, :, n) = single(dat(:, :, sl(n), ec(n), vw(n)));
    end
    source = 'tv6 P-file';
end

fprintf('%s: %d ADC events, %d samples, %d coils\n', source, size(d, 3), size(d, 1), size(d, 2));
outdir = fileparts(outfile);
if ~isempty(outdir) && ~exist(outdir, 'dir'), mkdir(outdir); end
if exist(outfile, 'file'), delete(outfile); end
h5create(outfile, '/adc_real', size(d), 'Datatype', 'single');
h5create(outfile, '/adc_imag', size(d), 'Datatype', 'single');
h5write(outfile, '/adc_real', real(d));
h5write(outfile, '/adc_imag', imag(d));
h5writeatt(outfile, '/', 'source', source);
h5writeatt(outfile, '/', 'rawfile', rawfile);
fprintf('wrote %s\n', outfile);
end


function [hasADC, loop] = read_tar(tarfile)
% Per-module ADC flag (modules.txt) and the scan loop (scanloop.txt) from a tv6 .tar
tmp = tempname; mkdir(tmp);
c = onCleanup(@() rmdir(tmp, 's'));
untar(tarfile, tmp);
% modules.txt: 3 header lines, then "fname dur hasRF hasADC trigpos" per module
lines = splitlines(strtrim(fileread(fullfile(tmp, 'modules.txt'))));
hasADC = zeros(numel(lines) - 3, 1);
for k = 4:numel(lines)
    f = strsplit(strtrim(lines{k}));
    hasADC(k - 3) = str2double(f{4});
end
loop = toppe.readloop(fullfile(tmp, 'scanloop.txt'));
end


function arg = vararg_pair_local(arg, opts)
for k = 1:2:numel(opts)
    if ~isfield(arg, opts{k}), error('unknown option %s', opts{k}); end
    arg.(opts{k}) = opts{k+1};
end
end
