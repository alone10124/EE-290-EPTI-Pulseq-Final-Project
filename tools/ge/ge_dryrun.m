function ge_dryrun(seqfiles, varargin)
% ge_dryrun - Convert .seq files to GE pge2 format and run PulCeq's hardware checks, before going to the scanner.
%
%   matlab -batch "cd tools/ge; ge_dryrun"                         % both GE EPTI scans
%   matlab -batch "cd tools/ge; ge_dryrun({'../../sequences/x.seq'})"
%
% For each .seq file: seq2ceq (the same conversion done at the scanner), pge2.check (block timing on GE rasters,
% peak B1 / gradient / slew per segment instance, GE PNS per segment), then writeceq -> <name>.pge next to the .seq.
% Hardware values are the GE MR750w / XRMW ones (pge2.opts table). Change them with name/value options once the
% scanner facility confirms them, e.g. ge_dryrun([], 'coil', 'xrm', 'g_max', 5, 'slew_max', 20).

here = fileparts(mfilename('fullpath'));
repo = fullfile(here, '..', '..');
addpath(fullfile(repo, 'legacy', 'libs', 'pulseq-master151', 'matlab'));   % Pulseq 1.5.1 (mr.*)
addpath(fullfile(here, 'PulCeq', 'matlab'), fullfile(here, 'PulCeq', 'matlab', 'DataHash'));                            % seq2ceq, writeceq, +pge2

if nargin < 1 || isempty(seqfiles)
    seqfiles = {fullfile(repo, 'sequences', 'epti_3mm_ge_calib.seq'), ...
                fullfile(repo, 'sequences', 'epti_3mm_ge_imaging.seq')};
end

arg.psd_rf_wait = 100e-6;   % s; scanner-specific RF/gradient delay, measured on site (Hannah sets this)
arg.psd_grd_wait = 100e-6;  % s; scanner-specific ADC/gradient delay
arg.b1_max = 0.25;          % G
arg.g_max = 3.3;            % G/cm   (MR750w XRMW: 33 mT/m)
arg.slew_max = 12.0;        % G/cm/ms (MR750w XRMW: 120 T/m/s)
arg.coil = 'xrmw';
arg.write = true;
arg.validate_row = 'all';   % pge2.validate: 'all' segments, or start at this .seq row
arg.validate = true;        % run pge2.validate (slow for the calibration scan)
arg = vararg_pair(arg, varargin);

sysGE = pge2.opts(arg.psd_rf_wait, arg.psd_grd_wait, arg.b1_max, arg.g_max, arg.slew_max, arg.coil);

allok = true;
for k = 1:numel(seqfiles)
    f = seqfiles{k};
    fprintf('\n===== %s\n', f);
    try
        seq = mr.Sequence();
        seq.read(f);
        ceq = seq2ceq(seq);
        fprintf('segments: %d, parent blocks: %d, rows: %d, readouts: %d, duration %.2f s\n', ...
            ceq.nSegments, ceq.nParentBlocks, ceq.nMax, ceq.nReadouts, ceq.duration);
        for i = 1:ceq.nSegments
            fprintf('  segment %d (TRID %d): %d blocks\n', i, ceq.segments(i).TRID, ceq.segments(i).nBlocksInSegment);
        end
        params = pge2.check(ceq, sysGE);
        fprintf('PASS  b1max %.4f G, gmax %.2f G/cm, smax %.2f G/cm/ms\n', params.b1max, params.gmax, params.smax);
        % Converted (Ceq) waveforms vs the original .seq: catches sign flips, conjugated RF and timing offsets.
        % Same function Hannah can run later against the scanner's scan.xml.* files (WTools / Pulse View output).
        if arg.validate
            if ~pge2.validate(ceq, sysGE, seq, [], 'row', arg.validate_row)
                error('pge2.validate: Ceq waveforms do not match the .seq file');
            end
            fprintf('PASS  pge2.validate (Ceq vs .seq waveforms)\n');
        end
        if arg.write
            [p, n] = fileparts(f);
            out = fullfile(p, [n '.pge']);
            pge2.writeceq(ceq, out, 'params', params);
            fprintf('wrote %s\n', out);
        end
    catch ME
        allok = false;
        fprintf(2, 'FAIL  %s\n', ME.message);
        for s = ME.stack(1:min(4, end))'
            fprintf(2, '    at %s:%d\n', s.name, s.line);
        end
    end
end
if ~allok
    error('ge_dryrun: at least one sequence failed');
end
end
