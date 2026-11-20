function ge_dryrun_tv6(seqfiles, maxView)
% ge_dryrun_tv6 - Convert .seq files with the older TOPPE v6 path (PulCeq branch tv6_v1.10: seq2ge -> .tar).
%
%   matlab -batch "cd tools/ge; ge_dryrun_tv6"
%
% Use this if the scanner runs the tv6 TOPPE interpreter rather than pge2. Run it in its own MATLAB session:
% the tv6 branch has its own seq2ceq.m that must not be mixed with the pge2 one (ge_dryrun.m).
% Writes <name>_tv6.tar next to each .seq file. maxView (default 2048) is how many ADCs go in one P-file 'slice' slot
% before the next slot starts; ge_raw_export.m reads the mapping back from the .tar, so any value works.

here = fileparts(mfilename('fullpath'));
repo = fullfile(here, '..', '..');
addpath(fullfile(repo, 'legacy', 'libs', 'pulseq-master151', 'matlab'));
addpath(fullfile(here, 'toppe'));                 % +toppe (current main)
addpath(fullfile(here, 'PulCeq_tv6', 'matlab'));  % seq2ge, ceq2ge (git worktree of branch tv6_v1.10)
addpath(fullfile(here, 'tv6_patch'));             % writeentryfile/writecoresfile: v6 entry format + LF line endings (first on path)

if nargin < 1 || isempty(seqfiles)
    seqfiles = {fullfile(repo, 'sequences', 'epti_3mm_ge_calib.seq'), ...
                fullfile(repo, 'sequences', 'epti_3mm_ge_imaging.seq')};
end

if nargin < 2 || isempty(maxView)
    maxView = 2048;
end

% MR750w / XRMW limits; psd_*_wait are site-specific (Hannah's values go here).
sysGE = toppe.systemspecs('maxGrad', 3.3, 'maxSlew', 12.0, 'gradient', 'xrmw', ...
    'psd_rf_wait', 100, 'psd_grd_wait', 100, ...
    'rfDeadTime', 100, 'rfRingdownTime', 60, 'adcDeadTime', 40, 'maxView', maxView);

allok = true;
for k = 1:numel(seqfiles)
    f = char(java.io.File(seqfiles{k}).getCanonicalPath());   % absolute: we cd into a temp folder below
    [p, n] = fileparts(f);
    out = fullfile(p, [n '_tv6.tar']);
    fprintf('\n===== %s\n', f);
    % ceq2ge writes its scan files into the current folder and tars them with a shell 'tar', which (Git's tar on
    % Windows) can't take a 'C:\...' path. So work in a temp folder with a relative tar name, then move it.
    work = fullfile(tempdir, ['tv6_' n]);
    if exist(work, 'dir'), rmdir(work, 's'); end
    mkdir(work);
    old = cd(work);
    try
        seq2ge(f, sysGE, 'scan.tar');
        cd(old);
        movefile(fullfile(work, 'scan.tar'), out, 'f');
        fprintf('PASS  wrote %s\n', out);
    catch ME
        cd(old);
        allok = false;
        fprintf(2, 'FAIL  %s\n', ME.message);
        for s = ME.stack(1:min(6, end))'
            fprintf(2, '    at %s:%d\n', s.name, s.line);
        end
    end
end
if ~allok
    error('ge_dryrun_tv6: at least one sequence failed');
end
end
