function [E, m, ky_idx] = generate_sampled_dft(mask_sample, do_check)
% GENERATE_SAMPLED_DFT  Generate the forward operator for sampled DFT
%
%  [E, m, ky_idx] = generate_sampled_dft(mask_sample)
%  [E, m, ky_idx] = generate_sampled_dft(mask_sample, do_check)
%
%  INPUT
%    mask_sample  Sampling mask. Size: [nx ny 1 nt 1 1 ng]
%    do_check     Runs sanity check. Optional. Default: true. 
%
%  OUTPUT
%    E            Sampled DFT matrix. Size: [1 ny 1 nt 1 1 ng 1 m]
%    m            Number of samples per (TE, SMS group). 
%    ky_idx       Sampled ky-indices. Size: [m nt ng].
%
%  EXAMPLE USAGE WITH BART
%    _E fmac <operand index> 2  % multiply by E and sum across BART dimension 1 (ky)
%    a E                        % apply the adjoint
%  
% Kaibo Tang, 2026 

if nargin < 2 || isempty(do_check)
    do_check = true;
end

% check mask_sample size
sz = size(mask_sample);
sz(end+1:8) = 1; 
ny = sz(2);  
nt = sz(4);  
ng = sz(7);
assert( ...
    sz(3)==1 && sz(5)==1 && sz(6)==1 && sz(8)==1, ...
    ['generate_sampled_dft: expect mask_sample to have singleton dimensions at dimensions 3, 5, 6, 8; got size %s.'], mat2str(sz) ...
);
mask_sample = (mask_sample ~= 0);
assert( ...
    isequal(any(mask_sample,1), all(mask_sample,1)), ...
    'generate_sampled_dft: expect mask_sample to be constant along the readout dimension.' ...
);
mk = reshape(mask_sample(1,:,1,:,1,1,:,1), [ny nt ng]);  % [ny nt ng]

% centered unitary IDFT/DFT matrices
Finv = ifftc(eye(ny),1);  % IDFT
Fc = Finv';  % DFT

% count the number of ky-lines per (TE, SMS group)
cnt = reshape(sum(mk,1), [nt ng]);  % [nt ng]
m = max(cnt(:));

% generate E matrix
E_mat  = complex(zeros(m, ny, nt, ng));  % [m ny nt ng]
ky_idx = zeros(m, nt, ng);  % [m nt ng]
for g = 1:ng
    for t = 1:nt
        idx = find(mk(:,t,g));  % find all sampled ky-indices
        for j = 1:numel(idx)
            E_mat(j,:,t,g) = Fc(idx(j),:);
            ky_idx(j,t,g)  = idx(j);
        end
    end
end

% run sanity check
if do_check
    u = randn(ny,1) + 1j*randn(ny,1);
    for g = 1:ng
        for t = 1:nt
            Et  = E_mat(:,:,t,g);
            ref = Finv * (double(mk(:,t,g)) .* (Fc*u));  % use full DFT/IDFT matrices
            err = norm(Et'*(Et*u) - ref) / max(norm(ref), eps);  % use E and E'
            assert( ...
                err < 1e-9, ...
                'generate_sampled_dft: sanity check failed at TE %d, group %d (rel. err. %g).', t, g, err ...
            );
        end
    end
end

% reshape E matrix for BART
E = reshape(permute(E_mat,[2 3 4 1]), [1 ny 1 nt 1 1 ng 1 m]);

end
