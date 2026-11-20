function test_raw_export(adcfile, tarfile, outfile, maxView)
% test_raw_export - Fake a tv6 P-file from ADC-ordered data and run it through ge_raw_export (used by
% tests/test_ge_raw.py). Run it through that test; it isn't meant to be called directly.
%
% adcfile holds /adc_real, /adc_imag ([nSamples nCoils nADC] in MATLAB) in ADC order. Each ADC is put in the P-file
% slot the tv6 converter assigns it (PulCeq tv6 ceq2ge.m: view = mod(n, maxView) + 1, slice = floor(n/maxView) + 1),
% which is independent of the scanloop.txt parsing that ge_raw_export uses to undo it.
d = complex(h5read(adcfile, '/adc_real'), h5read(adcfile, '/adc_imag'));
[ns, nc, nadc] = size(d);
n = 0:nadc-1;
view = mod(n, maxView) + 1;
sl = floor(n / maxView) + 1;
p = zeros(ns, nc, max(sl) + 1, 1, maxView, 'single');   % one spare slice, as the scanner allocates extra slots
p = complex(p, p);
for k = 1:nadc
    p(:, :, sl(k), 1, view(k)) = d(:, :, k);
end
ge_raw_export('fake.7', outfile, 'tar', tarfile, 'pfileData', p);
end
