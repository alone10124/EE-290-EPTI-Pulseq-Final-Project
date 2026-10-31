function warpI = warp_wjg(originalI,motion_param_file,BART_savepath)
% apply the motion parameters on  the originalI
% do not suppor the complex

folder = BART_savepath;
dir_nii = [folder,'tmp'];
mkdir(dir_nii);
mkdir([dir_nii,'/reg']);
moving = originalI;
save_nii(make_nii(double(squeeze((moving))), [1,1,1], [], 16), [folder,'tmp/moving.nii']);
fdata = fullfile(dir_nii,'moving.nii');

% apply  the transform
foutput = [dir_nii,'/after_register.nii'];
system(['rm -f ',foutput]);

[status, result] = system(['3dAllineate -source ',fdata,...
                     ' -1Dmatrix_apply ',motion_param_file...
                     ' -prefix ',foutput]);
 warpI = niftiread(foutput);
 % figure(1);imshow3(originalI,[-100,100]);title('original');colorbar;colormap jet
 % figure(2);imshow3(warpI,[-100,100]);title('my register');colorbar;colormap jet