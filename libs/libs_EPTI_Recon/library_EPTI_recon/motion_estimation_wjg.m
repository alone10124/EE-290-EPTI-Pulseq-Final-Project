function [motion_file,finalI,motion_matrix] = motion_estimation_wjg(lowresI1,lowresI2,BART_savepath)
%estimation the motion parameters
% warning("off")
% lowresI1: the reference image
% lowresI2: the moving image
folder = BART_savepath;
ref = abs(lowresI1);
ref = ref./max(ref(:));
dir_nii = [folder,'tmp'];
mkdir(dir_nii);
mkdir([dir_nii,'/reg']);
save_nii(make_nii(double(squeeze(abs(ref))), [1,1,1], [], 16), [folder,'tmp/ref.nii']);
fdatabase = fullfile(dir_nii,'ref');
system(['rm -f ',fdatabase,'+* ']);
system(['3dcopy ',fdatabase,'.nii ',fdatabase]);

moving = abs(lowresI2);
moving = moving./max(moving(:));
save_nii(make_nii(double(squeeze(abs(moving))), [1,1,1], [], 16), [folder,'tmp/moving.nii']);
fdata = fullfile(dir_nii,'moving');
system(['rm -f ',fdata,'+* ']);
system(['3dcopy ',fdata,'.nii ',fdata]);
foutput = fullfile(dir_nii,'registered');
foutput_mcparam = fullfile(dir_nii,'param');  
system(['rm -f ',foutput,'+* ']);
[status, result] = system(['3dvolreg -verbose -base ',fdatabase,'+orig',...
                     ' -prefix ',foutput...
                     ' -dfile ',foutput_mcparam,'_Param_mc',...
                     ' -1Dfile ',foutput_mcparam,'_Param_mc.1D',...
                     ' -1Dmatrix_save ',foutput_mcparam,'_Mat_mc.mat.aff12.1D',...
                     ' -float -Fourier ',...
                     fdata,'+orig']);
out_file = [folder,'1'];
system(['rm -f ',out_file,'.nii']);
system(['3dcopy ',foutput,'+orig ',out_file,'.nii']);
%get the motion parameters
lines = strsplit(result, newline);
for i = 1:length(lines)
    if contains(lines{i}, '++ Mean:')
        meanLine = strtrim(lines{i});
        break;
    end
end
meanLine = meanLine(10:end);
meanLineReplaced = regexprep(meanLine, '\s{2,}', ';');
eval(meanLineReplaced);
motion_matrix = [roll;pitch;yaw;dS;dL;dP];
% package the data 
files = [out_file,'.nii'];
finalI = niftiread(files);

motion_file = [foutput_mcparam,'_Mat_mc.mat.aff12.1D'];


% figure(1);imshow3(abs(moving));title('before register')
% figure(2);imshow3(abs(finalI));title('after register')
% figure(3);imshow3(abs(ref));title('reference')
%% apply  the transform
% foutput = [folder,'after_register.nii'];
% [status, result] = system(['3dAllineate -source ',fdata,'+orig',...
%                      ' -1Dmatrix_apply ',foutput_mcparam,'_Mat_mc.mat.aff12.1D'...
%                      ' -prefix ',foutput])
%  data = load_nii(foutput);
%  data = data.img;
%  figure(4);imshow3(data);title('my register')