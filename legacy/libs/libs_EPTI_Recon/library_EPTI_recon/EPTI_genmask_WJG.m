function msk = EPTI_genmask_WJG(im_GRE_comb)
% EPTI generate the mask
% By Jian Wu MGH 2025
% Fuyixue Wang & Zijing Dong, MGH, 2023
msk = abs(squeeze(im_GRE_comb(:,:,:,1)));
temp = msk;
maxinten =max(temp(:));
temp = temp./max(temp(:));
threshold = graythresh(temp)/3;
for loopi = 1:size(msk,3)
    temp = msk(:,:,loopi);
    temp = temp./maxinten;
    temp = imbinarize(temp,threshold);
    temp = bwareaopen(temp,10);
    % se = strel('disk', 4);
    % temp = imdilate(temp, se);
    temp=imfill(temp,'holes');
    msk(:,:,loopi) = temp;
end
