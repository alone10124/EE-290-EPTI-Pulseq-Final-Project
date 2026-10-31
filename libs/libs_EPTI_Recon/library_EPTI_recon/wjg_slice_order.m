function slice_order = wjg_slice_order(twix_obj)
%get slice order from the raw data
nslice = twix_obj.hdr.Config.NSlc;
SMS = twix_obj.hdr.Meas.MultiBandFactor;
rep_num = twix_obj.image.NRep;
Nseg = cell2mat(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree(34));
temp1 = twix_obj.image.slicePos;
temp1 = squeeze(sum(temp1(:,1:length(temp1)/rep_num/Nseg),1));
slice_order = combine_consecutive_values(temp1);
slice_order = slice_order(1:nslice/SMS);
[~,slice_order] = sort(slice_order);