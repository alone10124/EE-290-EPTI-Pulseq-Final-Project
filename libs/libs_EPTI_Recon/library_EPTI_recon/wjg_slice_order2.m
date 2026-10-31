function [slice_order,position,readDir] = wjg_slice_order2(twix_obj)
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


%% slice direction
readDir = zeros(3,nslice/SMS);
position = zeros(3,nslice/SMS);
hdr = twix_obj.hdr;
sl = hdr.MeasYaps.sSliceArray.asSlice;

for loopi = 1:nslice/SMS
    s = sl{loopi};

    if~isfield(s.sPosition, 'dSag')
        s.sPosition.dSag=0;
    end
    if~isfield(s.sPosition, 'dCor')
        s.sPosition.dCor=0;
    end
    if~isfield(s.sPosition, 'dTra')
        s.sPosition.dTra=0;
    end


    if~isfield(s.sNormal, 'dSag')
        s.sNormal.dSag=0;
    end
    if~isfield(s.sNormal, 'dCor')
        s.sNormal.dCor=0;
    end
    if~isfield(s.sNormal, 'dTra')
        s.sNormal.dTra=0;
    end
    if~isfield(s, 'dInPlaneRot')
        s.dInPlaneRot=0;
    end


    
    n = [s.sNormal.dSag;
         s.sNormal.dCor;
         s.sNormal.dTra];
    
    n = n / norm(n);       % slice normal in Siemens Sag/Cor/Tra coordinates
    
    pos = [s.sPosition.dSag;
           s.sPosition.dCor;
           s.sPosition.dTra];   % slice center / position, typically mm
    
    inPlaneRot = s.dInPlaneRot; % radians
    
    
    ref = [0; 0; 1];                  % choose a vector not parallel to n
    if abs(dot(ref,n)) > 0.95
        ref = [0; 1; 0];
    end
    
    % Arbitrary orthonormal in-plane basis before applying Siemens rotation
    e1 = cross(n, ref);
    e1 = e1 / norm(e1);
    
    e2 = cross(n, e1);
    e2 = e2 / norm(e2);
    
    % Apply the in-plane rotation
    R2 = [cos(inPlaneRot), -sin(inPlaneRot);
          sin(inPlaneRot),  cos(inPlaneRot)];
    
    inplane = [e1 e2] * R2;
    
    readDir_temp  = inplane(:,1);
    phaseDir_temp = inplane(:,2);
    position(:,loopi) = [s.sPosition.dSag,s.sPosition.dCor,s.sPosition.dTra];
    readDir(:,loopi) = readDir_temp;

end