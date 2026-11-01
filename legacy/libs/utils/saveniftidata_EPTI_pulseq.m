function saveniftidata_EPTI_pulseq(data,parameters,fndata)
% save NIFTI data with correct header info read from rawdata
% data size [nx,ny,nsl,ndyns]

% by Fuyixue Wang based on Jon's mrir
% v2, output NIFTI without mri convert from MGZ
    meas = parameters.nifti_header;
    disp(sprintf('==> [%s]: saving file "%s"...', mfilename, fndata));
    % meas.prot.tPatientPosition = 'HFS';

    meas.prot.tPatientPosition = 'HFS';

    M0_vox2ras = mrir_measdat_vox2ras_EPTI(meas.prot, meas.evp);

    res = parameters.resolution*1000;

    save_nii(make_nii_EPTI(double(data), [res',parameters.TotalDuration], [], 16,M0_vox2ras), fndata);

end