function [EPTIdata,data_phase,parameters] = ReadRaw_EPTI_imaging_pulseq(filename_EPTI,parameters,ISMRMRD_data)

if (ISMRMRD_data==1)
    [EPTIdata,data_phase,parameters] = ReadRaw_EPTI_pulseq_ISMRMRD(filename_EPTI,parameters);
else
    [twix_obj,EPTIdata,data_phase,parameters] = ReadRaw_EPTI_pulseq_Siemens(filename_EPTI,parameters);
end