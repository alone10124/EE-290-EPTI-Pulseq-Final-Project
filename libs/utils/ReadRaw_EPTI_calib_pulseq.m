function [data_temp,data_phasor1d_temp,parameters] = ReadRaw_EPTI_calib_pulseq(filename_calib,parameters,ISMRMRD_data)

if (ISMRMRD_data==1)
    [data_temp,data_phasor1d_temp,parameters] = ReadRaw_EPTI_calib_pulseq_ISMRMRD(filename_calib,parameters);
else
    [data_temp,data_phasor1d_temp,parameters] = ReadRaw_EPTI_calib_pulseq_Siemens(filename_calib,parameters);
end