function [EPTIdata,data_phase,parameters] = ReadRaw_EPTI_pulseq_ISMRMRD(filename_EPTI,parameters)

dset = ismrmrd.Dataset(filename_EPTI, 'dataset');
D = dset.readAcquisition();
hdr = ismrmrd.xml.deserialize(dset.readxml);
% Matrix size
Nx = hdr.encoding.encodedSpace.matrixSize.x;
try 
    nchannel = hdr.acquisitionSystemInformation.receiverChannels;
catch
    nchannel = 1;
end
echo_num = parameters.echo_num;
slice = parameters.nslice/parameters.SMS;
data = D.data;
rep_num = size(data,2)/slice/echo_num;
dims = [Nx,nchannel,echo_num,slice,rep_num];
data_unsorted = zeros(dims);
for loopi = 1:length(data)
    temp = data{loopi};
    data_unsorted(:,:,loopi) = temp;    
end
EPTIdata = reshape(data_unsorted,dims);
data_phase = [];

%update the parameters
parameters.readout = Nx;
parameters.nchannel = nchannel;
parameters.lPhaseEncodingLines = echo_num;
parameters.nechoGE = echo_num;
parameters.iEffectiveEpiEchoSpacing =  (parameters.rampUpTime+parameters.rampDownTime+parameters.flatTopTime)*1e6;
parameters.NslicesEX = parameters.SMS;