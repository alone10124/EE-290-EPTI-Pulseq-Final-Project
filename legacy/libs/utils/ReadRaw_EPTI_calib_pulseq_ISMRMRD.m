function [data_temp,data_phasor1d_temp,parameters] = ReadRaw_EPTI_calib_pulseq_ISMRMRD(filename_calib,parameters)

dset = ismrmrd.Dataset(filename_calib, 'dataset');
D = dset.readAcquisition();
hdr = ismrmrd.xml.deserialize(dset.readxml);
%% Encoding and reconstruction information
% Matrix size
Nx = hdr.encoding.encodedSpace.matrixSize.x;
try 
    nchannel = hdr.acquisitionSystemInformation.receiverChannels;
catch
    nchannel = 1;
end
%% Read all the data
echo_num = parameters.ETL;
dims  =[Nx,nchannel,echo_num,parameters.N_slices,parameters.pe_lines];
data = D.data;
data_unsorted = zeros(Nx,nchannel,echo_num*parameters.N_slices*parameters.pe_lines);
for loopi = 1:length(data)
    temp = data{loopi};
    data_unsorted(:,:,loopi) = temp;    
end
data = reshape(data_unsorted,dims);
data_phasor1d = data(:,:,1:3,:,:);
data_phasor1d = squeeze(data_phasor1d);
data_phasor1d = permute(data_phasor1d,[1,3,2,5,4]); 
dims = [Nx,3,nchannel,1,1,1,parameters.pe_lines,2,1,parameters.N_slices];
data_phasor1d_temp = zeros(dims);
data_phasor1d_temp(:,1:2:end,:,:,:,:,:,1,:,:) = squeeze(data_phasor1d(:,1:2:3,:,:,:));
data_phasor1d_temp(:,2:2:end,:,:,:,:,:,2,:,:) = data_phasor1d(end:-1:1,2,:,:,:);

data = data(:,:,4:end,:,:);
data = squeeze(data);
data = squeeze(permute(data,[1,3,2,5,4])); 
dims = [Nx,echo_num-3,nchannel,1,1,1,parameters.pe_lines,2,1,parameters.N_slices];
data_temp = zeros(dims);
data_temp(:,1:2:end,:,:,:,:,:,1,:,:) = squeeze(data(:,1:2:end,:,:,:));
data_temp(:,2:2:end,:,:,:,:,:,2,:,:) = squeeze(data(end:-1:1,2:2:end,:,:,:));

%update the parameters
parameters.readout = Nx;
parameters.nchannel = nchannel;
parameters.Nseg = 1;
parameters.iEffectiveEpiEchoSpacing =  (parameters.rampUpTime+parameters.rampDownTime+parameters.flatTopTime)*1e6;

