function [kdata] = putRawtoPattern_SMS_zigzag_PF_selfnav_WJ4(kdata,parameters,startpe,isge)
% put raw k-space EPTI data into EPTI encoding pattern
% change the Ngroup to be the number of echoes instead of echo group
% to allow arbitrary echo # for subspace reconstruction
% Fuyixue Wang 2018

% support EPTI with SMS
% Fuyixue Wang 2019

% put into pattern of zigzag single-shot EPTI data
% Fuyixue Wang 2021, MGH

% works for partial Fourier acquisition
% Zijing Dong 2022, MGH
% Jian WU 2025
[nt,nx,~,nc,nslice]=size(kdata);
MB = parameters.SMS;
npe = parameters.Nseg*parameters.Rseg;
if (isge)
    ky_acq = parameters.sampling_pattern_GE-parameters.PF_shift;
else
    ky_acq = parameters.sampling_pattern_SE-parameters.PF_shift;
end
%%  put data into pattern
kdata_EPTI=zeros(nt,nx,npe,nc,nslice);

for shot=1:parameters.Nseg
    for t=1:nt
        kdata_EPTI(t,:,ky_acq(shot,t),:,:)=kdata(t,:,shot,:,:);
    end
end
if (MB==3)
    kdata=zeros(t,nx,npe,2,nc,nslice);
    kdata(1:startpe-1,:,:,1,:,:)=kdata_EPTI(1:startpe-1,:,:,:,:);
    for i = 1:2
        kdata(startpe-1+i:2:end,:,:,i,:,:)=kdata_EPTI(startpe-1+i:2:end,:,:,:,:);
    end
elseif (MB==2)
    kdata=zeros(t,nx,npe,2,nc,nslice);
    kdata(1:startpe-1,:,:,1,:,:)=kdata_EPTI(1:startpe-1,:,:,:,:);
    for i = 1:MB
        kdata(startpe-1+i:2:end,:,:,i,:,:)=kdata_EPTI(startpe-1+i:2:end,:,:,:,:);
    end
elseif (MB==4)
    kdata=zeros(t,nx,npe,3,nc,nslice);
    kdata(1:startpe-1,:,:,1,:,:)=kdata_EPTI(1:startpe-1,:,:,:,:);
    for i = 1:3
        kdata(startpe-1+i:3:end,:,:,i,:,:)=kdata_EPTI(startpe-1+i:3:end,:,:,:,:);
    end
else
    kdata=reshape(kdata_EPTI,[nt,nx,npe,1,nc,nslice]);
end

end



