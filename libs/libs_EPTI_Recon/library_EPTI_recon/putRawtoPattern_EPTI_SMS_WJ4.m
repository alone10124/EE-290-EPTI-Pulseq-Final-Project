function [kdata] = putRawtoPattern_EPTI_SMS_WJ4(kdata,parameters,isge)
% put raw k-space EPTI data into EPTI encoding pattern
% change the Ngroup to be the number of echoes instead of echo group
% to allow arbitrary echo # for subspace reconstruction
% Fuyixue Wang 2018

% support EPTI with SMS
% Fuyixue Wang 2019
% Jian Wu 2025
MB = parameters.SMS;
[N_total,nx,~,nc,nslice]=size(kdata);

if (isge)
    ky_acq = parameters.sampling_pattern_GE;
else
    ky_acq = parameters.sampling_pattern_SE;
end

npe=parameters.Rseg*parameters.Nseg;
kdata_EPTI=zeros(N_total,nx,npe,nc,nslice);
for shot=1:parameters.Nseg
    for nt=1:N_total
        kdata_EPTI(nt,:,ky_acq(shot,nt),:,:)=kdata(nt,:,shot,:,:);
    end
end

% if (isinverse)
%     for shot=1:Nseg
%         for nt=1:N_total
%             kdata_EPTI(nt,:,ky_acq(shot,nt),:,:)=kdata(nt,:,Nseg+1-shot,:,:);
%         end
%     end
% else
%     for shot=1:Nseg
%         for nt=1:N_total
%             kdata_EPTI(nt,:,ky_acq(shot,nt),:,:)=kdata(nt,:,shot,:,:);
%         end
%     end
% end
% kdata_EPTI=kdata_EPTI(:,:,end:-1:1,:,:);

if MB>1
    kdata=zeros(N_total,nx,npe,MB,nc,nslice);
    for i = 1:MB
        kdata(i:MB:end,:,:,i,:,:)=kdata_EPTI(i:MB:end,:,:,:,:);
    end
else
    kdata=reshape(kdata_EPTI,[N_total,nx,npe,1,nc,nslice]);
end

end


