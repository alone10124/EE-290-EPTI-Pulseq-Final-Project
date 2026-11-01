function status = gradSpectrum_bb_v1(seq,sys,sys_type, fmax2use,esp)
% WJG
% pulseq\matlab\demoUnsorted\gradSpectrum.m
% 
% example script to plot gradient frequency spectrum
% expects: 
%  "seq" object to be already populated 
%  "sys" object to contain system specs
% furthermore, on Siemens you need the *.asc file for your gradient system

status = 1;
dt=sys.gradRasterTime; % time raster
fmax=10000; %10kHz
nwin=5000; % 0.05s
os=3; % frequency oversampling for prettier peaks

if strcmp(sys_type,'prisma')
        ascName = 'MP_GPA_K2309_2250V_951A_AS82.asc'; % prisma
else
    if strcmp(sys_type,'skyra')
        ascName = 'MP_GPA_K2309_2250V_793A_GC99.asc'; % skyra
    else
        if strcmp(sys_type,'Connectome2') || strcmp(sys_type,'C2_simulate_prisma')
            ascName = 'MP_GradSys_P034_c_CX600.asc';    
        else
            if strcmp(sys_type,'terra')
                ascName = 'MP_GPA_K2298_2250V_793A_SC72CD.asc'; 
            else
                warning('Unknown system :%s\nSetting asc filename to []', sys_type)
                ascName = [];
            end
        end
    end
end

% ascName=[]; % this disables the display of the system's resonance frequences
%ascName='idea/asc/MP_GPA_K2309_2250V_951A_AS82.asc'; % 3T prisma
% ascName='idea/asc/MP_GradSys_P034_X60.asc'; % 3T cima.X

if ischar(ascName)
    ascData=mr.Siemens.readasc(ascName);
end

faxis=(0:(nwin/2-1))/nwin/dt/os;

if nargin < 4
    nfmax=sum(faxis<=fmax);
else
    nfmax=sum(faxis<=fmax2use);
end

wave_data=seq.waveforms_and_times();
tmax=max([wave_data{1}(1,end) wave_data{2}(1,end) wave_data{3}(1,end)]);
nt=ceil(tmax/dt);
tmax=nt*dt;

gw=zeros(3,nt);
for i=1:3
    gw(i,:)=interp1(wave_data{i}(1,:),wave_data{i}(2,:),((1:nt)-0.5)*dt,'linear',0);
end

gs=[];

ng=size(gw,1);
for g=1:ng
    x=gw(g,:);
    nx = length(x);

    nx=ceil(nx/nwin)*nwin;
    if nx>length(x) 
        x=[x, zeros(1,nx-length(x))]; % zerofill
    end

    nseg1=nx/nwin; 
    xseg=zeros(nseg1*2-1,nwin*os); 

    xseg(1:2:end,1:nwin)=reshape(x,[nwin,nseg1])';
    if nseg1>1
        xseg(2:2:end,1:nwin)=reshape(x(1+nwin/2:end-nwin/2),[nwin,nseg1-1])';
    end

    xseg_dc=mean(xseg,2);
    xseg=xseg-xseg_dc(:,ones(1,nwin*os));

    if nseg1>1 % WARNING: thisintroduces inconsistency between short and long sequences in term os the peak amplitudes
        cwin=0.5*(1-cos(2*pi*(1:nwin)/nwin));
        xseg(:,1:nwin)=xseg(:,1:nwin).*cwin(ones(size(xseg,1),1),:);
    end

    fseg=abs(fft(xseg,[],2));
    fseg=fseg(:,1:end/2); 
        
    if nseg1>1 
        gs = [gs; mean(fseg.^2).^0.5]; % sos
        %figure; plot(faxis(1:nfmax),sum(fseg(:,1:nfmax).^2).^0.5);
    else
        gs = [gs; abs(fseg)]; % add abs
    end    
end

 
figure; plot(faxis(1:nfmax),gs(:,1:nfmax));
%%figure; plot(faxis(1:nfmax),sum(gs(:,1:nfmax))); % abs-sum
hold on; plot(faxis(1:nfmax),sum(gs(:,1:nfmax).^2).^0.5); % sos
xlabel('frequency / Hz');

disp(['forbidden frequencies: ', num2str(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency)])
forb_freq=[];
if ischar(ascName)
    hold on; 
    for i=1:length(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency)
        if ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i)>0
            start_fre = ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i)-ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceBandwidth(i)/2;
            end_fre = ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i)+ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceBandwidth(i)/2;
            disp(['forbidden frequency range: ', num2str(start_fre), ...
                ' to ', num2str(end_fre)])
            forb_freq=[forb_freq;start_fre,end_fre];
            xline(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i),'-');
            xline(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i)-ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceBandwidth(i)/2,'--');
            xline(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i)+ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceBandwidth(i)/2,'--');
        end
    end
end

legend({'Gx','Gy','Gz','Gtot'});
% check frequency
freq_main = 1 / (2 * esp);  % main epi frequency in Hz
disp(['esp: ', num2str(esp*1e3), ' main epi frequency: ', num2str(freq_main)])

for loopi = 1:size(forb_freq,1)
    if (freq_main>forb_freq(loopi,1))&&(freq_main<forb_freq(loopi,2))
        status = 0;
        break;
    end
end
if (status ==0)
    disp(['Please check the ESP!!!!!!!!'])
else
    disp(['ESP check passed'])
end