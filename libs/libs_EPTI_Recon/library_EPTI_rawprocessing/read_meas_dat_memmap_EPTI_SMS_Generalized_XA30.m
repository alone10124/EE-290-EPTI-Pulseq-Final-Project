function [meas] = read_meas_dat_memmap_EPTI_SMS_Generalized_XA30(filename,nRepToRead,BeginRep,ascendingSlice_acq)
% reading raw data for EPTI
%% Jan 2018
% Fuyixue Wang
% make it work for EPTI data
% GE-SE data with different number of echoes processed
% before put into kspace encoding pattern
% compatable to rawdata when ICE is turned OFF
%% June 2021
% Fuyixue, Zijing 
% raw data reading for GE-SE EPTI without phasecor navigator, works for SMS
%% Jan 2023
% Fuyixue Wang, Zijing Dong
% generalized for single-shot and multi-shot data with or without phasecor
% navigator, optimized parameter extraction for both 3T and 7T data
%%
if nargin == 4
    ascendingSlice_acq = 0; 
end

save_fname = [filename(1:end-4) '_FirstRep.dat'];

if exist([save_fname(1:end-4) '_Raw.mat'],'file') ~= 2

    opt.SMSRefscan = 0; 
    opt.ReturnStruct=1;    % Zijing settings to read data w/o phasecor navigator
    opt.ReadMultipleRepetitions = 0;
    opt.SqueezeChannels = 0;
    opt.ReturnCellArray = 1; 
    opt.SkipRTFeedback = 0;
    opt.OfflineSave = 0;

    %readin first rep and save prot,evp,(patrefscan,patrefscan_phascor)
    meas_first = read_meas_dat_nonav_XA30(filename,opt);
    meas_first.nsample = size(meas_first.data,9);
    if meas_first.prot.lAccelFact3D>1
        SMSdata = 1;
    else
        SMSdata = 0;
    end
    meas_first.prot.SMSdata = SMSdata;
    meas_first.evp.NChaMeas = size(meas_first.data,3);
 
    if isempty(meas_first.evp.RawCol)  % ICE is off when scaning
        if meas_first.prot.sWipMemBlock_alFree(36) ~=0
            meas_first.evp.NLinMeas=floor(meas_first.prot.sWipMemBlock_alFree(36)*meas_first.evp.NLinMeas/100); % new sequence with arbitray ETL
        else
            meas_first.evp.NLinMeas=size(meas_first.data,2);
        end
        meas_first.evp.NAFLin=meas_first.prot.lAccelFactPE;
        meas_first.evp.NRefLin = meas_first.prot.lRefLinesPE;
        meas_first.evp.NAFPar=1;
        meas_first.evp.NFirstPar=1;
        meas_first.evp.NParMeas=1;
        meas_first.evp.NFirstRefPar=0;
        meas_first.evp.NRefPar=1;
    end

    if size(meas_first.data,5)>1
        lPhaseDim_GE=floor(meas_first.prot.sWipMemBlock_alFree(37)/100*meas_first.evp.NLinMeas);
        lPhaseDim_GE = lPhaseDim_GE - mod(lPhaseDim_GE,2);
        meas_first.prot.pf_echo = (lPhaseDim_GE+0.01) / meas_first.evp.NLinMeas;
    else
        meas_first.prot.pf_echo = 0;
    end
    
    meas_first.evp.SegOneLines = 1:2:meas_first.evp.NLinMeas;
    meas_first.evp.SegTwoLines = 2:2:meas_first.evp.NLinMeas;

    if SMSdata == 1 % find non empty slice and grab the first SlcGroup slc position
        [SlcMask] = SlcMaskGenerator(meas_first);
        index = find(SlcMask == 1);
        FirstIndex = index(1);           
        Nslices = meas_first.evp.NSlcMeas;
        SlcsPerGroup = sum(SlcMask);
        NslicesEX = Nslices/SlcsPerGroup;
        meas_first.prot.SliceSep = sqrt((meas_first.prot.sSliceArray(1).sPosition_dSag- meas_first.prot.sSliceArray(1+SlcsPerGroup).sPosition_dSag)^2+...
        (meas_first.prot.sSliceArray(1).sPosition_dCor - meas_first.prot.sSliceArray(1+SlcsPerGroup).sPosition_dCor)^2+...
        (meas_first.prot.sSliceArray(1).sPosition_dTra - meas_first.prot.sSliceArray(1+SlcsPerGroup).sPosition_dTra)^2);   % distance between two slice

        meas_first.prot.sSliceArray = meas_first.prot.sSliceArray(1:SlcsPerGroup); % CONVENTION: slices position of the top slice group!!!
    else
        FirstIndex = 1;
    end

    if isfield(meas_first,'data_phascor1d') 
        meas_first.evp.NPhaseCorLines = size(meas_first.data_phascor1d,2);
        meas_first.prot.navigatorOn = 1;
    else
        meas_first.evp.NPhaseCorLines = 0;
        meas_first.prot.navigatorOn = 0;
    end
    
    if ascendingSlice_acq == 0
        deinterleave = strcmp(meas_first.prot.ucMultiSliceMode, 'MSM_INTERLEAVED');
    else
        meas_first.prot.ucMultiSliceMode = 'MSM_SEQUENTIAL';
        deinterleave = 2;
    end
    
    meas_first.data = []; 
    meas_first.data_phascor1d = []; 
    meas_first.smsrefscan_phascor = [];
    SaveRawData(meas_first,save_fname);

else
%     disp('Load First Rep Info')
%     tic
    meas_first = ReadRawData(save_fname);
%     toc
    if ascendingSlice_acq == 0
        deinterleave = strcmp(meas_first.prot.ucMultiSliceMode, 'MSM_INTERLEAVED');
    else
        meas_first.prot.ucMultiSliceMode = 'MSM_SEQUENTIAL';
        deinterleave = 2;
    end  
end
meas_first.evp.NChaMeas = meas_first.prot.iMaxNoOfRxChannels;
coil_index = 1:meas_first.evp.NChaMeas;

%% extract params from first rep

meas = meas_first;
SMSdata = meas.prot.SMSdata;
pf_echo = meas.prot.pf_echo;
clear meas_first;
if SMSdata == 1
    if meas.isvd == 1    
        if ~isempty(meas.prot.sWipMemBlock_adFree)
            meas.prot.sWipMemBlock_adFree(1) = meas.prot.sWipMemBlock_adFree(2);%*meas.prot.lAccelFactPE % FOV shift factor
            meas.prot.sWipMemBlock_adFree(2) = meas.prot.sWipMemBlock_adFree(3); % SMS factor
            
        else
            [prottemp] = read_meas_prot(filename);
            meas.prot.sWiPMemBlock_adFree(2) = prottemp.sWiPMemBlock_adFree(3); % shift in reduced FOV
            meas.prot.sWiPMemBlock_adFree(1) = prottemp.sWiPMemBlock_adFree(5); % SMS  !!%% original is (2)
            
        end
    else
       error('Error of the data: unrecognized data format');
    end
    if meas.prot.sWipMemBlock_adFree(1) > 10 %
        meas.prot.sWipMemBlock_adFree(1) = meas.prot.sWipMemBlock_adFree(2); % SMS factor
        meas.prot.sWipMemBlock_adFree(2) = meas.prot.sWipMemBlock_adFree(5); % FOV shift factor
    end
    NslicesEX =  meas.prot.sWipMemBlock_adFree(1);
else
    NslicesEX = 1;
end

if isempty(meas.evp.NColMeas)==0
    sData = [meas.evp.NColMeas, meas.evp.NLinMeas, meas.evp.NChaMeas,...
              1, meas.evp.NEcoMeas, 1, nRepToRead, 2, 1, meas.evp.NSlcMeas/NslicesEX ];
    sPhaseCor = sData; sPhaseCor(2) = meas.evp.NPhaseCorLines; sPhaseCor(5) = 1;

    nRead = sData(1);
    nPE = sData(2);
    nCoil = sData(3);
    nSlice = sData(10);
    nPhaseCor = sPhaseCor(2);
    nEco = sData(5);

    nPE_Raw = meas.evp.RawLin;
    nRep = meas.evp.RawRep;
    meas.evp.NRepMeas = nRepToRead;
else     
     sData = [meas.prot.lBaseResolution*2, meas.evp.NLinMeas, meas.evp.NChaMeas,...
     1, meas.evp.NEcoMeas, 1, nRepToRead, 2, 1, meas.evp.NSlcMeas/NslicesEX ];
    sPhaseCor = sData; sPhaseCor(2) = meas.evp.NPhaseCorLines;
    sPhaseCor(5)=1; % echo = 1
    if pf_echo==0
        sData(5)=1;
    end
    nRead = sData(1);
    nPE = sData(2);
    nCoil = meas.prot.iMaxNoOfRxChannels;
    nSlice = sData(10);
    nPhaseCor = sPhaseCor(2);
    nEco = sData(5);
    nPE_Raw = size(meas.evp.SegOneLines,2)+size(meas.evp.SegTwoLines,2);
    nRep=nRepToRead;
    meas.evp.NRepMeas = nRepToRead;
end

SegOneLines = meas.evp.SegOneLines;
SegTwoLines = meas.evp.SegTwoLines;


%% calculate parameters needed for mmap and readout and preallocate matrices
meas.data = single(zeros(sData));
if  meas.prot.navigatorOn == 1
    meas.data_phascor1d = single(zeros(sPhaseCor));
end

if pf_echo>0
    Npe_1echo= floor(pf_echo*nPE_Raw);
    Npe_2echo= nPE_Raw;
    linesperrep = round(Npe_2echo+Npe_1echo+nPhaseCor)*nCoil*nSlice;
else
    linesperrep=round(nPE_Raw+nPhaseCor)*nCoil*nSlice ;
end

%% mmap
nRep=max(nRep,BeginRep);
m = mmap_mdh_noheader_rep_v2(filename,linesperrep,nRep,nRead,nCoil,meas.data_begin,meas.isvd);
%% Readout

ICE_RAWDATA_SCALE       = 131072.0;  % 64 ^ 3 / 2
K_ICE_AMPL_SCALE_FACTOR = 1;%80 * 20 * ICE_RAWDATA_SCALE / 65536;
if meas.isvd == 1
    mdh_length_float32 = 48; % per PE line 
    mdh_ch_length_float32 = 8; % per chanel in each PE line
else
    mdh_length_float32 = 0;
    mdh_ch_length_float32 = 32;
end


for RepCount = 1:nRepToRead

    k = m.Data(BeginRep+(RepCount-1)).kdata;
    k = reshape(k,[],linesperrep/nCoil)*K_ICE_AMPL_SCALE_FACTOR;   
    k = k(1+mdh_length_float32:end,:);
    k = reshape(k,[],linesperrep);    
    k = k(1+mdh_ch_length_float32:2:end,:) + 1i*k(2+mdh_ch_length_float32:2:end,:); %32 headers
    

    if pf_echo>0
        k = reshape(k, [nRead, nCoil, (Npe_1echo+Npe_2echo)+ nPhaseCor, nSlice]);
    else
        k = reshape(k, [nRead, nCoil,  nPE_Raw + nPhaseCor, nSlice]);
    end
   
    if  meas.prot.navigatorOn == 1
        meas.data_phascor1d(:,1:2:end,:,:,:,:,RepCount,1,:,:) = permute(k(:,coil_index,1:2:nPhaseCor,:),[1 3 2 4]);
        meas.data_phascor1d(:,2:2:end,:,:,:,:,RepCount,2,:,:) = permute(k(end:-1:1,coil_index,2:2:nPhaseCor,:),[1 3 2 4]); % need to reverse k-line as assume EPI data
    end

    if pf_echo>0
        k0=zeros(nRead,nCoil,Npe_2echo,nEco,nSlice);
        k0(:,:,1:Npe_1echo,1,:) = reshape(k(:,:,nPhaseCor+1:nPhaseCor+Npe_1echo,:),nRead,nCoil,Npe_1echo,1,nSlice); 
        k0(:,:,:,2,:) = reshape(k(:,:,nPhaseCor+Npe_1echo+1:end,:),nRead,nCoil,Npe_2echo,1,nSlice); 
    else
        k0=zeros(nRead,nCoil,nPE_Raw,nEco,nSlice);
        k0(:,:,1:nPE_Raw,:,:) = reshape(k(:,:,nPhaseCor+1:end,:),nRead,nCoil,nPE_Raw,1,nSlice); 
    end
    k=k0;
    meas.data(:,SegOneLines,:,:,:,:,RepCount,1,:,:) = permute(k(:,coil_index,1:2:end,:,:),[1 3 2 4 5]);
    meas.data(:,SegTwoLines,:,:,:,:,RepCount,2,:,:) = permute(k(end:-1:1,coil_index,2:2:end,:,:),[1 3 2 4 5]); % need to reverse k-line as assume EPI data
end


if  meas.prot.navigatorOn == 1 
    if deinterleave == 1
        meas.data_phascor1d = mrir_image_slice_deinterleave(meas.data_phascor1d);
    elseif deinterleave == 2 % ascending need to reverse
        meas.data_phascor1d = meas.data_phascor1d(:,:,:,:,:,:,:,:,:,end:-1:1, :,:,:,:,:,:);
    end
end

function [evp,prot,patrefscan] = SiemensIPATConventionWorkAround(evp,prot,patrefscan) 
%% Modify MEAS file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% workaround for unusual Siemens convention #1:
if ( evp.NFirstRefLin == 0 ),
    evp.NFirstRefLin = mrir_ice_dimensions(patrefscan, 'lin') - evp.NRefLin + 1;
end;

% workaround for unusual Siemens convention #2:
% (possibly Siemens fills in with GRAPPA fewer lines than are in FFT, so
% last lines are effectively zero-padded; this could throw off SNR
% calculations, so by overriding this we force "mrir_epi_GRAPPA" to fill
% in same number of k-space lines as there are image lines.)
if ( ~isempty(evp.NAFLin) && (evp.NAFLin == 1) && (prot.ucPhasePartialFourier == 1) && (evp.NLinMeas < evp.NImageLins) ),
    jnotify
    keyboard
    evp.NLinMeas = evp.NImageLins;
end;

function [SlcMask] = SlcMaskGenerator(meas_first)

% figure out SMS factor and slc location

if  isfield(meas_first, 'patrefscan')
    R_inplane = meas_first.prot.lAccelFactPE;
else
    R_inplane = 1;
end

PElineSearchCount = 1; 
while PElineSearchCount <= R_inplane+1
    SlcMask = [];
    for ncoil = 1:size(meas_first.data,3)    % Zijing sreach coils in case there is empty coil
    a = meas_first.data(1,PElineSearchCount,ncoil,1,1,1,1,1,1,:);
    a = squeeze(a);
    SlcMask_tmp = length(a);
    for count = 1:length(a)
        SlcMask_tmp(count) = ~isempty(a{count});
    end
    if sum(SlcMask_tmp)>sum(SlcMask)
        SlcMask = SlcMask_tmp;
    end
    end
    
    if sum(SlcMask) ~= 0
        break; % found data!!
    else
        if (PElineSearchCount == R_inplane+1)
            disp('error!! cant find the data')
            keyboard
        end
        PElineSearchCount = PElineSearchCount+ 1; % keep looking!!
    end
end
    
    
    
     
  
  
  
  

