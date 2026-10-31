function ky_acq = calc_Pattern_MS5(Nseg,Npe,ETL,Rpe)
N_total = ETL;
dk_pe = 1;
Rseg2 = round(Npe/((Nseg-1)*2+1));   % the size of the center region
Rseg = round((Npe-Rseg2)/(Nseg-1));
dk_blipdown = dk_pe*Rpe;
Rpe2 = Rpe; %4, 6
ky_acq = zeros(Nseg,N_total);
blip = zeros(Nseg,N_total+1);
start_pe = 0;

start_pe_array=zeros(1,Nseg);
start_pe_array(1) = Rseg;
for loopi = 2:Nseg
    if (loopi ==(floor(Nseg/2)+1))
       start_pe_array(loopi) = start_pe_array(loopi-1)+round(Rseg2/2);
    elseif(loopi ==(floor(Nseg/2)+2))
        start_pe_array(loopi) = start_pe_array(loopi-1)+Rseg+round(Rseg2/2);
    else
        start_pe_array(loopi) = (start_pe_array(loopi-1)+Rseg);
    end
end
start_pe_array = floor(start_pe_array);
for diff_s = 1:Nseg
    if diff_s == (floor(Nseg/2)+1) %the center segment
        start_pe =start_pe_array(diff_s);
        PE_currect = start_pe;
        PE_currect = PE_currect+mod(PE_currect,2);
        Sign_PEdirection = -1;
        skip=0;
        pe_end1 = round(PE_currect-Rseg2);
        pe_end2 = round(PE_currect+Rseg2);
        for diff_n = 1:N_total
            if(((PE_currect>=Npe/2+round(Rseg2/2)-Rpe2)||(PE_currect<Npe/2-round(Rseg2/2)+Rpe2/2))&&(skip==1)) 
                Sign_PEdirection = Sign_PEdirection*-1;
                PE_currect = PE_currect+Sign_PEdirection*Rpe2/2;
                skip=0;
             elseif(((PE_currect>Npe/2+round(Rseg2/2)-Rpe2)||(PE_currect<Npe/2-round(Rseg2/2)+Rpe2))&&(skip==1))   
                PE_currect = PE_currect+Sign_PEdirection*Rpe2/2;
                Sign_PEdirection = Sign_PEdirection*-1;
                skip=0;
            else
                PE_currect = PE_currect+Sign_PEdirection*Rpe2;
                skip=1;
            end        
            ky_acq(diff_s,diff_n) = PE_currect;
        end
        start_pe = start_pe+round(Rseg2/2);

    else
        start_pe =start_pe_array(diff_s);
        start_pe =  start_pe;
        PE_currect = round(start_pe-round(Rseg/2));
        PE_currect = PE_currect+mod(PE_currect,2);

        Sign_PEdirection = -1;
        pe_end1 = max(round(PE_currect-round(Rseg/2)),1);
        pe_end2 = min(round(PE_currect+round(Rseg/2)),Npe);
        skip=0;
        for diff_n = 1:N_total
            if(((PE_currect>=pe_end2-Rpe/2)||(PE_currect<=Rpe/2+pe_end1))&&(skip==1)) 
                Sign_PEdirection = Sign_PEdirection*-1;
                PE_currect = PE_currect+Sign_PEdirection*Rpe/2;
                skip=0;
            elseif(((PE_currect>=pe_end2-Rpe)||(PE_currect<=Rpe))&&(skip==1)) 
                PE_currect = PE_currect+Sign_PEdirection*Rpe/2;
                Sign_PEdirection = Sign_PEdirection*-1;
                skip=0;
            else
                PE_currect = PE_currect+Sign_PEdirection*Rpe;
                skip=1;
            end  
            ky_acq(diff_s,diff_n) = PE_currect;
        end
    end

end
