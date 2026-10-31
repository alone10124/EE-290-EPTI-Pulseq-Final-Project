function ky_acq = calc_Pattern_MS_WJG(Nseg,Rseg,ETL,Rpe)
%more general original pattern
N_total = ETL;
dk_pe = 1;
N_blipdown = round(Rseg/Rpe);
dk_blipdown = dk_pe*Rpe;

dk_pre = dk_pe*Rseg;

ky_acq = zeros(Nseg,N_total);
ky_pre=zeros(1,Nseg);
blip = zeros(Nseg,N_total+1);

for diff_s = 1:Nseg
    if(mod(diff_s,2)==1) 
        ky_pre(diff_s) = diff_s*dk_pre;
        ky_pre(diff_s) = ky_pre(diff_s)+mod(ky_pre(diff_s),2);
        for diff_n = 1:N_total
            if( mod(diff_n,2*N_blipdown) == N_blipdown )
                blip(diff_s,diff_n+1)=(N_blipdown-1)*dk_blipdown-dk_pe*Rpe/2;
            elseif( mod(diff_n,2*N_blipdown) == 0)
                blip(diff_s,diff_n+1)=(N_blipdown-1)*dk_blipdown+dk_pe*Rpe/2;
            else
                blip(diff_s,diff_n+1) = -dk_blipdown;
            end
            ky_acq(diff_s,diff_n) = ky_pre(diff_s)+squeeze(sum(blip(diff_s,1:diff_n),2));
        end
    else
        ky_pre(diff_s) = (diff_s)*dk_pre-dk_pe*Rpe/2;
        ky_pre(diff_s) = ky_pre(diff_s)+mod(ky_pre(diff_s),2);
        for diff_n = 1:N_total
            if( mod(diff_n,2*N_blipdown) == N_blipdown )
                blip(diff_s,diff_n+1)=(N_blipdown-1)*dk_blipdown+dk_pe*Rpe/2;
            elseif( mod(diff_n,2*N_blipdown) == 0)
                blip(diff_s,diff_n+1)=(N_blipdown-1)*dk_blipdown-dk_pe*Rpe/2;
            else
                blip(diff_s,diff_n+1) = -dk_blipdown;
            end
            ky_acq(diff_s,diff_n) = ky_pre(diff_s)+squeeze(sum(blip(diff_s,1:diff_n),2));
        end
    end
end