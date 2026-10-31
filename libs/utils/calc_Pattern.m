function ky_acq = calc_Pattern(Nseg,Rseg,ETL,Rpe,PFshift)

ky_acq = zeros(Nseg,ETL);
Necho_each = round(Rseg/Rpe/2);
Necho_PFshift = floor(PFshift/Rpe/2);
if Nseg>1
    Start_point = (Rseg:Rseg:Nseg*Rseg) - Rseg/2;
else
    Start_point = (Rseg+PFshift)/2 - PFshift;
end
change_point = (1:2:50)*Necho_each-Necho_PFshift;
for shot_num = 1:Nseg
    PE_currect = Start_point(shot_num);
    Sign_PEdirection = -1;
    for i = 1:ETL
      ky_acq(shot_num,i) = PE_currect;
      if i>2
      if sum((i-2)==change_point(:))>0
          PE_currect = PE_currect+Sign_PEdirection*Rpe/2;
          Sign_PEdirection = -1*Sign_PEdirection;
      else
          PE_currect = PE_currect+Sign_PEdirection*Rpe;
      end
      end
    end
end