%#######################################################################
% multiple frequency in a waveguide environment
% sample has thickness delta, total length of the waveguide is d1+delta+d2
% 
% NIST Iterative
% 
% Advantages of NIST iterative method
% · Smooth permittivity results, no divergence.
% · Accurate.
% · Arbitrary length of samples can be used.
% · Robust for low loss and high loss materials.
% 
% Disadvantages of NIST iterative method
% · Applicable for permittivity measurement only.
% · Need initial guess of permittivity value.
% 
% algorithm details can be found in NIST Technical Notes 1355-R/1341 or R&S
% tech notes
% 
% how to choose algorithm
% +---------------------+---------+-------+----------+
% | Material length/    | Methods | Speed | Accuracy |
% | magnetic properties |         |       |          |
% +=====================+=========+=======+==========+
% | Lossy solids+       |   NRW   |  Fast |  Medium  |
% | short+              |         |       |          |
% | non-magnetics       |         |       |          |
% +---------------------+---------+-------+----------+
% | Lossy solids+       |   NRW   |  Fast |  Medium  |
% | short+              |         |       |          |
% | magnetics           |         |       |          |
% +---------------------+---------+-------+----------+
% | low loss solids+    |   NIST  |  Slow |   Good   |
% | long+               |         |       |          |
% | non-magnetics       |         |       |          |
% +---------------------+---------+-------+----------+
% | low loss solids+    |   NNI   |  Fast |   Good   |
% | long+               |         |       |          |
% | non-magnetics       |         |       |          |
% +---------------------+---------+-------+----------+
% 
% nanshu.wu@gmail.com
%#######################################################################
clear
clc
%#######################################################################
% tic
filename_1='dataset/FR4_d1_82_d2_81_delta_2.S2P';
filename_2='data/retrieve s21 ptfe lossy d 30 d1 75.15 d2 48.85.txt';
is_simu=0;% simulation, choose 1; experiment choose 0;
a=22.86e-3; % waveguide size
h=10.16e-3; % waveguide size
delta=2e-3; % sample thickness
d1=82e-3; % distance from port1 to sample front surface
d2=81e-3; % distance from port2 to sample back surface
n=-5:5;% branches
smooth_method=0;
maxe_k=8;
ave_n=10;
data_correction_method=1;
ds=1; % sampling interval
%#######################################################################
eps0=8.85e-12; 
mu0=4*pi*1e-7; 
c=1/sqrt(eps0*mu0);
f_cutoff=c/a/2;
legend_str=strings(length(n),1);
for ni=1:length(n)
legend_str{ni} = ['n=' num2str(n(ni))];
end
%#######################################################################
% import data
switch is_simu
    case 1
        [ f, S11_re, S11_img ]=readS( filename_1,3,ds);
        [ ~, S21_re, S21_img ]=readS( filename_2,3,ds);
        S11_mag=abs(S11_re+1i.*S11_img);
        S11_phase=angle(S11_re+1i.*S11_img);
        S21_mag=abs(S21_re+1i.*S21_img);
        S21_phase=angle(S21_re+1i.*S21_img);
    case 0
        [ f, S11_mag, S11_phase, S21_mag, S21_phase,S12_mag,S12_phase,S22_mag,S22_phase ]=readS_fromtest( filename_1,9,ds);
        f=f./1e9;
        S11_phase=deg2rad(S11_phase);
        S21_phase=deg2rad(S21_phase);
        S12_phase=deg2rad(S12_phase);
        S22_phase=deg2rad(S22_phase);
end
%#######################################################################
f=f.*1E9;
omega=2*pi*f;
beta0=sqrt((omega./c).^2-(pi/a).^2);
Z0=omega.*mu0./beta0;
lambda0=c./f;
lambdac=2*a;
lambda0g=lambda0./sqrt(1-(f_cutoff./f).^2);
group_delay_calculated=zeros(length(f),length(n));
R1=exp(-1i.*beta0.*d1);R2=exp(-1i.*beta0.*d2);
s11=(S11_mag.*exp(1i.*S11_phase))./(R1.^2);
s21=(S21_mag.*exp(1i.*S21_phase))./(R1.*R2);
X=(s11.^2-s21.^2+1)./2./s11;
ri_1=X+sqrt(X.^2-1);
ri_2=X-sqrt(X.^2-1);
if sum(abs(ri_1)<=1)==length(f)
    ri=ri_1;
elseif sum(abs(ri_1)<=1)==0
    ri=ri_2;
else
    ri(abs(ri_1)<=1)=ri_1(abs(ri_1)<=1);
    ri(abs(ri_1)>1)=ri_2(abs(ri_1)>1);
end
if size(ri,2)>length(n)
    ri=ri';
end
T=(s11+s21-ri)./(1-ri.*(s11+s21));
group_delay_measured=-gradient(unwrap(angle(T)))./gradient(f)/pi/2;% Group Delay measured
Kamp=abs(1./T);Kphase=unwrap(angle(1./T));
K=log(Kamp)+1i.*(Kphase+2.*pi.*n);

id_1=1i.*(K./2./pi./delta);
id_2=-1i.*(K./2./pi./delta);
id(real(id_1)>0)=id_1(real(id_1)>0);
id(real(id_1)<=0)=id_2(real(id_1)<=0);
if size(id,2)>length(n)
    id=id';
end
id=reshape(id,length(f),[]);

mu_ret=ones(length(f),length(n));
eps_ret=lambda0.^2.*(id.^2+(1./lambdac).^2)./mu_ret;

% Calculate Group Delay
for ni=1:length(n)
    partial_p=sqrt(((eps_ret(:,ni).*mu_ret(:,ni))./(lambda0.^2))-(1./lambdac).^2);
    group_delay_calculated(:,ni)=delta.*gradient(partial_p)./gradient(f);
end
%#######################################################################
%Evaluate goodness of fit
R_square=zeros(1,length(n));
for ni=1:length(n)
    ymean=mean(real(group_delay_measured));
    SStot=sum((ymean-real(group_delay_measured)).^2);
    SSreg=sum((real(group_delay_calculated(:,ni))-ymean).^2);
    SSres=sum((real(group_delay_measured)-real(group_delay_calculated(:,ni))).^2);
    R_square(ni)=1-(SSres)./(SStot);
end
%#######################################################################
% Branch Selection
selected_branch_base=find(R_square>0,1,'first');
selected_branch_all=find(R_square>0);
bss=find(min(abs(sum(diff(real(eps_ret(:,selected_branch_all))))))==abs(sum(diff(real(eps_ret(:,selected_branch_all))))));
selected_branch=selected_branch_base+bss-1;
%#######################################################################
% Data correction
switch data_correction_method
    case 1
        s21_nr=s21.*(R1.*R2);
        s11_nr=s11.*(R1.^2);
        if is_simu==1
            s22=s11_nr;
            s12=s21_nr;% for reciprocal system
        else
            s12=S12_mag.*exp(1i.*S12_phase);
            s22=S22_mag.*exp(1i.*S22_phase);
        end
        h=1e-8;
        
        for fi=1:length(f)
            eps_g=eps_ret(fi,selected_branch);
            for ii=1:50
            Errf=0;
            EffF=0;
            eps_gr=real(eps_g);
            eps_gi=imag(eps_g);
            propc0=1i*beta0(fi);
            propc=@(eps_g) 1i.*sqrt(eps_g*eps0*mu0*omega(fi)^2-(2*pi/lambdac)^2);
            Gam=@(propc) (propc0-propc)/(propc0+propc);
            transcoef=@(propc) exp(-1*propc*delta);

            fepsr1=@(transcoef,Gam) s21_nr(fi)*s12(fi)-s11_nr(fi)*s22(fi)-exp(-2*propc0*(d1+d2))*((transcoef^2-Gam^2)/(1-(Gam*transcoef)^2));
            
            eps_gr_pd=((eps_gr+h)+1i*eps_gi);
            eps_gr_nd=((eps_gr-h)+1i*eps_gi);
            eps_gi_pd=(eps_gr+1i*(eps_gi+h));
            eps_gi_nd=(eps_gr+1i*(eps_gi-h));

            Jrer=real((fepsr1(transcoef(propc(eps_gr_pd)),Gam(propc(eps_gr_pd)))-fepsr1(transcoef(propc(eps_gr_nd)),Gam(propc(eps_gr_nd))))/h/2);
            Jrei=real((fepsr1(transcoef(propc(eps_gi_pd)),Gam(propc(eps_gi_pd)))-fepsr1(transcoef(propc(eps_gi_nd)),Gam(propc(eps_gi_nd))))/h/2);
            Jier=imag((fepsr1(transcoef(propc(eps_gr_pd)),Gam(propc(eps_gr_pd)))-fepsr1(transcoef(propc(eps_gr_nd)),Gam(propc(eps_gr_nd))))/h/2);
            Jiei=imag((fepsr1(transcoef(propc(eps_gi_pd)),Gam(propc(eps_gi_pd)))-fepsr1(transcoef(propc(eps_gi_nd)),Gam(propc(eps_gi_nd))))/h/2);
            feps_g=[real(fepsr1(transcoef(propc(eps_g)),Gam(propc(eps_g))));imag(fepsr1(transcoef(propc(eps_g)),Gam(propc(eps_g))))];
            ErrF=abs(fepsr1(transcoef(propc(eps_g)),Gam(propc(eps_g))));
            J=[Jrer,Jrei;Jier,Jiei];
            injav=inv(J);
            deps=-injav*feps_g;
            if sum(abs(deps))<h
                eps_ret(fi,selected_branch)=eps_g;
                break;
            end
            eps_g=eps_gr+deps(1)+1i*(eps_gi+deps(2));
            end
        end
    case 2
        s21_nr=s21.*exp(-1i.*beta0.*(d1+d2));
        s11_nr=s11.*exp(-1i.*2.*beta0.*d1);
        if is_simu==1
            s22=s11_nr;
            s12=s21_nr;% for reciprocal system
        else
            s12=S12_mag.*exp(1i.*S12_phase);
            s22=S22_mag.*exp(1i.*S22_phase);
        end
        h=1e-8;

        for fi=1:length(f)
            eps_g=eps_ret(fi,selected_branch);
            for ii=1:50
            Errf=0;
            eps_gr=real(eps_g);
            eps_gi=imag(eps_g);
            propc0=1i*beta0(fi);
            propc=@(eps_g) 1i.*sqrt(eps_g*eps0*mu0*omega(fi)^2-(2*pi/lambdac)^2);
            Gam=@(propc) (propc0-propc)/(propc0+propc);
            transcoef=@(propc) exp(-1*propc*delta);
            fepsr2=@(transcoef,Gam) 0.5*(s21_nr(fi)+s12(fi))-exp(-1*propc0*(d1+d2))*((transcoef*(1-Gam^2))/(1-(transcoef*Gam)^2));
            eps_gr_pd=((eps_gr+h)+1i*eps_gi);
            eps_gr_nd=((eps_gr-h)+1i*eps_gi);
            eps_gi_pd=(eps_gr+1i*(eps_gi+h));
            eps_gi_nd=(eps_gr+1i*(eps_gi-h));

            Jrer=real((fepsr2(transcoef(propc(eps_gr_pd)),Gam(propc(eps_gr_pd)))-fepsr2(transcoef(propc(eps_gr_nd)),Gam(propc(eps_gr_nd))))/h/2);
            Jrei=real((fepsr2(transcoef(propc(eps_gi_pd)),Gam(propc(eps_gi_pd)))-fepsr2(transcoef(propc(eps_gi_nd)),Gam(propc(eps_gi_nd))))/h/2);
            Jier=imag((fepsr2(transcoef(propc(eps_gr_pd)),Gam(propc(eps_gr_pd)))-fepsr2(transcoef(propc(eps_gr_nd)),Gam(propc(eps_gr_nd))))/h/2);
            Jiei=imag((fepsr2(transcoef(propc(eps_gi_pd)),Gam(propc(eps_gi_pd)))-fepsr2(transcoef(propc(eps_gi_nd)),Gam(propc(eps_gi_nd))))/h/2);

            feps_g=[real(fepsr2(transcoef(propc(eps_g)),Gam(propc(eps_g))));imag(fepsr2(transcoef(propc(eps_g)),Gam(propc(eps_g))))];
            
            ErrF=abs(fepsr2(transcoef(propc(eps_g)),Gam(propc(eps_g))));
            Errf=Errf+(abs(feps_g(1))+abs(feps_g(2))); 
            J=[Jrer,Jrei;Jier,Jiei];
            injav=inv(J);
            deps=-injav*feps_g;
            if sum(abs(deps))<h
                eps_ret(fi,selected_branch)=eps_g;
                break;
            end
                eps_g=eps_gr+deps(1)+1i*(eps_gi+deps(2));       
            end
        end
    otherwise
end
%#######################################################################
% Smoothing
switch smooth_method
    case 1 % The method of maximum entropy
        for ni=1:length(n)
            eps_ret_r=real(eps_ret(:,ni));
            eps_ret_i=imag(eps_ret(:,ni));
    %         smooth_k=8; % k power least squares fit
            aki=ones(length(f),maxe_k);
            Moments=zeros(1,maxe_k);
            % Moments(1)=1;
            for akii=1:maxe_k
               aki(:,akii)=aki(:,akii).*omega.^(akii);
            end
            aki=aki';
            for ki=1:maxe_k
                for fi=1:length(f)
                    Moments(ki)=Moments(ki)+aki(ki,fi).*eps_ret_r(fi);
                end
            end

            Atran=aki';
            matC=aki*Atran;
            matD=inv(matC);
            matE=matD*Moments';
            eps_sr=Atran*matE;

            Moments=zeros(1,maxe_k);
            for ki=1:maxe_k
                for fi=1:length(f)
                    Moments(ki)=Moments(ki)+aki(ki,fi).*eps_ret_i(fi);
                end
            end
            Atran=aki';
            matC=aki*Atran;
            matD=inv(matC);
            matE=matD*Moments';
            eps_si=Atran*matE;

            eps_ret(:,ni)=eps_sr+1i*eps_si;
        end
    case 2 % A general smooth method 
        eps_ret=smoothdata(eps_ret,'movmean',ave_n);
end
%#######################################################################
% Figures
fg1=figure(1);
subplot(151);
plot(f/1e9,real(mu_ret(:,selected_branch))); 
xlabel('Frequency in GHz')
ylabel('Re(\mu)') 
ylim([-10 10])
subplot(152) 
plot(f/1e9,imag(conj(mu_ret(:,selected_branch)))); 
xlabel('Frequency in GHz') 
ylabel('Im(\mu)')
ylim([-10 10])
subplot(153); 
plot(f/1e9,real(eps_ret(:,selected_branch))); 
xlabel('Frequency in GHz') 
ylabel('Re(\epsilon)') 
ylim([-10 10])
subplot(154) 
plot(f/1e9,imag(conj(eps_ret(:,selected_branch)))); 
xlabel('Frequency in GHz') 
ylabel('Im(\epsilon)')
ylim([-10 10])
subplot(155)
plot(f/1e9,imag(conj(eps_ret(:,selected_branch)))./real(eps_ret(:,selected_branch))); 
xlabel('Frequency in GHz') 
ylabel('tan(\delta)')
ylim([0 1])
fg1.Position = [600 600 1500 600];
