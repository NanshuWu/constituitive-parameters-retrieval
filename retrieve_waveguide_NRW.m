%#######################################################################
% multiple frequency in a waveguide environment
% sample has thickness delta, total length of the waveguide is d1+delta+d2
% 
% Nicholson-Ross-Weir Conversion Method（NRW）
% 
% Advantages of NRW method
% ·Fast, non-iterative
% ·Applicable to waveguides and coaxial line
% 
% Disadvantages of NRW method
% ·Divergence at frequencies corresponding to multiples of one-half wavelength
% ·Short sample should be used （less than a quarter of wavelength）
% ·Not suitable for low loss materials
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
% | low loss solids+    |   ITER  |  Slow |  Medium  |
% | short+              |         |       |          |
% | magnetics/non-mag.  |         |       |          |
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
filename_2='data/wr90 s21 fr-4 d1 0 d2 0 delta 73.txt';
is_simu=0;% simulation, choose 1; experiment choose 0;
a=22.86e-3; % waveguide size
h=10.16e-3; % waveguide size
delta=2e-3; % sample thickness
d1=82e-3; % distance from port 1 to sample front face
d2=81e-3; % distance from port 2 to sample end face
n=-5:5; % branched
non_magnetic=0;
smooth_method=0;
maxe_k=8;
ave_n=10;
gap_correction=0;
hg=0e-3;% gap in height
ds=1;%data sample span
plot_type=1; % plot type
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
end
%#######################################################################
f=f.*1E9;
omega=2*pi*f;
beta0=sqrt((omega./c).^2-(pi/a).^2);
Z0=omega.*mu0./beta0;
lambda0=c./f;% working frequency wavelength
lambdac=2*a;
lambda0g=lambda0./sqrt(1-(f_cutoff./f).^2);
group_delay_calculated=zeros(length(f),length(n));
R1=exp(-1i.*beta0.*d1);R2=exp(-1i.*beta0.*d2);
s11=(S11_mag.*exp(1i.*S11_phase))./(R1.^2);
s21=(S21_mag.*exp(1i.*S21_phase))./(R1.*R2);
V1=s21+s11;V2=s21-s11;
X=(1-V1.*V2)./(V1-V2);
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
group_delay_measured=-gradient(unwrap(angle(T)))./gradient(f)/pi/2;% group delay measured
Kamp=abs(1./T);Kphase=unwrap(angle(1./T));
K=log(Kamp)+1i.*(Kphase+2.*pi.*n);
id_1=1i.*(K./2./pi./delta);
id_2=-1i.*(K./2./pi./delta);
id(real(id_1)>0)=id_1(real(id_1)>0);
id(real(id_1)<=0)=id_2(real(id_1)<=0);
id=reshape(id(:),length(f),[]);
switch non_magnetic
    case 0
        mu_ret=lambda0g.*id.*((1+ri)./(1-ri));
        eps_ret=lambda0.^2.*(id.^2+(1/lambdac).^2)./mu_ret;
    case 1
        mu_ret=ones(length(f),length(n));
        eps_ret=lambda0.^2.*(id.^2+(1/lambdac).^2)./mu_ret;
end
% group delay calculated 
for ni=1:length(n)
partial_p=real(sqrt((eps_ret(:,ni).*mu_ret(:,ni)./(lambda0.^2))-(1/lambdac).^2));
group_delay_calculated(:,ni)=delta.*gradient(partial_p)./gradient(f);
end
%#######################################################################
% Curve fitting evaluation
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
% Air gap correction
if gap_correction==1
    hm=h-hg;
    % for eps
    eps_cr=real(eps_ret).*(hm./(h-(h-hm).*real(eps_ret)));
    ltanm=-imag(eps_ret)./real(eps_ret);
    ltanc=ltanm.*(h./(h-(h-hm).*real(eps_ret)));
    eps_ci=-1.*ltanc.*eps_cr;
    eps_ret=eps_cr+1i.*eps_ci;
    % for mu
    mu_cr=real(mu_ret).*(h./hm)-((h-hm)./hm);
    mu_ci=imag(mu_ret).*h./hm;
    mu_ret=mu_cr+1i.*mu_ci;
end
%#######################################################################
% Smoothing 
switch smooth_method
    case 1 %The method of maximum entropy
        for ni=1:length(n)
            eps_ret_r=real(eps_ret(:,ni));
            eps_ret_i=imag(eps_ret(:,ni));
            maxe_k; % k power least squares fit
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
    case 2 % A general method
        eps_ret=movmean(eps_ret,ave_n);
        mu_ret=movmean(mu_ret,ave_n);
    otherwise
end
%#######################################################################
% correct store format
eps_ret=conj(eps_ret);
mu_ret=conj(mu_ret);
%#######################################################################
% figures
switch plot_type
    case 1
        fg1=figure(1);
        subplot(161); 
        plot(f/1e9,real(eps_ret(:,selected_branch)),'Color','k'); 
        xlabel('Frequency in GHz') 
        ylabel('Re(\epsilon)') 
        ylim([-10 10])
        subplot(162) 
        plot(f/1e9,imag(eps_ret(:,selected_branch)),'Color','k'); 
        xlabel('Frequency in GHz') 
        ylabel('Im(\epsilon)')
        ylim([-10 10])
        subplot(163);
        plot(f/1e9,real(mu_ret(:,selected_branch)),'Color','k'); 
        xlabel('Frequency in GHz')
        ylabel('Re(\mu)') 
        ylim([-10 10])
        subplot(164) 
        plot(f/1e9,imag(mu_ret(:,selected_branch)),'Color','k'); 
        xlabel('Frequency in GHz') 
        ylabel('Im(\mu)')
        ylim([-10 10])
        subplot(165)
        plot(f/1e9,imag(eps_ret(:,selected_branch))./real(eps_ret(:,selected_branch)),'Color','k'); 
        xlabel('Frequency in GHz') 
        ylabel('\epsilon tan(\delta)')
        subplot(166)
        plot(f/1e9,imag(mu_ret(:,selected_branch))./real(mu_ret(:,selected_branch)),'Color','k'); 
        xlabel('Frequency in GHz') 
        ylabel('\mu tan(\delta)')
        fg1.Position = [500 500 1600 500];
        fg1.Color='w';
    case 2
        legend_str={'Re(\epsilon)','Im(\epsilon)','\epsilon tan(\delta)','Re(\mu)','Im(\mu)','\mu tan(\delta)'};
        fg1=figure(1);
        plot(f/1e9,real(eps_ret(:,selected_branch)),...
            f/1e9,imag(eps_ret(:,selected_branch)),...
            f/1e9,imag(eps_ret(:,selected_branch))./real(eps_ret(:,selected_branch)),...
            f/1e9,real(mu_ret(:,selected_branch)),...
            f/1e9,imag(mu_ret(:,selected_branch)),...
            f/1e9,imag(mu_ret(:,selected_branch))./real(mu_ret(:,selected_branch))); 
        legend(legend_str,'Location','northeastoutside');
        xlabel('Frequency in GHz')
        fg1.Position = [500 500 800 500];
        fg1.Color='w';
end
% toc
