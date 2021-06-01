%#######################################################################
% multiple frequency in a waveguide environment
% sample has thickness delta, total length of the waveguide is d1+delta+d2
% 
% New Non-Iterative Conversion Method
% 
% suitable for permittivity calculation for the case of non-magnetic material
% the method has the advantage of being stable over a whole range of
% frequencies for an arbitrary sample length.
% 
% Advantages of new non-iterative method
% · Smooth permittivity results, no divergence.
% · Accurate.
% · Arbitrary length of samples can be used.
% · Fast, non-iterative.
% · No initial guess needed.
% 
% Disadvantages of new non-iterative method
% • Applicable for permittivity measurement only.
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
filename_1='FR4_may30/FR4-SAMPLE3-D1=82-D2=81-DELTA=2.S2P';
filename_2='data/wr90 s21 fr-4 d1 0 d2 0 delta 73.txt';
% tic
a=22.86e-3; % waveguide size
h=10.16e-3; % waveguide size
d1=82e-3; % distance from port 1 to sample front surface 
d2=81e-3;% distance from port 2 to sample back surface
delta=2e-3; % sample thickness
is_simu=0; % for simulation choose 1, for experiment choose 0;
n=-5:5; % branches
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
s21=(S21_mag.*exp(1i.*S21_phase)).*exp(1i.*beta0.*(d1+d2));
s11=(S11_mag.*exp(1i.*S11_phase)).*exp(1i.*2.*beta0.*d1);
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
% New method calculate K
Kamp=abs(1./T);Kphase=unwrap(angle(1./T));
% Kamp=abs(1./T);Kphase=angle(1./T);
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
% %#######################################################################
% group delay measured
group_delay_measured=-gradient(unwrap(angle(T)))./gradient(f)/pi/2;
% group delay calculated
for ni=1:length(n)
partial_p=sqrt((eps_ret(:,ni).*mu_ret(:,ni)./lambda0./lambda0)-(1/lambdac).^2);
group_delay_calculated(:,ni)=delta.*gradient(partial_p)./gradient(f);
end
%#######################################################################
% evaluate goodness of fit
R_square=zeros(1,length(n));
for ni=1:length(n)
    ymean=mean(real(group_delay_measured));
    SStot=sum((ymean-real(group_delay_measured)).^2);
    SSreg=sum((real(group_delay_calculated(:,ni))-ymean).^2);
    SSres=sum((real(group_delay_measured)-real(group_delay_calculated(:,ni))).^2);
    R_square(ni)=1-(SSres)./(SStot);
end
%#######################################################################
% branch selection
selected_branch_base=find(R_square>0,1,'first');
selected_branch_all=find(R_square>0);
bss=find(min(abs(sum(diff(real(eps_ret(:,selected_branch_all))))))==abs(sum(diff(real(eps_ret(:,selected_branch_all))))));
selected_branch=selected_branch_base+bss-1;
%#######################################################################
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
% toc
