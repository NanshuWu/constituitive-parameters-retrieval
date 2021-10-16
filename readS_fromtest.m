function [f,S11_mag,S11_phase,S21_mag,S21_phase,S12_mag,S12_phase,S22_mag,S22_phase]=readS_fromtest(filename,p,s)
%% Read S params and export as linear amplitude and radical phase angle
fileID = fopen(filename);
c2 = textscan(fileID,'# %s %s %s %s %f',1,'CommentStyle', '!');
c1 = textscan(fileID,'%f','headerlines',1);
freq_unit = c2{1,1};
data_type = c2{1,3};
f = c1{1,1}(1:s*p:size(c1{1,1}));
switch lower(char(freq_unit))
    case 'hz'
    case 'khz'
        f=f.*1e3;
    case 'mhz'
        f=f.*1e6;
    case 'ghz'
        f=f.*1e9;
    case 'thz'
        f=f.*1e12;
    case 'phz'
        f=f.*1e15;
end
switch lower(char(data_type))
    case 'db'
        S11_mag = c1{1,1}(2:s*p:size(c1{1,1}));
        S11_phase = c1{1,1}(3:s*p:size(c1{1,1}));
        S21_mag = c1{1,1}(4:s*p:size(c1{1,1}));
        S21_phase = c1{1,1}(5:s*p:size(c1{1,1}));
        S12_mag = c1{1,1}(6:s*p:size(c1{1,1}));
        S12_phase = c1{1,1}(7:s*p:size(c1{1,1}));
        S22_mag = c1{1,1}(8:s*p:size(c1{1,1}));
        S22_phase = c1{1,1}(9:s*p:size(c1{1,1}));
        S11_mag=10.^(S11_mag./20);% dB to linear
        S21_mag=10.^(S21_mag./20);% dB to linear
        S12_mag=10.^(S12_mag./20);% dB to linear
        S22_mag=10.^(S22_mag./20);% dB to linear
        S11_phase=deg2rad(S11_phase);% degree to rad
        S21_phase=deg2rad(S21_phase);% degree to rad
        S12_phase=deg2rad(S12_phase);% degree to rad
        S22_phase=deg2rad(S22_phase);% degree to rad
    case 'ma'
        S11_mag = c1{1,1}(2:s*p:size(c1{1,1}));
        S11_phase = c1{1,1}(3:s*p:size(c1{1,1}));
        S21_mag = c1{1,1}(4:s*p:size(c1{1,1}));
        S21_phase = c1{1,1}(5:s*p:size(c1{1,1}));
        S12_mag = c1{1,1}(6:s*p:size(c1{1,1}));
        S12_phase = c1{1,1}(7:s*p:size(c1{1,1}));
        S22_mag = c1{1,1}(8:s*p:size(c1{1,1}));
        S22_phase = c1{1,1}(9:s*p:size(c1{1,1}));
        S11_phase=deg2rad(S11_phase);% degree to rad
        S21_phase=deg2rad(S21_phase);% degree to rad
        S12_phase=deg2rad(S12_phase);% degree to rad
        S22_phase=deg2rad(S22_phase);% degree to rad
    case 'ri'
        S11_re = c1{1,1}(2:s*p:size(c1{1,1}));
        S11_img = c1{1,1}(3:s*p:size(c1{1,1}));
        S21_re = c1{1,1}(4:s*p:size(c1{1,1}));
        S21_img = c1{1,1}(5:s*p:size(c1{1,1}));
        S12_re = c1{1,1}(6:s*p:size(c1{1,1}));
        S12_img = c1{1,1}(7:s*p:size(c1{1,1}));
        S22_re = c1{1,1}(8:s*p:size(c1{1,1}));
        S22_img = c1{1,1}(9:s*p:size(c1{1,1}));
        S11_mag=abs(S11_re+1i.*S11_img);
        S11_phase=angle(S11_re+1i.*S11_img);
        S21_mag=abs(S21_re+1i.*S21_img);
        S21_phase=angle(S21_re+1i.*S21_img);
        S12_mag=abs(S12_re+1i.*S12_img);
        S12_phase=angle(S12_re+1i.*S12_img);
        S22_mag=abs(S22_re+1i.*S22_img);
        S22_phase=angle(S22_re+1i.*S22_img);
    end
end
