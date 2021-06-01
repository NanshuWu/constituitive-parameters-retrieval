function [f,S11_re,S11_img,S21_re,S21_img,S12_re,S12_img,S22_re,S22_img]=readS_fromtest(filename, p,s)
%% Read S params as linear form 
    fileID = fopen(filename);

    c1 = textscan(fileID,'%f','headerlines',8);
    f = c1{1,1}(1:s*p:size(c1{1,1}));
    S11_re = c1{1,1}(2:s*p:size(c1{1,1}));
    S11_img = c1{1,1}(3:s*p:size(c1{1,1}));
    S21_re = c1{1,1}(4:s*p:size(c1{1,1}));
    S21_img = c1{1,1}(5:s*p:size(c1{1,1}));
    S12_re = c1{1,1}(6:s*p:size(c1{1,1}));
    S12_img = c1{1,1}(7:s*p:size(c1{1,1}));
	S22_re = c1{1,1}(8:s*p:size(c1{1,1}));
    S22_img = c1{1,1}(9:s*p:size(c1{1,1}));
end
% function [ f, S11, S21 ] = readSparams( filename )
%     fid = fopen(filename);
%     c = textscan(fid,'%f %f %f','headerlines',2);
%     f = c{1,1};
%     S11 = (c{1,2} .* exp(1i*c{1,3}*pi/180));
%     c = textscan(fid,'%f %f %f','headerlines',2);
%     S21 = c{1,2} .* exp(1i*c{1,3}*pi/180);
%     fclose(fid);
% end