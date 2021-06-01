function [f,S_re,S_img]=readS(filename, p,s)
%% Read S params as linear form 
    fileID = fopen(filename);
%   exclude the first Engligh legend and Delimiter '---'
    textscan(fileID,'%[^-]');
    textscan(fileID,'%[--]');
    c1 = textscan(fileID,'%f');
    f = c1{1,1}(1:s*p:size(c1{1,1}));
    S_re = c1{1,1}(2:s*p:size(c1{1,1}));
    S_img = c1{1,1}(3:s*p:size(c1{1,1}));
end