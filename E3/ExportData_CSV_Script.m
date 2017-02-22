clear all;
close all;
clc;

datafolder1 = '/Users/user/Desktop/E1';
datafolder2 = '/Users/user/Desktop/E2';
datafolder3 = '/Users/user/Desktop/E3/datafiles';
csvfolder3  = '/Users/user/Desktop/E3/csv';

%cd(datafolder1)
%cd(datafolder2)
cd(datafolder3)

sublist1 = [2 3 5 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28]; %experiment 1
sublist2 = [1 2 3 5 6 7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 24 28 29 30]; %experiment 2
sublist3 = [1]; %experiment 3, within subjects
nsub = length(sublist3);
%%
for isub = 1:length(sublist3)

fprintf('\nWorking on S%02d/%02d.',isub,length(sublist3));
fname = sprintf('%s/S%02d.mat',datafolder3,sublist3(isub));
load(fname,'sequence','response');
iblocks = 1:length(response);


resp      = cat(1,response(iblocks).resp);
responded = all(~isnan(resp) & ~isinf(resp),2); resp(~responded,:) = [];

theta     = cat(1,response(iblocks).theta); theta(~responded,:) = [];
if size(theta,1) < size(resp,1), theta(end+1,:) = 0; end
theta     = resp*pi/180;
time      = cat(1,response(iblocks).time);  time(~responded,:)  = [];

dist      = cat(1,response(iblocks).dist);  time(~responded,:)  = [];
prec      = cat(1,response(iblocks).prec);  prec(~responded,:)  = [];
pstartang = cat(1,response(iblocks).probestartangle); pstartang(~responded,:) = [];


pord     = cat(1,sequence(iblocks).probeorder); pord(~responded,:)  = [];
angs     = cat(1,sequence(iblocks).angle);      angs(~responded,:)  = [];
locs     = cat(1,sequence(iblocks).loc);        locs(~responded,:)  = [];
cues     = cat(1,sequence(iblocks).retrocue);   cues(~responded,:)  = [];
btype    = cat(1,sequence(iblocks).showorder);  btype(~responded,:) = [];


subID    = zeros(length(resp),1);
subID(:) = sublist3(isub);

dat = cat(2, subID,   btype,   resp,   theta,   time,   pord,   angs,   locs,   cues,   dist,   prec,   pstartang);
names =    {'subID', 'btype', 'resp', 'theta', 'time', 'pord', 'angs', 'locs', 'cues', 'dist', 'prec', 'pstartang'};

data = dataset(subID, btype, resp, theta, time, pord, angs, locs, cues, dist, prec, pstartang, 'VarNames', names);
subj = int2str(sublist3(isub));
subject.filename = sprintf('sub%s',subj);
filename = sprintf('%s/DoubleProbe_Within_S%02d.csv', csvfolder3, sublist3(isub));

export(data, 'File', filename, 'delimiter', ',');
clear data
end