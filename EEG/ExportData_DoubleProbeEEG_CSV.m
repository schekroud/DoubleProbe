clear all;
close all;
clc;

datafolder = '/Users/user/Desktop/Experiments/Nick/DoubleProbe/EEG/data';
csvfolder  = '/Users/user/Desktop/Experiments/Nick/DoubleProbe/EEG/data/csv';
cd(datafolder);
if ~exist(csvfolder); mkdir(csvfolder); end
  

sublist = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
nsubs = length(sublist);
%%
for isub = 1:nsubs;
    fprintf('\nWorking on S%02d/%02d.',isub, nsubs);
    fname = sprintf('%s/S%02d.mat',datafolder, sublist(isub));
    load(fname, 'sequence', 'response');
    iblocks = 1:length(response);
    
    %sequence: length, probe, probeorder, loc, colors, angle, rate,
    %retrocue, memdelay, stimecc, locangle, 
    
    %response: probestartangle, resp, dist, prec, corr, time, starttime,
    %theta
    
    resp      = cat(1, response(iblocks).resp);
    responded = all(~isnan(resp) & ~isinf(resp),2)           ; resp(~responded,:)  = [];
    
    theta     = cat(1,response(iblocks).theta)               ; theta(~responded,:) = [];
    if size(theta,1) < size(resp,1); theta(end+1,:) = 0; end
    theta     = resp*pi/180;
    time = cat(1,response(iblocks).time)                     ; time(~responded,:) = [];
    
    dist = cat(1,response(iblocks).dist)                     ; dist(~responded,:) = [];
    prec = cat(1,response(iblocks).prec)                     ; prec(~responded,:) = [];
    pstartang = cat(1,response(iblocks).probestartangle)     ; pstartang(~responded,:) = [];
    
    pord = cat(1,sequence(iblocks).probeorder)               ; pord(~responded,:) = [];
    angs = cat(1,sequence(iblocks).angle)                    ; angs(~responded,:) = [];
    locs = cat(1,sequence(iblocks).loc)                      ; locs(~responded,:) = [];
    cues = cat(1,sequence(iblocks).retrocue)                 ; cues(~responded,:) = [];
    pids = cat(1,response(iblocks).probeitems); pids(~responded,:) = [];
    prepeat = pids(:,1) == pids(:,2);
    
    
    subID    = zeros(length(resp),1);
    subID(:) = sublist(isub);
    
    dat = cat(2, subID,   resp,   theta,   time,   pord,   angs,   locs,   cues,   dist,   prec,   pstartang);
    names = {   'subID', 'resp', 'theta', 'time', 'pord', 'angs', 'locs', 'cues', 'dist', 'prec', 'pstartang'};
    
    data = dataset(subID, resp, theta, time, pord, angs, locs, cues, dist, prec, pstartang, 'VarNames', names);
    filename = sprintf('%s/DoubleProbeEEG_S%02d.csv', csvfolder, sublist(isub));
    export(data, 'File', filename, 'delimiter', ',');
    clear data;
end

    
    
    