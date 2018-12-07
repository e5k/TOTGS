% Program for the determination of the total grainsize distribution of tephra-fall deposits
%
% Usage: >> TOTGS
%
% Code: TOTGS
% By: Costanza Bonadonna (SOEST, University of Hawaii) and Giacomo Marani (Autonomous Systems Laboratory, University of Hawaii)
% Copyright (C) 2004 C. Bonadonna and G. Marani
%
% Data in the examples is from:
% ﻿ Alfano, F., Bonadonna, C., Watt, S., Connor, C., Volentik, A., Pyle, D.M., 2016.
%       Reconstruction of total grain size distribution of the climactic phase
%       of a long-lasting eruption: the example of the 2008-2013 Chaiten
%       eruption. Bull. Volcanol. 78, 46. doi:10.1007/s00445-016-1040-5
%
% Updates by Sebastien Biass (Departement of Earth Sciences, University of Geneva, Switzerland, 2012)
% Novembre 2015:    version 2
% January 2016:     Added parameters of Inman (1952), minor bug fixes 
% July 2016:        The delimitation of the 0-mass contour is now interactive
% October 2017:     Updated plotting functions
% Novemeber 2018:   Replaced dependencies for conversion between ll-utm
%                   Deleted plot_google_map
% December 2018:    Fixed bug in the computation of the CDF
%
% Email contact: costanza.bonadonna _AT_ unige.ch, sbiasse _AT_ hawaii.edu
%
% This program is free software; 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation. 
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%
% This program uses the following script:
%   ll2utm and utm2ll by Francois Beauducel (https://www.mathworks.com/matlabcentral/fileexchange/45699-ll2utm-and-utm2ll)
%

function TOTGS
global t

addpath(genpath('Dependencies/'));

scr = get(0,'ScreenSize');
w   = 400;
h   = 300;

t.fig = figure(...
    'position', [scr(3)/2-w/2 scr(4)/2-h/2 w h],...
    'Color', [.25 .25 .25],...
    'Resize', 'off',...
    'Toolbar', 'none',...
    'Menubar', 'none',...
    'Name', 'TOTGS v2.1',...
    'NumberTitle', 'off');

        t.main = uipanel(...
            'parent', t.fig,...
            'units', 'normalized',...
            'position', [.025 .025 .95 .95],...
            'title', '',...
            'BackgroundColor', [.25 .25 .25],...
            'ForegroundColor', [.9 .5 0],...
            'HighlightColor', [.9 .5 0],...
            'BorderType', 'line',...
            'Title', 'TOTGS');    

        t.coorP = uipanel(...
            'units', 'normalized',...
            'position', [.05 .6750 .425 .25],...
            'title', 'Coordinates',...
            'BackgroundColor', [.25 .25 .25],...
            'ForegroundColor', [.5 .5 .5],...
            'HighlightColor', [.5 .5 .5],...
            'BorderType', 'line');
          
                t.coor_type = uibuttongroup(...
                    'parent', t.coorP,...
                    'units', 'normalized',...
                    'position', [.05 .05 .9 .9],...
                    'BackgroundColor', [.25 .25 .25],...
                    'ForegroundColor', [.25 .25 .25],...
                    'HighlightColor', [.25 .25 .25],...
                    'BorderType', 'line');
                    
                        s1 = sprintf('Requires columns organized as:\n1. Lat (decimal degrees)\n2. Lon (decimal degrees)');
                        t.coor_ll = uicontrol(...
                           'style', 'radiobutton',...
                           'parent', t.coor_type,...
                           'units', 'normalized',...
                           'position', [.0 .5 .95 .45],...
                           'string', 'Lat/Lon',...
                           'ForegroundColor', [1 1 1],...
                           'BackgroundColor', [.25 .25 .25],...
                           'TooltipString', s1);
                        
                       s2 = sprintf('Requires columns organized as:\n1. Easting (m)\n2. Northing (m)\n3. Zone (eg 19 M)');
                        t.coor_utm = uicontrol(...
                           'style', 'radiobutton',...
                           'parent', t.coor_type,...
                           'units', 'normalized',...
                           'position', [.0 .0 .95 .45],...
                           'string', 'UTM',...
                           'ForegroundColor', [1 1 1],...
                           'BackgroundColor', [.25 .25 .25],...
                           'TooltipString', s2);  

    
        t.diamP = uipanel(...
            'units', 'normalized',...
            'position', [.525 .6750 .425 .25],...
            'title', 'Diameter',...
            'BackgroundColor', [.25 .25 .25],...
            'ForegroundColor', [.5 .5 .5],...
            'HighlightColor', [.5 .5 .5],...
            'BorderType', 'line');
          
                t.diam_type = uibuttongroup(...
                    'parent', t.diamP,...
                    'units', 'normalized',...
                    'position', [.05 .05 .9 .9],...
                    'BackgroundColor', [.25 .25 .25],...
                    'ForegroundColor', [.25 .25 .25],...
                    'HighlightColor', [.25 .25 .25],...
                    'BorderType', 'line');
                    
                        s1 = sprintf('Diameter in millimeters\nColumns should be order from fine to coarse (i.e. increasing mm)');
                        t.diam_mm = uicontrol(...
                           'style', 'radiobutton',...
                           'parent', t.diam_type,...
                           'units', 'normalized',...
                           'position', [.0 .5 .95 .45],...
                           'string', 'Millimetre',...
                           'ForegroundColor', [1 1 1],...
                           'BackgroundColor', [.25 .25 .25],...
                           'TooltipString', s1);
                        
                       s2 = sprintf('Diameter in phi units\nColumns should be order from coarse to fine (i.e. increasing phi units)');
                        t.diam_phi = uicontrol(...
                           'style', 'radiobutton',...
                           'parent', t.diam_type,...
                           'units', 'normalized',...
                           'position', [.0 .0 .95 .45],...
                           'string', 'Phi',...
                           'ForegroundColor', [1 1 1],...
                           'BackgroundColor', [.25 .25 .25],...
                           'TooltipString', s2);  

    
        t.fileP = uipanel(...
            'units', 'normalized',...
            'position', [.05 0.375 .9 .25],...
            'title', 'Input file',...
            'BackgroundColor', [.25 .25 .25],...
            'ForegroundColor', [.5 .5 .5],...
            'HighlightColor', [.5 .5 .5],...
            'BorderType', 'line');
        

               t.inp_p = uicontrol(...
                    'parent', t.fileP,...
                    'Style', 'pushbutton',...
                    'units', 'normalized',...
                    'position', [.05 .2 .3 .6],...
                    'BackgroundColor', [.3 .3 .3],...
                    'ForegroundColor', [.9 .5 .0],...
                    'String', 'Load file',...
                    'TooltipString', 'Supported formats: .txt, .xls, .xlsx');
                             
               t.inp_e = uicontrol(...
                    'parent', t.fileP,...
                    'style', 'edit',...
                    'unit', 'normalized',...
                    'position', [.4 .2 .55 .6],...
                    'HorizontalAlignment', 'right',...
                    'ForegroundColor', [1 1 1],...
                    'BackgroundColor', [.35 .35 .35]);

    
        t.zeroP = uipanel(...
            'units', 'normalized',...
            'position', [.05 .075 .425 .25],...
            'title', 'Zero contour',...
            'BackgroundColor', [.25 .25 .25],...
            'ForegroundColor', [.5 .5 .5],...
            'HighlightColor', [.5 .5 .5],...
            'BorderType', 'line');
        
            t.zero_p = uicontrol(...
                'parent', t.zeroP,...
                'Style', 'pushbutton',...
                'units', 'normalized',...
                'position', [.1 .2 .8 .6],...
                'BackgroundColor', [.3 .3 .3],...
                'ForegroundColor', [.9 .5 .0],...
                'String', 'Add',...
                'Enable', 'off',...
                'TooltipString', 'Interactively define a line of 0 deposit from a map');
                
                t.ok_p = uicontrol(...
                    'parent', t.main,...
                    'Style', 'pushbutton',...
                    'units', 'normalized',...
                    'position', [.525 .05 .425 .25],...
                    'BackgroundColor', [.3 .3 .3],...
                    'ForegroundColor', [.9 .5 .0],...
                    'String', 'Ok',...
                    'Enable', 'off',...
                    'TooltipString', 'Go for it!');

set(t.inp_p, 'callback', @txtread)
set(t.zero_p, 'callback', @add_zero) 
set(t.ok_p, 'callback', @Voronoi_TOTGS) 

% Read input file     
function txtread(~, ~)
global t tr

tr = struct;                                            % Main storage structure

%% Get input file
[fl,pt] = uigetfile({'*.txt;*.xls;*.xlsx', 'Compatible files';...
    '*.txt', 'Tab-delimited text files (*.txt)';...
    '*.xls;*.xlsx', 'Excel files'});
if isequal(fl,0)
    return
end
txt2read = [pt,fl];                                     % Path to the file to read

%% Read input file
% If text file
if strcmp(txt2read(end-3:end), '.txt')

    raw         = importdata(txt2read);                 % Load file

    % If geographic coordinates
    if get(t.coor_ll, 'Value') == 1 
        if size(raw.textdata, 2) ~= 3                  % Check file
             errordlg('The input file does not fit the coordinate system. Please double check it!');
             return
        end
        lon     = str2double(raw.textdata(2:end,2));    % Get longitude
        lat     = str2double(raw.textdata(2:end,1));    % Get latitude
        mass    = str2double(raw.textdata(2:end,3));    % Get mass
    % If projected coordinates
    else   
        if size(raw.textdata, 2) ~= 4                  % Check file
             errordlg('The input file does not fit the coordinate system. Please double check it!');
             return
        end
        east    = str2double(raw.textdata(2:end,1));    % Get easting
        nort    = str2double(raw.textdata(2:end,2));    % Get get northing
        zone    = str2double(raw.textdata(2:end,3));                % Get zone 
        mass    = str2double(raw.textdata(2:end,4));    % Get mass          
    end
    
    gClass      = raw.data(1,2:end);
    gwt         = raw.data(2:end,2:end);

% If excel file
else 
    [~, ~, raw] = xlsread(txt2read);                    % Load file
    
    % If geographic coordinates
    if get(t.coor_ll, 'Value') == 1
        if ~isnan(raw{1,4})                            % Check file
             errordlg('The input file does not fit the coordinate system. Please double check it!');
             return  
        end 
        lon     = cell2mat(raw(2:end,2));               % Get longitude
        lat     = cell2mat(raw(2:end,1));               % Get latitude
        mass    = cell2mat(raw(2:end,3));               % Get mass
        % Get data and order it by ascending bins (i.e. in mm: fine -> coarse; in phi: coarse -> fine)
        tmp     = sortrows([cell2mat(raw(1,5:end)); cell2mat(raw(2:end,5:end))]',1)';          
        gClass  = tmp(1,:);                             % Get bins
        gwt     = tmp(2:end,:);                         % Get grainsize data      
        
    % If projected coordinates
    else
        if ~isnan(raw{1,5})                            % Check file
             errordlg('The input file does not fit the coordinate system. Please double check it!');
             return  
        end 
        east    = cell2mat(raw(2:end,1));               % Get easting 
        nort    = cell2mat(raw(2:end,2));               % Get northing
        zone    = cell2mat(raw(2:end,3));               % Get zone
        mass    = cell2mat(raw(2:end,4));               % Get mass
        % Get data and order it by ascending bins (i.e. in mm: fine -> coarse; in phi: coarse -> fine)
        tmp     = sortrows([cell2mat(raw(1,6:end)); cell2mat(raw(2:end,6:end))]',1)';  
        gClass  = tmp(1,:);                             % Get bins
        gwt     = tmp(2:end,:);                         % Get grainsize data     
    end
end


%% Process data
% In case some Nan were read, remove them
idxNan = ~isnan(mass);
mass    = mass(idxNan,:);
gwt     = gwt(idxNan,:);
% Flip matrices if in PHI units
if get(t.diam_phi, 'Value') == 1
    gwt     = fliplr(gwt);
    gClass  = fliplr(gClass);
end

% Check that bins are equally-spaced in phi units
if get(t.diam_mm, 'Value') == 1
    p    = -log2(gClass);
    pdif = p(1:end-1)-p(2:end);
else
    pdif = gClass(1:end-1)-gClass(2:end);
end

if mean([abs(mean(pdif)-min(pdif)), abs(max(pdif)-mean(pdif))])*100 > 10
    choice = questdlg(sprintf('It seems that the bin spacing varies by more than 10% (Max:  %2.2f)\n Do you want to continue?',max([abs(mean(pdif)-min(pdif)), abs(max(pdif)-mean(pdif))])*100),...
        '', 'No', 'Yes', 'Yes');
    switch choice
        case 'No'
            return
        case 'Yes'
    end
end

% Check there are no missing data
if ~isempty(gwt(isnan(gwt)))
    choice = questdlg(sprintf('It seems there are missing data. Would you like to feel them with 0?'),...
        '', 'No', 'Yes', 'Yes');
    switch choice
        case 'No'
            return
        case 'Yes'
            gwt(isnan(gwt)) = 0;
    end
end

% Check that sum of every row is equal to 100
nDigits  = 2;
for i = 1:size(mass,1)
    kSum = sum(gwt(i,:));
    if round(kSum*10^nDigits)/10^nDigits ~= 100 %|| round(kSum*10^nDigits)/10^nDigits > 0
        choice = questdlg(sprintf('Warning: sum of gwt on row %i is not 100 (actual value: %.02f).\n Would you like to continue?', i, kSum),...
            '', 'No', 'Yes', 'Yes to all', 'Yes');
        switch choice
            case 'No'
                return
            case 'Yes'
                continue
            case 'Yes to all'
                break
        end
    end
end
        

% Creates and corrects coordinates
% Coordinates in Lat/Lon
if get(t.coor_ll, 'Value') == 1   
    % Remove nan
    lat     = lat(idxNan, :);
    lon     = lon(idxNan, :);
    [east, nort, zone]  = ll2utm(lat,lon);
   
% Coordinates in UTM  
else
    % Remove nan
    east        = east(idxNan, :);
    nort        = nort(idxNan, :);
    zone        = zone(idxNan, :);   
    [lat, lon]  = utm2ll(east,nort,zone);
end

% Check if different zones
if length(unique(zone)) > 1
    zns = unique(zone);
    lst = cellstr(num2str(zns)); 
    indx = listdlg('PromptString','Choose a reference UTM zone:',...
                           'SelectionMode','single',...
                           'ListString',lst);
    tr.ref_zone = zns(indx);
    zone = ones(size(zone)).*tr.ref_zone;
    [east, nort] = ll2utm(lat,lon,zone);
else
    tr.ref_zone = unique(zone);
end


%% Save data to structure
tr.east  = east;
tr.nort  = nort;
tr.zone  = zone;
tr.lat   = lat;
tr.lon   = lon;
tr.m     = mass;
tr.gClass= gClass;
tr.gwt   = gwt;

% Prepare GS binning for later display
if get(t.diam_mm, 'Value') == 1     % If input in mm
    tr.pClass = -log2(tr.gClass');
    tr.mClass = tr.gClass';
else                                % If input in phi
    tr.pClass = gClass';
    tr.mClass = 2.^-gClass';
end

set(t.inp_e, 'String', txt2read); 
set(t.zero_p, 'Enable', 'on');
set(t.ok_p, 'Enable', 'on');

% Interactive map to add the zero line
function add_zero(~, ~)                        
global t tr map

if isempty(get(t.inp_e, 'String'))
    errordlg('Please first load a file');
    return
end

scr = get(0,'ScreenSize');
h   = 500;
w   = 600;
map.fig = figure(...
    'position', [scr(3)/2-w/2 scr(4)/2-h/2 w h],...
    'Resize', 'off',...
    'Name', 'Plot zero line',...
    'NumberTitle', 'off',...
    'KeyPressFcn', @process_zero);

map.add = uicontrol(...
    'parent', map.fig,...
    'Style', 'togglebutton',...
    'units', 'normalized',...
    'position', [.1 .03 .25 .08],...
    'String', 'Add points',...
    'Callback', @process_zero);

map.reset = uicontrol(...
    'parent', map.fig,...
    'Style', 'pushbutton',...
    'units', 'normalized',...
    'position', [.375 .03 .25 .08],...
    'String', 'Reset',...
    'Enable', 'off',...
    'Callback', @process_zero);

map.go = uicontrol(...
    'parent', map.fig,...
    'Style', 'pushbutton',...
    'units', 'normalized',...
    'position', [.65 .03 .25 .08],...
    'String', 'Ok',...
    'Enable', 'off',...
    'Callback', @process_zero);

map.a  = axes('Parent', map.fig, 'Position', [.1 .2 .8 .75], 'Units', 'normalized', 'Box', 'on');
map.p  = plot(map.a, tr.lon, tr.lat, '.r','MarkerSize', 15); axis([min(tr.lon)-10 max(tr.lon+10) min(tr.lat)-10 max(tr.lat)+10]);
xlabel('Longitude'); ylabel('Latitude');
%plot_google_map('Maptype', 'terrain')
set(map.a,'layer','top')
hold on

% Handle results of the zero-line map
function process_zero(h,~)
global map tr

%% Add point
if strcmp(get(h, 'String'), 'Add points') && get(h, 'Value') == 1
    
    warndlg(sprintf('Central mouse click to cancel the last point\nRight click to exit'));
    uiwait
    
    set(map.reset, 'Enable', 'off'); set(map.go, 'Enable', 'off'); set(map.add, 'String', 'Right click to exit!');
    hold on
    while get(h, 'Value') == 1        
        [x,y,k] = ginput(1);
        plot_voronoi(x,y,k);
    end

%% Reset point
elseif strcmp(get(h, 'String'), 'Reset')
    set(map.add, 'Enable', 'on');  
    delete(map.h1);
    delete(map.h2);
%% Get data
elseif strcmp(get(h, 'String'), 'Ok')
    zpoint     = findobj(map.h1, 'Marker', 'x');

    lt  = get(zpoint, 'YData')';
    ln  = get(zpoint, 'XData')';
    [east, nort,z] = ll2utm(lt,ln,tr.ref_zone);

    tr.idx   = [ones(size(tr.east)); zeros(size(lt))]; % Added index for plotting
    tr.east  = [tr.east; east];
    tr.nort  = [tr.nort; nort];
    tr.zone  = [tr.zone; z];
    tr.lat   = [tr.lat; lt];
    tr.lon   = [tr.lon; ln];
    tr.m     = [tr.m; zeros(length(lt),1)];
    tr.gwt   = [tr.gwt; zeros(length(lt), length(tr.gClass))];
    
    close(map.fig);
end
  
% Voronoi functions             
function Voronoi_TOTGS(~, ~)
global tr

% Voronoi
[v,c]    = voronoin([tr.east,tr.nort]);
nPoints  = size(tr.east,1);
nClass   = size(tr.gClass',1);
% Mass computation
PolyMass = zeros(nPoints,1);
gm       = zeros(nPoints, nClass);
gmTot    = zeros(nClass,1);
vorWt    = zeros(nClass,1);

for i = 1:size(c,1)
    PolyMass(i) = tr.m(i) * polyarea(v(c{i}',1), v(c{i}',2));
    if isnan(PolyMass(i))
        PolyMass(i) = 0;
    end
    
    for j=1:nClass
        gm(i,j) = tr.gwt(i,j)*PolyMass(i)/100;
    end
end


for i=1:nClass
    gmTot(i) = sum(gm(:,i));
end

TotalGran = sum(gmTot);

for i=1:nClass
    vorWt(i) = gmTot(i)/TotalGran*100;
end

% Plot voronoi cells
[vg,cg]  = voronoin([tr.lon,tr.lat]);                           % Recalculate Voronoi with geographic coordinates for plotting
figure;
voronoi(tr.lon,tr.lat);                                         % Plot Voronoi cells
axis([min(tr.lon)-2, max(tr.lon)+2, min(tr.lat)-2, max(tr.lat)+2])
hold on
% Plot patch with color as total mass in the cell
for i = 1:length(cg)
    if all(cg{i}~=1) && PolyMass(i) ~= 0
        %patch(vg(cg{i},1), vg(cg{i},2), log10(PolyMass(i)/1000));
        patch(vg(cg{i},1), vg(cg{i},2), polyarea(v(c{i}',1), v(c{i}',2)));
    end
end
alpha .7        % Transparency
c = colorbar;   % Colorbar
ylabel(c, 'Log10 Mass (kg in the cell)');
% Plot points
if isfield(tr, 'idx')
    plot(tr.lon(tr.idx==1),tr.lat(tr.idx==1), '.r')
    plot(tr.lon(tr.idx==0),tr.lat(tr.idx==0), 'ok', 'MarkerFaceColor', 'm', 'MarkerSize',5)
else
    plot(tr.lon,tr.lat, '.r')
end

xlabel('Longitude');
ylabel('Latitude');
%plot_google_map('maptype', 'terrain')
set(gca,'layer','top')

res_voron(vorWt);

% Calculations and results                        
function res_voron(vorWt) 
global tr tmp r

% Cumulative


cum = ones(size(tr.pClass));
for i = 2:length(cum)
    cum(i) = 1-sum(vorWt(1:i-1))/sum(vorWt);
end

% Prepare result table
tmp = [tr.mClass, tr.pClass, vorWt, cum];
data_table = cell(size(tmp,1),4);
for i = 1:size(tmp,1)
    for j = 1:4
        data_table{i,j} = sprintf('%.3f', tmp(i,j));
    end
end


% Calculate parameters according to Table 1 of Inman (1952)
% Percentile values are interpolated, not fitted
gs_50       = interp1(cum, tr.pClass, .50);       % Median GS
gs_5        = interp1(cum, tr.pClass, .05);       % 5th percentile
gs_16       = interp1(cum, tr.pClass, .16);       % 16th percentile
gs_84       = interp1(cum, tr.pClass, .84);       % 84th percentile
gs_95       = interp1(cum, tr.pClass, .95);       % 95th percentile

gs_std_new  = (gs_84 - gs_16)/2;               % std
gs_mean     = (gs_16 + gs_84)/2;               % Mean
gs_skw      = (gs_mean - gs_50)/gs_std_new;    % Skewness
gs_kurt     = (.5*(gs_95-gs_5)-gs_std_new)/gs_std_new; % Kurtosis

% Prepare summary table
sum_table       = cell(11,2);
sum_table(1,:)  = {sprintf('%.3f', gs_50), sprintf('%.4f', 2^-gs_50)};
sum_table(2,:)  = {sprintf('%.3f', gs_std_new), sprintf('%.4f', 2^-gs_std_new)};
sum_table(3,:)  = {'', ''};
sum_table(4,:)  = {sprintf('%.3f', gs_mean), sprintf('%.4f', 2^-gs_mean)};
sum_table(5,:)  = {sprintf('%.3f', gs_skw), sprintf('%.4f', 2^-gs_skw)};
sum_table(6,:)  = {sprintf('%.3f', gs_kurt), sprintf('%.4f', 2^-gs_kurt)};
sum_table(7,:)  = {'', ''};
sum_table(8,:)  = {sprintf('%.3f', gs_5), sprintf('%.4f', 2^-gs_5)};
sum_table(9,:)  = {sprintf('%.3f', gs_16), sprintf('%.4f', 2^-gs_16)};
sum_table(10,:) = {sprintf('%.3f', gs_84), sprintf('%.4f', 2^-gs_84)};
sum_table(11,:) = {sprintf('%.3f', gs_95), sprintf('%.4f', 2^-gs_95)};


scr = get(0,'ScreenSize');
w   = 750;
h   = 550;

r.fig = figure(...
    'position', [scr(3)/2-w/2 scr(4)/2-h/2 w h],...
    'Color', [.25 .25 .25],...
    'Resize', 'off',...
    'Toolbar', 'none',...
    'Menubar', 'none',...
    'Name', 'Results',...
    'NumberTitle', 'off');

        r.main = uipanel(...
            'parent', r.fig,...
            'units', 'normalized',...
            'position', [.025 .025 .95 .95],...
            'title', '',...
            'BackgroundColor', [.25 .25 .25],...
            'ForegroundColor', [.9 .5 0],...
            'HighlightColor', [.9 .5 0],...
            'BorderType', 'line'); 
        
        r.tab = uitable(...
            'parent', r.main,...
            'units', 'normalized',...
            'position', [.025 .02 .26 .95],...
            'BackgroundColor', [.3 .3 .3; .2 .2 .2],...
            'ColumnName', {'mm', 'phi', 'PDF', 'CDF'},...
            'ColumnWidth', {40, 40, 40, 40},...
            'ColumnFormat', {'char', 'char', 'char', 'char'},...
            'ColumnEditable', [false false false false],...
            'Data', data_table,...
            'RowStriping', 'on',...
            'RowName', [],...
            'ForegroundColor', [1 1 1]);
        
        r.ax = axes(...
            'parent', r.main,...
            'position', [.365 .425 .55 .525],...
            'units', 'normalized',...
            'XColor', [1 1 1],...
            'YColor', [1 1 1],...
            'FontSize', 8);  
        
        r.tab2 = uitable(...
            'parent', r.main,...
            'units', 'normalized',...
            'position', [.37 .125 .55 .225],...
            'BackgroundColor', [.3 .3 .3; .2 .2 .2],...
            'ColumnName', {'phi', 'mm'},...
            'ColumnWidth', {160, 160},...
            'ColumnFormat', {'char', 'char'},...
            'ColumnEditable', [false false],...
            'Data', sum_table,...
            'RowStriping', 'on',...
            'RowName', {'Median', 'Sigma', '', 'Mean', 'Skewness', 'Kurtosis', '', '5th pctl', '16th pctl', '84th pctl', '95th pctl'},...
            'ForegroundColor', [1 1 1]);       
                
       r.close_p = uicontrol(...
            'parent', r.main,...
            'Style', 'pushbutton',...
            'units', 'normalized',...
            'position', [.835 .02 .15 .075],...
            'BackgroundColor', [.3 .3 .3],...
            'ForegroundColor', [.9 .5 .0],...
            'String', 'Close',...
            'TooltipString', 'Close window');
                
       r.exp_cdf_p = uicontrol(...
            'parent', r.main,...
            'Style', 'pushbutton',...
            'units', 'normalized',...
            'position', [.670 .02 .15 .075],...
            'BackgroundColor', [.3 .3 .3],...
            'ForegroundColor', [.9 .5 .0],...
            'String', 'Export CDF',...
            'TooltipString', 'Export the cumulative density function to new plot');
                
       r.exp_pdf_p = uicontrol(...
            'parent', r.main,...
            'Style', 'pushbutton',...
            'units', 'normalized',...
            'position', [.505 .02 .15 .075],...
            'BackgroundColor', [.3 .3 .3],...
            'ForegroundColor', [.9 .5 .0],...
            'String', 'Export PDF',...
            'TooltipString', 'Export the probability density function to new plot');
                
       r.exp_table_p = uicontrol(...
            'parent', r.main,...
            'Style', 'pushbutton',...
            'units', 'normalized',...
            'position', [.340 .02 .15 .075],...
            'BackgroundColor', [.3 .3 .3],...
            'ForegroundColor', [.9 .5 .0],...
            'String', 'Export Table',...
            'TooltipString', 'Export the table to a text file');
              
        
        [AX, H1, H2] = plotyy(tmp(:,2), tmp(:,3), tmp(:,2), tmp(:,4)*100, 'bar', 'plot', 'Parent', r.ax);
        set(get(AX(1), 'Ylabel'), 'String', 'PDF (wt. %)');
        set(get(AX(1), 'Ylabel'), 'FontSize', 8);
        set(get(AX(1), 'Ylabel'), 'FontWeight', 'bold');
        set(get(AX(2), 'Ylabel'), 'String', 'CDF (wt. %)');
        set(get(AX(2), 'Ylabel'), 'FontSize', 8);
        set(get(AX(2), 'Ylabel'), 'FontWeight', 'bold');
        title('TOTGS', 'FontWeight', 'bold', 'Color', [1 1 1]);
        xlabel('Diameter (phi)', 'FontWeight', 'bold', 'Color', [1 1 1]);
        
        set(H1, 'FaceColor', [.8 .8 .8]);
        set(H2, 'Color', [.9 .5 0]);
        set(H2, 'LineStyle', '-');
        set(AX(2), 'XTickLabel', []);
        set(AX(1), 'YColor', [1 1 1]);
        set(AX(1), 'XColor', [1 1 1]);
        set(AX(2), 'YColor', [.9 .5 0]);
        set(AX(2), 'FontSize', 9);

        set(r.exp_table_p, 'callback', @exp_table)
        set(r.exp_cdf_p, 'callback', @exp_cdf)
        set(r.exp_pdf_p, 'callback', @exp_pdf)
        set(r.close_p, 'callback', @p_close)

% Callbacks to the results GUI        
function p_close(~, ~)
    close(gcf)
    
function exp_table(~, ~)
    global tmp
    [fl, pt] = uiputfile('TOTGS.txt', 'Save TOTGS');
    if isequal(fl,0)
       return
    end
    % Write TOTGS file (modified 2013/03/29 after R. Sulpizio'problems)
    fid = fopen([pt, fl], 'w');
    fprintf(fid, '%s\t%s\t%s\t%s\n', 'mm', 'phi', 'PDF (wt%)', 'CDF');
    for i = 1:size(tmp,1)
        fprintf(fid, '%7.3g\t%4.2g\t%6.3f\t%6.3f\n', tmp(i,1), (round(tmp(i,2)*100))/100, tmp(i,3), tmp(i,4));
    end
    fclose(fid);
    
function exp_pdf(~,~)  
    global tmp
    figure
    bar(tmp(:,2), tmp(:,3));
    colormap([.8 .8 .8]);
    xlabel('Diameter (phi)');
    ylabel('wt. %');
    
function exp_cdf(~,~)  
    global tmp
    figure
    plot(tmp(:,2), tmp(:,4), ':.k');
    xlabel('Diameter (phi)');
    ylabel('wt. %');

function [lt, ln] = plot_voronoi(x,y,k)
    global map tr
    
    if isfield(map, 'h1')
        zpoint     = findobj(map.h1, 'Marker', 'x'); % Retrieve added points
        lt  = get(zpoint, 'YData')';
        ln  = get(zpoint, 'XData')';
    else
        lt = [];
        ln = [];
    end
    
    % Add and delete points
    if k < 3
        % If voronoi already exists, delete it
        if isfield(map, 'h2'); delete(map.h2); end
        if isfield(map, 'h1'); delete(map.h1); end

        if k == 1
            ln = [x;ln];
            lt = [y;lt];
            map.h1 = plot(map.a, ln, lt,'xm');
        elseif k == 2
            lt = lt(2:end);
            ln = ln(2:end);
            map.h1 = plot(map.a,ln,lt,'xm');
        end

        % Plot voronoi
        [tX, tY]   = voronoi([tr.lon;ln], [tr.lat;lt]);
        map.h2     = plot(map.a, [tr.lon;ln], [tr.lat;lt], '.m', tX, tY, 'k-', 'LineWidth', .2);

        % Reorder objects
        uistack(map.p, 'top');
        uistack(map.a, 'top');
    
    % Finish adding points
    else
        set(map.reset, 'Enable', 'on'); set(map.go, 'Enable', 'on');
        set(map.add, 'Value', 0, 'String', 'Add points', 'Enable', 'off');
        return
    end


