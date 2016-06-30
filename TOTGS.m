% Program for the determination of the total grainsize distribution of tephra-fall deposits
%
% Usage: >> TOTGS
%
% Code: TOTGS
% By: Costanza Bonadonna (SOEST, University of Hawaii) and Giacomo Marani (Autonomous Systems Laboratory, University of Hawaii)
% Copyright (C) 2004 C. Bonadonna and G. Marani
%
% Updates by Sebastien Biass (Departement of Earth Sciences, University of Geneva, Switzerland, 2012)
% Novembre 2015: version 2
%                           
%
% Email contact: costanza.bonadonna@unige.ch, sebastien.biasse@unige.ch
%
% This program is free software; 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation. 
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%
% This program uses the following script:
% utm2deg           by Rafael Palacios (http://www.mathworks.com/matlabcentral/fileexchange/10914-utm2deg)
% deg2utm           by Rafael Palacios (http://www.mathworks.com/matlabcentral/fileexchange/10915-deg2utm)
% plot_google_map   by Zohar Bar-Yehunda (http://www.mathworks.com/matlabcentral/fileexchange/27627-plotgooglemap)

function TOTGS
global t

scr = get(0,'ScreenSize');
w   = 400;
h   = 300;

t.fig = figure(...
    'position', [scr(3)/2-w/2 scr(4)/2-h/2 w h],...
    'Color', [.25 .25 .25],...
    'Resize', 'off',...
    'Toolbar', 'none',...
    'Menubar', 'none',...
    'Name', 'TOTGS v2',...
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
[fl,pt] = uigetfile({'*.txt'; '*.xls'; '*.xlsx'});
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
        if size(raw.textdata, 2) ~= 3;                  % Check file
             errordlg('The input file does not fit the coordinate system. Please double check it!');
             return
        end
        lon     = str2double(raw.textdata(2:end,2));    % Get longitude
        lat     = str2double(raw.textdata(2:end,1));    % Get latitude
        mass    = str2double(raw.textdata(2:end,3));    % Get mass
    % If projected coordinates
    else   
        if size(raw.textdata, 2) ~= 4;                  % Check file
             errordlg('The input file does not fit the coordinate system. Please double check it!');
             return
        end
        east    = str2double(raw.textdata(2:end,1));    % Get easting
        nort    = str2double(raw.textdata(2:end,2));    % Get get northing
        zone    = raw.textdata(2:end,3);                % Get zone 
        mass    = str2double(raw.textdata(2:end,4));    % Get mass          
    end
    
    gClass      = raw.data(1,2:end);
    gwt         = raw.data(2:end,2:end);

% If excel file
else 
    [~, ~, raw] = xlsread(txt2read);                    % Load file
    
    % If geographic coordinates
    if get(t.coor_ll, 'Value') == 1
        if ~isnan(raw{1,4});                            % Check file
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
        if ~isnan(raw{1,5});                            % Check file
             errordlg('The input file does not fit the coordinate system. Please double check it!');
             return  
        end 
        east    = cell2mat(raw(2:end,1));               % Get easting
        nort    = cell2mat(raw(2:end,2));               % Get northing
        zone    = raw(2:end,3);                         % Get zone
        mass    = cell2mat(raw(2:end,4));               % Get mass
        % Get data and order it by ascending bins (i.e. in mm: fine -> coarse; in phi: coarse -> fine)
        tmp     = sortrows([cell2mat(raw(1,6:end)); cell2mat(raw(2:end,6:end))]',1)';  
        gClass  = tmp(1,:);                             % Get bins
        gwt     = tmp(2:end,:);                         % Get grainsize data     
    end
end


%% Process data
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
    [east, nort, zone]  = deg2utm(lat, lon);
    %zone                = cellstr(zone);    
% Coordinates in UTM  
else   
    [lat, lon]          = utm2deg(east,nort,char(zone));
end

% Check if different zones
if length(unique(cellstr(zone))) > 1
   dummy            = unique(cellstr(zone));
   dummy1           = sscanf(dummy{1}, '%d');
   dummy2           = sscanf(dummy{2}, '%d');

   % Check if zones vary laterally
   if dummy1 ~= dummy2
       choice = questdlg('It seems your data crosses over several UTM zones. Which one would you like to use as a reference projection?', ...
            'UTM', ...
            dummy{1},dummy{2},dummy{1});
       switch choice
            case dummy{1}
                ref_zone = dummy{1};		
            case dummy{2}
                ref_zone = dummy{2};
       end
       [east, zone] = correct_utm(east, zone, ref_zone, lat);
   end  
end


%% Save data to structure
tr.east  = east;
tr.nort  = nort;
tr.zone  = zone;
tr.lat   = lat;
tr.lon   = lon;
if isfield(tr, 'ref_zone');
    tr.ref_zone = ref_zone;
end
tr.m     = mass;
tr.gClass= gClass;
tr.gwt   = gwt;

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
plot(map.a, tr.lon, tr.lat, '.r','MarkerSize', 15); axis([min(tr.lon)-10 max(tr.lon+10) min(tr.lat)-10 max(tr.lat)+10]);
xlabel('Longitude'); ylabel('Latitude');
plot_google_map
hold on

% Handle results of the zero-line map
function process_zero(h,~)
global map tr

%% Add point
if strcmp(get(h, 'String'), 'Add points') && get(h, 'Value') == 1
    set(map.reset, 'Enable', 'off'); set(map.go, 'Enable', 'off'); set(map.add, 'String', 'Right click to exit!');
    hold on
    while get(h, 'Value') == 1        
        [x,y,k] = ginput(1);
        if k == 1
            map.h1      = plot(map.a,x,y,'xr');
        else
            zpoint     = findobj(map.a, 'Marker', 'x');
            lt         = zeros(size(zpoint));
            ln         = zeros(size(zpoint));
            for i = 1:length(zpoint)
                lt(i)  = get(zpoint(i), 'YData');
                ln(i)  = get(zpoint(i), 'XData');
            end
            map.h2     = voronoi(map.a, [tr.lon;ln], [tr.lat;lt]);
            set(map.reset, 'Enable', 'on'); set(map.go, 'Enable', 'on'); 
            set(map.add, 'Value', 0, 'String', 'Add points', 'Enable', 'off');  
            return
        end
    end

%% Reset point
elseif strcmp(get(h, 'String'), 'Reset')
    set(map.add, 'Enable', 'on');  
    delete(findobj(map.a, 'Marker', 'x'));
    delete(map.h2);
%% Get data
elseif strcmp(get(h, 'String'), 'Ok')
    zpoint     = findobj(map.a, 'Marker', 'x');
    lt         = zeros(size(zpoint));
    ln         = zeros(size(zpoint));
    for i = 1:length(zpoint)
        lt(i)  = get(zpoint(i), 'YData');
        ln(i)  = get(zpoint(i), 'XData');
    end
    [dummyE, dummyN, dummyZ] = deg2utm(lt, ln);
    dummyZ  = num2cell(dummyZ,2);
    % If zero contour is over more than 1 §zone
    if length(unique(dummyZ)) > 1
        dummy =  unique(dummyZ);
        dummy1= sscanf(dummy{1}, '%d');
        dummy2= sscanf(dummy{2}, '%d');
        % Check if zones vary laterally
        if dummy1 ~= dummy2 && ~isfield(tr, 'ref_zone')
            choice = questdlg('It seems your data crosses over several UTM zones. Which one would you like to use as a reference projection?', ...
                'UTM', ...
                dummy{1},dummy{2},dummy{1});
            switch choice
                case dummy{1}
                    ref_zone = dummy{1};
                case dummy{2}
                    ref_zone = dummy{2};
            end
            [dummyE, dummyZ] = correct_utm(dummyE, dummyZ, ref_zone, lt);
        end
        % If zero contour is over one zone, but which is different from
        % reference zone§
    elseif length(unique(dummyZ)) == 1 && isfield(tr, 'ref_zone') && ~strcmp(unique(dummyZ), tr.ref_zone)
        [dummyE, dummyZ] = correct_utm(dummyE, dummyZ, tr.ref_zone, lt);
    end
    
    % Update structure
    if exist('ref_zone', 'var')
        tr.ref_zone = ref_zone;
    end
    tr.east  = [tr.east; dummyE];
    tr.nort  = [tr.nort; dummyN];
    tr.zone  = [tr.zone; dummyZ];
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
        patch(vg(cg{i},1), vg(cg{i},2), log10(PolyMass(i)/1000));
    end
end
alpha .7        % Transparency
c = colorbar;   % Colorbar
ylabel(c, 'Log10 Mass (kg in the cell)');
% Plot points
plot(tr.lon,tr.lat, '.r')
xlabel('Longitude');
ylabel('Latitude');
plot_google_map


res_voron(tr.gClass, vorWt);

% Calculations and results                        
function res_voron(mClass, vorWt) 
global tr tmp t r

% Convert mm 2 phi
if get(t.diam_mm, 'Value') == 1
    pClass = -log2(tr.gClass');
else
    pClass = tr.gClass';
end

% Cumulative
cum = zeros(size(pClass));
for i = 1:length(cum)
    cum(i) = 1-sum(vorWt(1:i))/sum(vorWt);
end

tmp = [mClass', pClass, vorWt, cum];
data_table = cell(size(tmp,1),4);
for i = 1:size(tmp,1)
    for j = 1:4
        data_table{i,j} = sprintf('%.2f', tmp(i,j));
    end
end

gs_50       = get_perc(.5, cum, pClass);    % Median GS
gs_16       = get_perc(.16, cum, pClass);   % 14th percentile
gs_84       = get_perc(.84, cum, pClass);   % 86th percentile
gs_std_new  = (gs_84 - gs_16)/2;            % std

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
            'position', [.475 .175 .323 .112],...
            'BackgroundColor', [.3 .3 .3; .2 .2 .2],...
            'ColumnName', {'phi', 'mm'},...
            'ColumnWidth', {60, 60},...
            'ColumnFormat', {'numeric', 'numeric'},...
            'ColumnEditable', [false false],...
            'Data', [gs_50, 2^-(gs_50); gs_std_new, 2^-(gs_std_new)],...
            'RowStriping', 'on',...
            'RowName', {'Median', 'Sigma'},...
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

    
% Additional functions used throughout the code
function [east_c, zone_c] = correct_utm(east, zone, ref_zone, lat)
dummy  = unique(cellstr(zone));
dummy1 = sscanf(dummy{1}, '%d');
dummy2 = sscanf(dummy{2}, '%d');
dummy3 = dummy{1}(length(dummy{1}));
min_z  = min([dummy1, dummy2]);
max_z  = max([dummy1, dummy2]);
lim_W  = utmzone([num2str(min_z), dummy3]);
lim_E  = utmzone([num2str(max_z), dummy3]);

east_c = east;
zone_c = cellstr(zone);

if strcmp([num2str(min_z), ' ', dummy3], ref_zone)
    idx = find(~strcmp(zone, ref_zone));
    for i = 1:length(idx)
        east_c(idx(i)) = east(idx(i)) + deg2utm(lat(idx(i)), lim_W(4)-0.001) - deg2utm(lat(idx(i)), lim_E(3)+0.001);
        zone_c{idx(i)} = ref_zone;
    end
elseif strcmp([num2str(max_z), ' ', dummy3], ref_zone)
    idx = find(~strcmp(zone, ref_zone));
    for i = 1:length(idx)
        east_c(idx(i)) = east(idx(i)) - deg2utm(lat(idx(i)), lim_W(4)-0.001) + deg2utm(lat(idx(i)), lim_E(3)+0.001);
        zone_c{idx(i)} = ref_zone;
    end
end

zone_c = char(zone_c);

function out = get_perc(P, cum, pClass)
        % Median GS
[~, idx] = min(abs(cum-P));
if cum(idx) < P
    cum_low  = cum(idx);
    phi_low  = pClass(idx);
    cum_high = cum(idx+1);
    phi_high = pClass(idx+1);
else
    cum_low  = cum(idx-1);
    phi_low  = pClass(idx-1);
    cum_high = cum(idx);
    phi_high = pClass(idx);
end

out = ((phi_high - phi_low) / (cum_high - cum_low) * (P - cum_low)) + phi_low;

%% END OF TOTGS CODE
function varargout = plot_google_map(varargin)
% function h = plot_google_map(varargin)
% Plots a google map on the current axes using the Google Static Maps API
%
% USAGE:
% h = plot_google_map(Property, Value,...)
% Plots the map on the given axes. Used also if no output is specified
%
% Or:
% [lonVec latVec imag] = plot_google_map(Property, Value,...)
% Returns the map without plotting it
%
% PROPERTIES:
%    Axis           - Axis handle. If not given, gca is used. (LP)
%    Height (640)   - Height of the image in pixels (max 640)
%    Width  (640)   - Width of the image in pixels (max 640)
%    Scale (2)      - (1/2) Resolution scale factor. Using Scale=2 will
%                     double the resulotion of the downloaded image (up
%                     to 1280x1280) and will result in finer rendering,
%                     but processing time will be longer.
%    MapType        - ('roadmap') Type of map to return. Any of [roadmap, 
%                     satellite, terrain, hybrid]. See the Google Maps API for
%                     more information. 
%    Alpha (1)      - (0-1) Transparency level of the map (0 is fully
%                     transparent). While the map is always moved to the
%                     bottom of the plot (i.e. will not hide previously
%                     drawn items), this can be useful in order to increase
%                     readability if many colors are plotted 
%                     (using SCATTER for example).
%    ShowLabels (1) - (0/1) Controls whether to display city/street textual labels on the map
%    Language       - (string) A 2 letter ISO 639-1 language code for displaying labels in a 
%                     local language instead of English (where available).
%                     For example, for Chinese use:
%                     plot_google_map('language','zh')
%                     For the list of codes, see:
%                     http://en.wikipedia.org/wiki/List_of_ISO_639-1_codes
%    Marker         - The marker argument is a text string with fields
%                     conforming to the Google Maps API. The
%                     following are valid examples:
%                     '43.0738740,-70.713993' (default midsize orange marker)
%                     '43.0738740,-70.713993,blue' (midsize blue marker)
%                     '43.0738740,-70.713993,yellowa' (midsize yellow
%                     marker with label "A")
%                     '43.0738740,-70.713993,tinyredb' (tiny red marker
%                     with label "B")
%    Refresh (1)    - (0/1) defines whether to automatically refresh the
%                     map upon zoom/pan action on the figure.
%    AutoAxis (1)   - (0/1) defines whether to automatically adjust the axis
%                     of the plot to avoid the map being stretched.
%                     This will adjust the span to be correct
%                     according to the shape of the map axes.
%    FigureResizeUpdate (1) - (0/1) defines whether to automatically refresh the
%                     map upon resizing the figure. This will ensure map
%                     isn't stretched after figure resize.
%    APIKey         - (string) set your own API key which you obtained from Google: 
%                     http://developers.google.com/maps/documentation/staticmaps/#api_key
%                     This will enable up to 25,000 map requests per day, 
%                     compared to a few hundred requests without a key. 
%                     To set the key, use:
%                     plot_google_map('APIKey','SomeLongStringObtaindFromGoogle')
%                     You need to do this only once to set the key.
%                     To disable the use of a key, use:
%                     plot_google_map('APIKey','')
%
% OUTPUT:
%    h              - Handle to the plotted map
%
%    lonVect        - Vector of Longidute coordinates (WGS84) of the image 
%    latVect        - Vector of Latidute coordinates (WGS84) of the image 
%    imag           - Image matrix (height,width,3) of the map
%
% EXAMPLE - plot a map showing some capitals in Europe:
%    lat = [48.8708   51.5188   41.9260   40.4312   52.523   37.982];
%    lon = [2.4131    -0.1300    12.4951   -3.6788    13.415   23.715];
%    plot(lon,lat,'.r','MarkerSize',20)
%    plot_google_map
%
% References:
%  http://www.mathworks.com/matlabcentral/fileexchange/24113
%  http://www.maptiler.org/google-maps-coordinates-tile-bounds-projection/
%  http://developers.google.com/maps/documentation/staticmaps/
%
% Acknowledgements:
%  Val Schmidt for his submission of get_google_map.m
%
% Author:
%  Zohar Bar-Yehuda
%
% Version 1.6 - 12/11/2015
%       - Use system temp folder for writing image files (with fallback to current dir if missing write permissions)
% Version 1.5 - 20/11/2014
%       - Support for MATLAB R2014b
%       - several fixes complex layouts: several maps in one figure, 
%         map inside a panel, specifying axis handle as input (thanks to Luke Plausin)
% Version 1.4 - 25/03/2014
%       - Added the language parameter for showing labels in a local language
%       - Display the URL on error to allow easier debugging of API errors
% Version 1.3 - 06/10/2013
%       - Improved functionality of AutoAxis, which now handles any shape of map axes. 
%         Now also updates the extent of the map if the figure is resized.
%       - Added the showLabels parameter which allows hiding the textual labels on the map.
% Version 1.2 - 16/06/2012
%       - Support use of the "scale=2" parameter by default for finer rendering (set scale=1 if too slow).
%       - Auto-adjust axis extent so the map isn't stretched.
%       - Set and use an API key which enables a much higher usage volume per day.
% Version 1.1 - 25/08/2011

persistent apiKey useTemp
if isnumeric(apiKey)
    % first run, check if API key file exists
    if exist('api_key.mat','file')
        load api_key
    else
        apiKey = '';
    end
end

if isempty(useTemp)
    % first run, check we we have wrtie access to the temp folder
    try 
        tempfilename = tempname;
        fid = fopen(tempfilename, 'w');
        if fid > 0
            fclose(fid);
            useTemp = true;
            delete(tempfilename);
        else
            % Don't have write access to temp folder or it doesn't exist, fallback to current dir
            useTemp = false;
        end
    catch
        % in case tempname fails for some reason
        useTemp = false;
    end
end

hold on

% Default parametrs
axHandle = gca;
height = 640;
width = 640;
scale = 2;
maptype = 'roadmap';
alphaData = 1;
autoRefresh = 1;
figureResizeUpdate = 1;
autoAxis = 1;
showLabels = 1;
language = '';
markeridx = 1;
markerlist = {};

% Handle input arguments
if nargin >= 2
    for idx = 1:2:length(varargin)
        switch lower(varargin{idx})
            case 'axis'
                axHandle = varargin{idx+1};
            case 'height'
                height = varargin{idx+1};
            case 'width'
                width = varargin{idx+1};
            case 'maptype'
                maptype = varargin{idx+1};
            case 'alpha'
                alphaData = varargin{idx+1};
            case 'refresh'
                autoRefresh = varargin{idx+1};
            case 'showlabels'
                showLabels = varargin{idx+1};
            case 'figureresizeupdate'
                figureResizeUpdate = varargin{idx+1};
            case 'language'
                language = varargin{idx+1};
            case 'marker'
                markerlist{markeridx} = varargin{idx+1};
                markeridx = markeridx + 1;
            case 'autoaxis'
                autoAxis = varargin{idx+1};
            case 'apikey'
                apiKey = varargin{idx+1}; % set new key
                % save key to file
                funcFile = which('plot_google_map.m');
                pth = fileparts(funcFile);
                keyFile = fullfile(pth,'api_key.mat');
                save(keyFile,'apiKey')
            otherwise
                error(['Unrecognized variable: ' varargin{idx}])
        end
    end
end
if height > 640
    height = 640;
end
if width > 640
    width = 640;
end

% Store paramters in axis handle (for auto refresh callbacks)
ud = get(axHandle, 'UserData');
if isempty(ud)
    % explicitly set as struct to avoid warnings
    ud = struct;
end
ud.gmap_params = varargin;
set(axHandle, 'UserData', ud);

curAxis = axis(axHandle);
% Enforce Latitude constraints of EPSG:900913 
if curAxis(3) < -85
    curAxis(3) = -85;
end
if curAxis(4) > 85
    curAxis(4) = 85;
end
% Enforce longitude constrains
if curAxis(1) < -180
    curAxis(1) = -180;
end
if curAxis(1) > 180
    curAxis(1) = 0;
end
if curAxis(2) > 180
    curAxis(2) = 180;
end
if curAxis(2) < -180
    curAxis(2) = 0;
end

if isequal(curAxis,[0 1 0 1]) % probably an empty figure
    % display world map
    curAxis = [-200 200 -85 85];
    axis(curAxis)
end


if autoAxis
    % adjust current axis limit to avoid strectched maps
    [xExtent,yExtent] = latLonToMeters(curAxis(3:4), curAxis(1:2) );
    xExtent = diff(xExtent); % just the size of the span
    yExtent = diff(yExtent); 
    % get axes aspect ratio
    drawnow
    org_units = get(axHandle,'Units');
    set(axHandle,'Units','Pixels')
    ax_position = get(axHandle,'position');        
    set(axHandle,'Units',org_units)
    aspect_ratio = ax_position(4) / ax_position(3);
    
    if xExtent*aspect_ratio > yExtent        
        centerX = mean(curAxis(1:2));
        centerY = mean(curAxis(3:4));
        spanX = (curAxis(2)-curAxis(1))/2;
        spanY = (curAxis(4)-curAxis(3))/2;
       
        % enlarge the Y extent
        spanY = spanY*xExtent*aspect_ratio/yExtent; % new span
        if spanY > 85
            spanX = spanX * 85 / spanY;
            spanY = spanY * 85 / spanY;
        end
        curAxis(1) = centerX-spanX;
        curAxis(2) = centerX+spanX;
        curAxis(3) = centerY-spanY;
        curAxis(4) = centerY+spanY;
    elseif yExtent > xExtent*aspect_ratio
        
        centerX = mean(curAxis(1:2));
        centerY = mean(curAxis(3:4));
        spanX = (curAxis(2)-curAxis(1))/2;
        spanY = (curAxis(4)-curAxis(3))/2;
        % enlarge the X extent
        spanX = spanX*yExtent/(xExtent*aspect_ratio); % new span
        if spanX > 180
            spanY = spanY * 180 / spanX;
            spanX = spanX * 180 / spanX;
        end
        
        curAxis(1) = centerX-spanX;
        curAxis(2) = centerX+spanX;
        curAxis(3) = centerY-spanY;
        curAxis(4) = centerY+spanY;
    end            
    % Enforce Latitude constraints of EPSG:900913
    if curAxis(3) < -85
        curAxis(3:4) = curAxis(3:4) + (-85 - curAxis(3));
    end
    if curAxis(4) > 85
        curAxis(3:4) = curAxis(3:4) + (85 - curAxis(4));
    end
    axis(axHandle, curAxis); % update axis as quickly as possible, before downloading new image
    drawnow
end

% Delete previous map from plot (if exists)
if nargout <= 1 % only if in plotting mode
    curChildren = get(axHandle,'children');
    map_objs = findobj(curChildren,'tag','gmap');
    bd_callback = [];
    for idx = 1:length(map_objs)
        if ~isempty(get(map_objs(idx),'ButtonDownFcn'))
            % copy callback properties from current map
            bd_callback = get(map_objs(idx),'ButtonDownFcn');
        end
    end
    delete(map_objs)
    
end

% Calculate zoom level for current axis limits
[xExtent,yExtent] = latLonToMeters(curAxis(3:4), curAxis(1:2) );
minResX = diff(xExtent) / width;
minResY = diff(yExtent) / height;
minRes = max([minResX minResY]);
tileSize = 256;
initialResolution = 2 * pi * 6378137 / tileSize; % 156543.03392804062 for tileSize 256 pixels
zoomlevel = floor(log2(initialResolution/minRes));

% Enforce valid zoom levels
if zoomlevel < 0 
    zoomlevel = 0;
end
if zoomlevel > 19 
    zoomlevel = 19;
end

% Calculate center coordinate in WGS1984
lat = (curAxis(3)+curAxis(4))/2;
lon = (curAxis(1)+curAxis(2))/2;

% Construct query URL
preamble = 'http://maps.googleapis.com/maps/api/staticmap';
location = ['?center=' num2str(lat,10) ',' num2str(lon,10)];
zoomStr = ['&zoom=' num2str(zoomlevel)];
sizeStr = ['&scale=' num2str(scale) '&size=' num2str(width) 'x' num2str(height)];
maptypeStr = ['&maptype=' maptype ];
if ~isempty(apiKey)
    keyStr = ['&key=' apiKey];
else
    keyStr = '';
end
markers = '&markers=';
for idx = 1:length(markerlist)
    if idx < length(markerlist)
        markers = [markers markerlist{idx} '%7C'];
    else
        markers = [markers markerlist{idx}];
    end
end
if showLabels == 0
    labelsStr = '&style=feature:all|element:labels|visibility:off';
else
    labelsStr = '';
end
if ~isempty(language)
    languageStr = ['&language=' language];
else
    languageStr = '';
end
    
if ismember(maptype,{'satellite','hybrid'})
    filename = 'tmp.jpg';
    format = '&format=jpg';
    convertNeeded = 0;
else
    filename = 'tmp.png';
    format = '&format=png';
    convertNeeded = 1;
end
sensor = '&sensor=false';
url = [preamble location zoomStr sizeStr maptypeStr format markers labelsStr languageStr sensor keyStr];

% Get the image
if useTemp
    filepath = fullfile(tempdir, filename);
else
    filepath = filename;
end

try
    urlwrite(url,filepath);
catch % error downloading map
    warning(['Unable to download map form Google Servers.\n' ...
        'Matlab error was: %s\n\n' ...
        'Possible reasons: missing write permissions, no network connection, quota exceeded, or some other error.\n' ...
        'Consider using an API key if quota problems persist.\n\n' ...
        'To debug, try pasting the following URL in your browser, which may result in a more informative error:\n%s'], lasterr, url);
    varargout{1} = [];
    varargout{2} = [];
    varargout{3} = [];
    return
end
[M Mcolor] = imread(filepath);
M = cast(M,'double');
delete(filepath); % delete temp file
width = size(M,2);
height = size(M,1);

% Calculate a meshgrid of pixel coordinates in EPSG:900913
centerPixelY = round(height/2);
centerPixelX = round(width/2);
[centerX,centerY] = latLonToMeters(lat, lon ); % center coordinates in EPSG:900913
curResolution = initialResolution / 2^zoomlevel/scale; % meters/pixel (EPSG:900913)
xVec = centerX + ((1:width)-centerPixelX) * curResolution; % x vector
yVec = centerY + ((height:-1:1)-centerPixelY) * curResolution; % y vector
[xMesh,yMesh] = meshgrid(xVec,yVec); % construct meshgrid 

% convert meshgrid to WGS1984
[lonMesh,latMesh] = metersToLatLon(xMesh,yMesh);

% We now want to convert the image from a colormap image with an uneven
% mesh grid, into an RGB truecolor image with a uniform grid.
% This would enable displaying it with IMAGE, instead of PCOLOR.
% Advantages are:
% 1) faster rendering
% 2) makes it possible to display together with other colormap annotations (PCOLOR, SCATTER etc.)

% Convert image from colormap type to RGB truecolor (if PNG is used)
if convertNeeded
    imag = zeros(height,width,3);
    for idx = 1:3
        imag(:,:,idx) = reshape(Mcolor(M(:)+1+(idx-1)*size(Mcolor,1)),height,width);
    end
else
    imag = M/255;
end

% Next, project the data into a uniform WGS1984 grid
sizeFactor = 1; % factoring of new image
uniHeight = round(height*sizeFactor);
uniWidth = round(width*sizeFactor);
latVect = linspace(latMesh(1,1),latMesh(end,1),uniHeight);
lonVect = linspace(lonMesh(1,1),lonMesh(1,end),uniWidth);
[uniLonMesh,uniLatMesh] = meshgrid(lonVect,latVect);
uniImag = zeros(uniHeight,uniWidth,3);

% old version (projection using INTERP2)
% for idx = 1:3
%      % 'nearest' method is the fastest. difference from other methods is neglible
%          uniImag(:,:,idx) =  interp2(lonMesh,latMesh,imag(:,:,idx),uniLonMesh,uniLatMesh,'nearest');
% end
uniImag =  myTurboInterp2(lonMesh,latMesh,imag,uniLonMesh,uniLatMesh);

if nargout <= 1 % plot map
    % display image
    hold(axHandle, 'on');
    h = image(lonVect,latVect,uniImag, 'Parent', axHandle);
    set(axHandle,'YDir','Normal')
    set(h,'tag','gmap')
    set(h,'AlphaData',alphaData)
    
    % add a dummy image to allow pan/zoom out to x2 of the image extent
    h_tmp = image(lonVect([1 end]),latVect([1 end]),zeros(2),'Visible','off', 'Parent', axHandle);
    set(h_tmp,'tag','gmap')
    
    % older version (display without conversion to uniform grid)
    % h =pcolor(lonMesh,latMesh,(M));
    % colormap(Mcolor)
    % caxis([0 255])
    % warning off % to avoid strange rendering warnings
    % shading flat
   
    uistack(h,'bottom') % move map to bottom (so it doesn't hide previously drawn annotations)
    axis(axHandle, curAxis) % restore original zoom
    if nargout == 1
        varargout{1} = h;
    end
    
    % if auto-refresh mode - override zoom callback to allow autumatic 
    % refresh of map upon zoom actions.
    figHandle = axHandle;
    while ~strcmpi(get(figHandle, 'Type'), 'figure')
        % Recursively search for parent figure in case axes are in a panel
        figHandle = get(figHandle, 'Parent');
    end
    
    zoomHandle = zoom(axHandle);   
    panHandle = pan(figHandle); % This isn't ideal, doesn't work for contained axis    
    if autoRefresh        
        set(zoomHandle,'ActionPostCallback',@update_google_map);          
        set(panHandle, 'ActionPostCallback', @update_google_map);        
    else % disable zoom override
        set(zoomHandle,'ActionPostCallback',[]);
        set(panHandle, 'ActionPostCallback',[]);
    end
    
    % set callback for figure resize function, to update extents if figure
    % is streched.
    if figureResizeUpdate &&isempty(get(figHandle, 'ResizeFcn'))
        % set only if not already set by someone else
        set(figHandle, 'ResizeFcn', @update_google_map_fig);       
    end    
    
    % set callback properties 
    set(h,'ButtonDownFcn',bd_callback);
else % don't plot, only return map
    varargout{1} = lonVect;
    varargout{2} = latVect;
    varargout{3} = uniImag;
end


% Coordinate transformation functions

function [lon,lat] = metersToLatLon(x,y)
% Converts XY point from Spherical Mercator EPSG:900913 to lat/lon in WGS84 Datum
originShift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
lon = (x ./ originShift) * 180;
lat = (y ./ originShift) * 180;
lat = 180 / pi * (2 * atan( exp( lat * pi / 180)) - pi / 2);

function [x,y] = latLonToMeters(lat, lon )
% Converts given lat/lon in WGS84 Datum to XY in Spherical Mercator EPSG:900913"
originShift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
x = lon * originShift / 180;
y = log(tan((90 + lat) * pi / 360 )) / (pi / 180);
y = y * originShift / 180;

function ZI = myTurboInterp2(X,Y,Z,XI,YI)
% An extremely fast nearest neighbour 2D interpolation, assuming both input
% and output grids consist only of squares, meaning:
% - uniform X for each column
% - uniform Y for each row
XI = XI(1,:);
X = X(1,:);
YI = YI(:,1);
Y = Y(:,1);

xiPos = nan*ones(size(XI));
xLen = length(X);
yiPos = nan*ones(size(YI));
yLen = length(Y);
% find x conversion
xPos = 1;
for idx = 1:length(xiPos)
    if XI(idx) >= X(1) && XI(idx) <= X(end)
        while xPos < xLen && X(xPos+1)<XI(idx)
            xPos = xPos + 1;
        end
        diffs = abs(X(xPos:xPos+1)-XI(idx));
        if diffs(1) < diffs(2)
            xiPos(idx) = xPos;
        else
            xiPos(idx) = xPos + 1;
        end
    end
end
% find y conversion
yPos = 1;
for idx = 1:length(yiPos)
    if YI(idx) <= Y(1) && YI(idx) >= Y(end)
        while yPos < yLen && Y(yPos+1)>YI(idx)
            yPos = yPos + 1;
        end
        diffs = abs(Y(yPos:yPos+1)-YI(idx));
        if diffs(1) < diffs(2)
            yiPos(idx) = yPos;
        else
            yiPos(idx) = yPos + 1;
        end
    end
end
ZI = Z(yiPos,xiPos,:);

function update_google_map(obj,evd)
% callback function for auto-refresh
drawnow;
try
    axHandle = evd.Axes;
catch ex
    % Event doesn't contain the correct axes. Panic!
    axHandle = gca;
end
ud = get(axHandle, 'UserData');
if isfield(ud, 'gmap_params')
    params = ud.gmap_params;
    plot_google_map(params{:});
end

function update_google_map_fig(obj,evd)
% callback function for auto-refresh
drawnow;
axes_objs = findobj(get(gcf,'children'),'type','axes');
for idx = 1:length(axes_objs)
    if ~isempty(findobj(get(axes_objs(idx),'children'),'tag','gmap'));
        ud = get(axes_objs(idx), 'UserData');
        if isfield(ud, 'gmap_params')
            params = ud.gmap_params;
        else
            params = {};
        end
        
        % Add axes to inputs if needed
        if ~sum(strcmpi(params, 'Axis'))
            params = [params, {'Axis', axes_objs(idx)}];
        end
        plot_google_map(params{:});
    end
end

function  [x,y,utmzone] = deg2utm(Lat,Lon)
% -------------------------------------------------------------------------
% [x,y,utmzone] = deg2utm(Lat,Lon)
%
% Author: 
%   Rafael Palacios
%   Universidad Pontificia Comillas
%   Madrid, Spain
% Version: Apr/06, Jun/06, Aug/06, Aug/06
% Aug/06: fixed a problem (found by Rodolphe Dewarrat) related to southern 
%    hemisphere coordinates. 
% Aug/06: corrected m-Lint warnings
%-------------------------------------------------------------------------

% Argument checking
%
error(nargchk(2, 2, nargin));  %2 arguments required
n1=length(Lat);
n2=length(Lon);
if (n1~=n2)
   error('Lat and Lon vectors should have the same length');
end


% Memory pre-allocation
%
x=zeros(n1,1);
y=zeros(n1,1);
utmzone(n1,:)='60 X';

% Main Loop
%
for i=1:n1
   la=Lat(i);
   lo=Lon(i);

   sa = 6378137.000000 ; sb = 6356752.314245;
         
   %e = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sa;
   e2 = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sb;
   e2cuadrada = e2 ^ 2;
   c = ( sa ^ 2 ) / sb;
   %alpha = ( sa - sb ) / sa;             %f
   %ablandamiento = 1 / alpha;   % 1/f

   lat = la * ( pi / 180 );
   lon = lo * ( pi / 180 );

   Huso = fix( ( lo / 6 ) + 31);
   S = ( ( Huso * 6 ) - 183 );
   deltaS = lon - ( S * ( pi / 180 ) );

   if (la<-72), Letra='C';
   elseif (la<-64), Letra='D';
   elseif (la<-56), Letra='E';
   elseif (la<-48), Letra='F';
   elseif (la<-40), Letra='G';
   elseif (la<-32), Letra='H';
   elseif (la<-24), Letra='J';
   elseif (la<-16), Letra='K';
   elseif (la<-8), Letra='L';
   elseif (la<0), Letra='M';
   elseif (la<8), Letra='N';
   elseif (la<16), Letra='P';
   elseif (la<24), Letra='Q';
   elseif (la<32), Letra='R';
   elseif (la<40), Letra='S';
   elseif (la<48), Letra='T';
   elseif (la<56), Letra='U';
   elseif (la<64), Letra='V';
   elseif (la<72), Letra='W';
   else Letra='X';
   end

   a = cos(lat) * sin(deltaS);
   epsilon = 0.5 * log( ( 1 +  a) / ( 1 - a ) );
   nu = atan( tan(lat) / cos(deltaS) ) - lat;
   v = ( c / ( ( 1 + ( e2cuadrada * ( cos(lat) ) ^ 2 ) ) ) ^ 0.5 ) * 0.9996;
   ta = ( e2cuadrada / 2 ) * epsilon ^ 2 * ( cos(lat) ) ^ 2;
   a1 = sin( 2 * lat );
   a2 = a1 * ( cos(lat) ) ^ 2;
   j2 = lat + ( a1 / 2 );
   j4 = ( ( 3 * j2 ) + a2 ) / 4;
   j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat) ) ^ 2) ) / 3;
   alfa = ( 3 / 4 ) * e2cuadrada;
   beta = ( 5 / 3 ) * alfa ^ 2;
   gama = ( 35 / 27 ) * alfa ^ 3;
   Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
   xx = epsilon * v * ( 1 + ( ta / 3 ) ) + 500000;
   yy = nu * v * ( 1 + ta ) + Bm;

   if (yy<0)
       yy=9999999+yy;
   end

   x(i)=xx;
   y(i)=yy;
   utmzone(i,:)=sprintf('%02d %c',Huso,Letra);
end

function  [Lat,Lon] = utm2deg(xx,yy,utmzone)
% -------------------------------------------------------------------------
% [Lat,Lon] = utm2deg(x,y,utmzone)
%
% Description: Function to convert vectors of UTM coordinates into Lat/Lon vectors (WGS84).
% Some code has been extracted from UTMIP.m function by Gabriel Ruiz Martinez.
%
% Inputs:
%    x, y , utmzone.
%
% Outputs:
%    Lat: Latitude vector.   Degrees.  +ddd.ddddd  WGS84
%    Lon: Longitude vector.  Degrees.  +ddd.ddddd  WGS84
%
% Example 1:
% x=[ 458731;  407653;  239027;  230253;  343898;  362850];
% y=[4462881; 5126290; 4163083; 3171843; 4302285; 2772478];
% utmzone=['30 T'; '32 T'; '11 S'; '28 R'; '15 S'; '51 R'];
% [Lat, Lon]=utm2deg(x,y,utmzone);
% fprintf('%11.6f ',lat)
%    40.315430   46.283902   37.577834   28.645647   38.855552   25.061780
% fprintf('%11.6f ',lon)
%    -3.485713    7.801235 -119.955246  -17.759537  -94.799019  121.640266
%
% Example 2: If you need Lat/Lon coordinates in Degrees, Minutes and Seconds
% [Lat, Lon]=utm2deg(x,y,utmzone);
% LatDMS=dms2mat(deg2dms(Lat))
%LatDMS =
%    40.00         18.00         55.55
%    46.00         17.00          2.01
%    37.00         34.00         40.17
%    28.00         38.00         44.33
%    38.00         51.00         19.96
%    25.00          3.00         42.41
% LonDMS=dms2mat(deg2dms(Lon))
%LonDMS =
%    -3.00         29.00          8.61
%     7.00         48.00          4.40
%  -119.00         57.00         18.93
%   -17.00         45.00         34.33
%   -94.00         47.00         56.47
%   121.00         38.00         24.96
%
% Author: 
%   Rafael Palacios
%   Universidad Pontificia Comillas
%   Madrid, Spain
% Version: Apr/06, Jun/06, Aug/06
% Aug/06: corrected m-Lint warnings
%-------------------------------------------------------------------------

% Argument checking
%
error(nargchk(3, 3, nargin)); %3 arguments required
n1=length(xx);
n2=length(yy);
n3=size(utmzone,1);
if (n1~=n2 || n1~=n3)
   error('x,y and utmzone vectors should have the same number or rows');
end
c=size(utmzone,2);
if (c~=4)
   error('utmzone should be a vector of strings like "30 T"');
end

   
 
% Memory pre-allocation
%
Lat=zeros(n1,1);
Lon=zeros(n1,1);


% Main Loop
%
for i=1:n1
   if (utmzone(i,4)>'X' || utmzone(i,4)<'C')
      fprintf('utm2deg: Warning utmzone should be a vector of strings like "30 T", not "30 t"\n');
   end
   if (utmzone(i,4)>'M')
      hemis='N';   % Northern hemisphere
   else
      hemis='S';
   end

   x=xx(i);
   y=yy(i);
   zone=str2double(utmzone(i,1:2));

   sa = 6378137.000000 ; sb = 6356752.314245;
  
%   e = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sa;
   e2 = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sb;
   e2cuadrada = e2 ^ 2;
   c = ( sa ^ 2 ) / sb;
%   alpha = ( sa - sb ) / sa;             %f
%   ablandamiento = 1 / alpha;   % 1/f

   X = x - 500000;
   
   if hemis == 'S' || hemis == 's'
       Y = y - 10000000;
   else
       Y = y;
   end
    
   S = ( ( zone * 6 ) - 183 ); 
   lat =  Y / ( 6366197.724 * 0.9996 );                                    
   v = ( c / ( ( 1 + ( e2cuadrada * ( cos(lat) ) ^ 2 ) ) ) ^ 0.5 ) * 0.9996;
   a = X / v;
   a1 = sin( 2 * lat );
   a2 = a1 * ( cos(lat) ) ^ 2;
   j2 = lat + ( a1 / 2 );
   j4 = ( ( 3 * j2 ) + a2 ) / 4;
   j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat) ) ^ 2) ) / 3;
   alfa = ( 3 / 4 ) * e2cuadrada;
   beta = ( 5 / 3 ) * alfa ^ 2;
   gama = ( 35 / 27 ) * alfa ^ 3;
   Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
   b = ( Y - Bm ) / v;
   Epsi = ( ( e2cuadrada * a^ 2 ) / 2 ) * ( cos(lat) )^ 2;
   Eps = a * ( 1 - ( Epsi / 3 ) );
   nab = ( b * ( 1 - Epsi ) ) + lat;
   senoheps = ( exp(Eps) - exp(-Eps) ) / 2;
   Delt = atan(senoheps / (cos(nab) ) );
   TaO = atan(cos(Delt) * tan(nab));
   longitude = (Delt *(180 / pi ) ) + S;
   latitude = ( lat + ( 1 + e2cuadrada* (cos(lat)^ 2) - ( 3 / 2 ) * e2cuadrada * sin(lat) * cos(lat) * ( TaO - lat ) ) * ( TaO - lat ) ) * ...
                    (180 / pi);
   
   Lat(i)=latitude;
   Lon(i)=longitude;
   
end