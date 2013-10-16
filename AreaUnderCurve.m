%   AREA UNDER THE CURVE
% minarikova.lenka@gmail.com
% v 1.0

function [] = AreaUnderCurve(directory, w_ppm, bdwtd, trnct, control) %, CSI_shft_ud, CSI_shft_lr)
%  directory = '~/Desktop/Current/phantoms/phantom_2013_10_10/';
%  w_ppm = 4.74;
%  bdwtd = 120;
%  trnct = 0;
%  control = 1;


% !!!!!!!!!!!!!!!!!! readme !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% for working you need my other function called read_ascconv_lenk.m 
%   for reading parameters from the dicom
% directory = '~/Patient_name' - where is a directory called "Spec" with 
%   a dicom 3D CSI file
% cho_ppm = 3.2 - exact position of Choline peak if the spectrum is set
%   to begin at 8.76 ppm and ends at 0.64, with 1000 Hz bandwidth
% bdwtd = 50 - bandwidth of the peak
% trnct = number of points to truncate at the end of FID
% control = 1 - at first you need to control if the baseline correction 
%   is working properly and everything is set all right
% shft_ud: for shifted CSI: up -0.% down +0.%
% shft_lr: for shifted CSI: left -0.% right +0.%

% the output is maximal, mean value of all SNRs of Cho and a table 
%   with all SNRs in one row, all saved in txt files in Spec directory

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% search variables in spectroscopy file:
tic;

directory_spec = strcat(directory,'Spec/');
cd(directory_spec);
disp(strcat('Processing:',directory));
%% <<<<<<<<<<<<<<<<< Parameters >>>>>>>>>>>>>>>>>>
snr.cho = w_ppm;
press_big = 0; % choose if you want more voxel (even that no totally inside the press box)
%snr.noise = 7.4; % Aus diesem Bereich wird das Rausch Siganl genommen +0.5 bis -0.5
% Anfangs und Endewert fuer die ppm Skala eintragen
anfang = 8.76; % Anfang der ppm Skala 8.76 normally
ende = 0.64; % Ende der ppm Skala 0.64 
% truncate, replace the last few points in fid with 0s:

% read the .txt header from spectroscopic file made of the text part
% beginnig by "### ASCCONV BEGIN ###" and ending by "### ASCCONV END ###"
slices = dir(pwd);
for k = length(slices):-1:1 % find .IMA file
    fname = slices(k).name;
    if fname(1) == '.'
        slices(k) = [];
    elseif strcat(fname(end-2),fname(end-1),fname(end)) ~= strcat('IMA')
        slices(k) = [];
    end
end

disp('Searching for spec file...');
spect = dir(pwd);
for k = length(spect):-1:1
    fname = spect(k).name;
    if fname(1) == '.' %|| strcat(spect(k).name) == strcat('Header.txt') || spect(k).name == strcat('Signal.txt')
        spect(k) = [];
    end
    % check filetype
    strg = (spect(k).name); % name of the file
    points = strfind(strg,'.');
    last_point = max(points);
    filetype_end = (strg(last_point:end));
    end_IMA=strcmp(filetype_end,'.IMA');    
    if end_IMA == 0;
        spect(k) = [];
        continue
    end    
end
voxel = read_ascconv_lenk(slices(1).name); % read parrameter from spectroscopy dicom

%% determine CSI-in-press parameters:
vecSize = voxel.vecSize;
voxel.size_z = voxel.FoV_z / voxel.number_z;
voxel.size_y = voxel.FoV_y / voxel.number_y;
voxel.size_x = voxel.FoV_x / voxel.number_x;
% number of steps in each direction:
if press_big == 1
    voxel.step_x = floor((voxel.number_x  - ceil((voxel.FoV_x - voxel.p_fov_x) / voxel.size_x)) / 2); %include also voxels touching the pressbox fringes
    voxel.step_y = floor((voxel.number_y  - ceil((voxel.FoV_y - voxel.p_fov_y) / voxel.size_y)) / 2);
    voxel.step_z = floor((voxel.number_z  - ceil((voxel.FoV_z - voxel.p_fov_z) / voxel.size_z)) / 2);
elseif press_big == 0
    voxel.step_x = voxel.number_x / 2 - fix((voxel.FoV_x - voxel.p_fov_x) / 2 / voxel.size_x) - 1; %include only voxels inside the PRESS box
    voxel.step_y = voxel.number_y / 2 - fix((voxel.FoV_y - voxel.p_fov_y) / 2 / voxel.size_y) - 1;
    voxel.step_z = voxel.number_z / 2 - fix((voxel.FoV_z - voxel.p_fov_z) / 2 / voxel.size_z) - 1;
else
    voxel.step_x = voxel.number_x / 2 - fix((voxel.FoV_x - voxel.p_fov_x) / 2 / voxel.size_x); %include voxels with PRESSbox fringes inside
    voxel.step_y = voxel.number_y / 2 - fix((voxel.FoV_y - voxel.p_fov_y) / 2 / voxel.size_y);
    voxel.step_z = voxel.number_z / 2 - fix((voxel.FoV_z - voxel.p_fov_z) / 2 / voxel.size_z);
end


%% Slices v0.3
% by Matthias Riha, modified by minarikova.lenka@gmail.com

% Nur fuor Traversal Ebenen zu gebrauchen!!!!
% Das File bf.m muess im selben Order sein!!!!!!!

%% DEFINITIONS

ROW = voxel.number_x; % x
COL = voxel.number_y; % y
SLC = voxel.number_z; % z

bbox.col_start = COL / 2 - voxel.step_y + 1; %Sagital von
bbox.col_end = COL / 2 + voxel.step_y; %bis
bbox.row_start = ROW / 2 - voxel.step_x + 1; %Coronal von
bbox.row_end = ROW / 2 + voxel.step_x; %bis
bbox.slc_start = SLC / 2 - voxel.step_z + 1; %Diese beiden Werte muessen immer gleich sein
bbox.slc_end = SLC / 2 + voxel.step_z;

%% Einlesen der Daten
csi.file_in = strcat(spect(1,1).name); %Pfad immer an die Datei anpassen

% read CSI data
fid = fopen(csi.file_in,'r');
fseek(fid,-((vecSize*ROW*COL*SLC*2*4)),'eof');
csi.data = fread(fid,'float32');
fclose(fid);

csi.real = csi.data(1:2:end);
csi.imag = csi.data(2:2:end);
csi.complex = complex(csi.real,csi.imag);
%% create 4D matrix of the whole csi grid
csi.mat_complex = zeros(ROW,COL,SLC,vecSize);
csi.mat_real = zeros(ROW,COL,SLC,vecSize);

o = -1;
for z = 1:SLC
    for y = 1:COL
        for x = 1:ROW
            o = o + 1;
            csi.mat_cmplx(x,y,z,:) = csi.complex(o * vecSize + 1:(o + 1) * vecSize);	
            csi.mat_rl(x,y,z,:) = csi.real(o * vecSize + 1:(o + 1) * vecSize);
        end
    end
end
%%
i = 0;
SNR.mid = 0;
SNR.tab = 0;
vecSize = 2 * vecSize;
anfang = anfang * (-1);
c = 0;

% Define the ppm scale
div = (anfang) + ende;
int = div / vecSize;
e = anfang;
F = 0;
for i = 1:vecSize
    F(1,i) = e;
    e = e - int;
end
F = F * (-1);
Vecsize = (1:vecSize);
Vecsize = Vecsize.';
F_col = (F.'); % frequency
tabulk = [Vecsize,F_col];

for xx = 1:vecSize % define the values for baseline correction
    if round(100 * real(tabulk(xx,2))) == round(snr.cho * 100)
        r_sd = real(tabulk(xx - round(bdwtd / 2),1)); % right minimum next to choline peak
        l_sd = real(tabulk(xx + round(bdwtd / 2),1)); % left minimum next to Cho peak
        
        break
    end
end
SNR.main = 0;
lala = 0;
%% PRESSbox Daten FFT und in Txt File schreiben

for z = bbox.slc_start:bbox.slc_end
    b = 0;
    c = c + 1;
    for y = bbox.col_start:bbox.col_end
        a = 0;
        b = b + 1;
        for x = bbox.row_start:bbox.row_end
            lala = lala + 1;
            i = i + 1;
            a = a + 1;
            % 1 Voxel aus der 4D Matrix auslesen
            csi.rshpd_cmplx = reshape(csi.mat_cmplx(x,y,z,:),[],1);
            % FFT eines Voxels
            x1 = csi.rshpd_cmplx;
            x2 = csi.rshpd_cmplx(1:vecSize / 2 - trnct,1);
            rest1 = zeros(trnct,1);
            x1 = [x2;rest1];
            % ff = apod('Lorentz',vecSize / 2,1,250);
            % ff = ff';
            % x3 = x1.*ff;
            N = vecSize;
            % X = fft(x3,N);
            X = fft(x1,N);
            X = abs(fftshift(X)); % take the absolute value
            %Xnixfit = X;
            diff_ = var([X(r_sd - 28:r_sd);X(l_sd:l_sd + 28)]); % variance between two means of points at the edge of peak

            % Baseline FIT !!!!!!
            if control == 1
                disp(SNR.main);
                X = bf(X,[r_sd,l_sd],5,'linear','confirm');
                plot(F,X);
                set(gca,'XDir','reverse');
            else
                X = bf(X,[r_sd,l_sd],5,'linear');
            end
             
            % ppm scale + values from X
            test_snr = [Vecsize,F_col,X];
            yy = 0;
            zz = 0;
            w = 0;
            snr.min = 0;
            snr.max = 0;
            snr.base = 0;
            for xx = 1:vecSize
                snr.schleife = test_snr(xx,:,:);
                ppm_wert = snr.schleife(:,2);
                
                if ppm_wert < snr.cho + (bdwtd / 128) %approx for 3T
                    if ppm_wert > snr.cho - (bdwtd / 128)
                        %write all the values
                        zz = zz + 1;
                        snr.peak_amp(zz) = real(snr.schleife(:,3));
                    end
                end
            end
            %sum the whole peak
            SNR.main = sum(snr.peak_amp);
            SNR.tab(b,a,c) = SNR.main;
            
            % write it all to:
            SNR.reshaped(lala) = SNR.main;
        end
    end
end

%% treshold for all values:
%reshape and add coordinates as are in siemens:
c.pix_width = voxel.step_y * 2; % number of voxels in the pressbox - x axis
c.pix_height = voxel.step_x * 2; % -||- - y axis
c.pix_depth = voxel.step_z * 2; % -||- - z axis
ooo = 0;
for k = 1:c.pix_depth
    for j = 1:c.pix_width
        for i = 1:c.pix_height
            ooo = ooo + 1;
            SNR.w_coor(ooo,1) = i + (voxel.number_x / 2 - voxel.step_x);
            SNR.w_coor(ooo,2) = j + (voxel.number_y / 2 - voxel.step_y);
            SNR.w_coor(ooo,3) = k + (voxel.number_z / 2 - voxel.step_z);
            SNR.w_coor(ooo,4) = SNR.reshaped(1,ooo);
        end
    end
end

%% saving important things:
path = sprintf('Output_choSNR.txt');
fid_write = fopen(path,'w');
fprintf(fid_write,'%d\n',SNR.reshaped);
fclose(fid_write);

    pvc_mean = mean(SNR.reshaped);
    path = sprintf('Output_cho_mean.txt');
    fid_write = fopen(path,'w');
    fprintf(fid_write,'%d\n',pvc_mean);
    fclose(fid_write);
    disp(pvc_mean);

    pvc_max = max(SNR.reshaped);
    path = sprintf('Output_cho_max.txt');
    fid_write = fopen(path,'w');
    fprintf(fid_write,'%d\n',pvc_max);
    fclose(fid_write);
    disp(pvc_max);

dlmwrite('Output_choSNR_w_coor.txt', SNR.w_coor, 'delimiter', '\t', ...
         'precision', 6);
 
toc