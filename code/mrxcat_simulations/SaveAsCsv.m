%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================
function [data, fname] = SaveAsCsv(fname, plot_result )


cd('../../data/mrxcat_simulations');
foldernames = glob ("snr*");

for i = 1:numel(foldernames);
%files.name
%files.folder
  [~, foldername] = fileparts (foldernames{i});
  cd(foldername);
  files = glob('*.cpx');
  global file;
  [~, file] = fileparts (files{1});
  count = size(files)(1);
  %mkdir(foldername + "/csvs");
  global dest;
  dest = strcat(pwd, "/csvs");

    if nargin < 2 #|| ~exist(fname,'file')
        % get MRXCAT output (.cpx file)
        #[filename,pathname] = uigetfile({'*.cpx'})
        global fname;
        fname = strcat(pwd, "/",  file);
    end
    if nargin < 2
        plot_result = 1;
    end

    if ~isequal(fname(2),0) % check if cancel was pressed!
        %-------------------------------------------------------------------------
        % load MRXCAT parameter file (acquisition matrix, coils, ...) and update
        %-------------------------------------------------------------------------
        load( [fname(1:end) '_par.mat'] );


        %-------------------------------------------------------------------------
        % get MRXCAT data
        %-------------------------------------------------------------------------
        fid  = fopen( [fname '.cpx'] );
        data = fread(fid,inf,'float','l');
        fclose(fid);
        %-------------------------------------------------------------------------
        % get sensitivity maps
        %-------------------------------------------------------------------------
        fid  = fopen( [fname(1:end) '.sen'] );
        sen  = fread(fid,inf,'float','l');
        fclose(fid);

        %-------------------------------------------------------------------------
        % reformat to complex data (image and sensitivities)
        %-------------------------------------------------------------------------
        data = reshape( data, 2, []);
        data = data(1,:)+1i*data(2,:);
        data = reshape( data, Par.acq_matrix(1), Par.acq_matrix(2), Par.acq_matrix(3), Par.ncoils, []);
        data = permute( data, [1 2 3 5 4] );
        sen  = reshape( sen, 2, []);
        sen  = sen(1,:)+1i*sen(2,:);
        sen  = reshape( sen, Par.acq_matrix(1), Par.acq_matrix(2), Par.acq_matrix(3), Par.ncoils, []);
        sen  = permute( sen, [1 2 3 5 4] );

        %-------------------------------------------------------------------------
        % "coil-combine" data
        %-------------------------------------------------------------------------
        sos = sum(data./sen,5);
        clear sen;

        %-------------------------------------------------------------------------
        % Display 4-panel figure showing time frame, slices, coil images and (for perf) example signal-time curves
        %-------------------------------------------------------------------------
        if plot_result
        SaveData(sos);
        end


    else
        warndlg('No file loaded');
        data = [];
        fname = '';
    end
cd("..");
end
end


%=========================================================================


%=========================================================================
%	Save data as csv
%=========================================================================
function SaveData(data)
global dest;
global filename1;
global dest2;
global file;
numberimages = numel(data)/(size(data,1)*size(data,2));
imageposno = 1;
if isreal(data),
        RealData = 1;       % real data
    elseif isreal(data*i),
        RealData = -1;      % imaginary data
    else
        RealData = 0;       % complex data
    end
if ~exist(dest)
       mkdir(dest);
    end
for imageno=1:numberimages
            if RealData==1,         % real data
                csvwrite([dest "/" file num2str(imageno) '.csv'], (  (data(:,:,imageno)) )); #transpose(  (data(:,:,imageno)) ));
            elseif RealData==-1,    % imaginary data
                csvwrite([dest "/" file num2str(imageno) '.csv'], (  imag(data(:,:,imageno)) )); #transpose(  imag(data(:,:,imageno)) ))
            else                    % complex data
                csvwrite([dest "/" file num2str(imageno) '.csv'], (  abs(data(:,:,imageno)) )); #transpose(  abs(data(:,:,imageno)) ));
            end
      end
end

