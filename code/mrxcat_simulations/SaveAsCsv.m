%=========================================================================
%    Function to save all simulation .cpx files in subfolders of
%    the ../../data/mrxcat_simulations folder to csv files
%    This file is derived from DisplayMRXCAT.m of the MRXCAT software
%=========================================================================
function [data, fname] = SaveAsCsv()
  cd('../../data/mrxcat_simulations');
  foldernames = glob ("snr*");

  for i = 1:numel(foldernames);
    [~, foldername] = fileparts (foldernames{i});
    cd(foldername);
    files = glob('*.cpx');
    [~, file] = fileparts (files{1});
    count = size(files)(1);
    dest = strcat(pwd, "/csvs");
  
    fname = strcat(pwd, "/",  file);
  
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
  
    SaveData(sos, dest, file);

    cd("..");
  end
end



%=========================================================================
%	Save data as csv
%=========================================================================
function SaveData(data, dest, file)
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

