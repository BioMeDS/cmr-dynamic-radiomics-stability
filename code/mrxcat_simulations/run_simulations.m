% Run this script in the main MRXCAT directory
mkdir("simulations");

simFiles = dir(fullfile('@MRXCAT_CMR_CINE/simulations',"*.m"));
for i = 1:length(simFiles)
    f = simFiles(i).name;
    f = f(1:length(f)-2);
    filePath = fullfile(simFiles(i).folder, simFiles(i).name)
    copyfile(filePath, "@MRXCAT_CMR_CINE/CINEpar.m", "f");
    % force reload of CINEpar file
    clear CINEpar;
    for rep = 1:5
        outdir = fullfile("simulations",f + "_" + rep);
        disp(outdir);
        indir = "cine_breathhold";
        mkdir(outdir);
        filename = convertStringsToChars(fullfile(indir,'cine_act_1.bin'));
        MRXCAT_CMR_CINE(filename);
        movefile(fullfile(indir,"cine_1x1*"), outdir);
        copyfile(filePath, fullfile(outdir,"CINEpar.m"));
    end
end
