 close all; clc; clear;

% Define paths
ospreyPath = 'C:/Users/ishad/Documents/MATLAB/osprey'; % Osprey toolbox
spmPath = 'C:/Users/ishad/Documents/MATLAB/spm'; % SPM toolbox
projectPath = 'C:/Users/ishad/Downloads/MRS_Data/Big_Gaba/GE'; % Root folder for data

% Add & verify SPM
addpath(spmPath); addpath(genpath(spmPath));
if exist('spm_check_version','file')~=2
    error('SPM not on path; fix spmPath.');
end

% Add Osprey
addpath(genpath(ospreyPath));

% Find subject folders
subs = dir(fullfile(projectPath,'G*_S*'));
subs = subs([subs.isdir]);
nSub = numel(subs);

% Initialize table to store QC metrics
logTbl = table( ...
    cell(nSub,1), cell(nSub,1), nan(nSub,1), nan(nSub,1), ...
    nan(nSub,1), nan(nSub,1), nan(nSub,1), ...
  'VariableNames',{'Subject','File','SNR','FWHM','FreqShift','ResidualWater','RelativeResidual'} );
row = 0;

% Loop over subjects
for i = 1:nSub
    subj = subs(i).name;
    subjPath = fullfile(projectPath, subj);
    fprintf(i, nSub, subj);

    % Find all the [.7] files
    D = dir(fullfile(subjPath,'*.7'));
    if isempty(D)
        warning(' No .7 in %s → skipping\n', subj);
        continue;
    end
    idx = find(contains({D.name},'PRESS'),1); % Prefer press (they are all press)
    if isempty(idx), idx=1; end
    filename = D(idx).name;
    full7    = fullfile(subjPath, filename);
    fprintf(' Found: %s\n', filename);

    % Create output folder
    outFld = fullfile(projectPath,'OspreyOutput_GE_Gaba',subj);
    if ~exist(outFld,'dir'), mkdir(outFld); end

    % Create temporary job file
    jobFile = fullfile(subjPath,'job_temp.m');
    fid = fopen(jobFile,'w');
    fprintf(fid, 'files = {''%s''};\n', full7);
    fprintf(fid, 'files_ref = {};\nfiles_w = {};\nfiles_nii = {};\n\n'); % No water, reference or NIfTI
    fprintf(fid, 'outputFolder = ''%s'';\n\n', outFld);
    fprintf(fid, 'seqType = ''unedited'';\n'); % Unedited data
    fprintf(fid, 'editTarget = {''none''};\n');
    fprintf(fid, 'MultiVoxel = ''SVS'';\n\n'); % Single voxel
    fprintf(fid, 'opts = struct();\n');

    % Fitting settings
    fprintf(fid, 'opts.fit = struct();\n');
    fprintf(fid, 'opts.fit.method = ''Osprey'';\n');
    fprintf(fid, 'opts.fit.style = ''Concatenated'';\n');
    fprintf(fid, 'opts.fit.range = [0.5 4.2];\n');   % Fitting range (ppm)
    fprintf(fid, 'opts.fit.rangeWater = [2.0 7.4];\n');  % Water fitting range
    fprintf(fid, 'opts.fit.bLineKnotSpace = 0.4;\n');
    fprintf(fid, 'opts.fit.fitMM = 1;\n');  % Include macromolecules
    fprintf(fid, 'opts.fit.coMM3 = ''none'';\n\n');

    % Processing settings
    fprintf(fid, 'opts.SpecReg = ''RobSpecReg'';\n');
    fprintf(fid, 'opts.SubSpecAlignment = ''L2Norm'';\n\n');  % Use L2 norm for alignment

    % Output settings
    fprintf(fid, 'opts.saveProc = 1;   %% save processed.mat\n');
    fprintf(fid, 'opts.savePDF = 1;   %% PDF figures\n');
    fprintf(fid, 'opts.saveLCM = 1;   %% LCM .RAW exports\n');
    fprintf(fid, 'opts.savejMRUI = 0;\n');
    fprintf(fid, 'opts.saveVendor = 0;\n');
    fprintf(fid, 'opts.saveNII = 0;\n\n');

    % Export fitting parameters
    fprintf(fid, 'opts.exportParams.flag = 1;\n');
    fprintf(fid, 'opts.exportParams.path = ''%s'';\n', fullfile(outFld,'exportParams'));
    fclose(fid);

    % Run Osprey pipeline
    try
        fprintf(' ▶ OspreyJob\n');  % Load job file
        MRSCont = OspreyJob(jobFile,0,'11');

        fprintf(' ▶ OspreyLoad\n');  % Load data
        MRSCont = OspreyLoad(MRSCont);

        fprintf(' ▶ OspreyProcess\n');  % Preprocessing
        MRSCont = OspreyProcess(MRSCont);

        fprintf(' ▶ OspreyFit\n');  % Fit spectra
        MRSCont = OspreyFit(MRSCont);

        % Save containers (manually to make sure)
        save(fullfile(outFld,'processed.mat'),'MRSCont','-v7.3');
        fitMat = MRSCont.fit;
        save(fullfile(outFld,'fit.mat'),'fitMat');

        % Get QC metrics and store
        row = row+1;
        logTbl.Subject{row} = subj;
        logTbl.File{row} = filename;
        logTbl.SNR(row) = MRSCont.QM.SNR.metab(1);
        logTbl.FWHM(row) = MRSCont.QM.FWHM.metab(1);
        logTbl.FreqShift(row) = MRSCont.QM.freqShift.metab(1);
        logTbl.ResidualWater(row) = MRSCont.QM.res_water_amp.metab(1);
        if isfield(MRSCont.QM,'relAmpl') && isfield(MRSCont.QM.relAmpl,'metab_A')
            logTbl.RelativeResidual(row) = MRSCont.QM.relAmpl.metab_A(1);
        else
            logTbl.RelativeResidual(row) = NaN;
        end

        fprintf('%s done\n\n', subj);
    catch ME
        warning('%s FAILED: %s\n\n', subj, ME.message);
    end

    delete(jobFile);
end

% Save QC logs to CSV
logTbl = logTbl(1:row,:);
outCSV = fullfile(projectPath,'OspreyOutput_GE_Gaba','processing_log.csv');
writetable(logTbl,outCSV);
fprintf('Finished. QC log → %s\n', outCSV);