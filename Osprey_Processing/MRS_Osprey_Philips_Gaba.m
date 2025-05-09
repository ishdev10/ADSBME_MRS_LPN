close all; clc; clear;

% Define paths
ospreyPath = 'C:/Users/ishad/Documents/MATLAB/osprey'; % Osprey toolbox
spmPath = 'C:/Users/ishad/Documents/MATLAB/spm'; % SPM toolbox
projectPath = 'C:/Users/ishad/Downloads/MRS_Data/Big_Gaba/Philips'; % Root folder for data

% Add & verify SPM
addpath(spmPath); addpath(genpath(spmPath));
if exist('spm_check_version','file')~=2
    error('SPM not on path; fix spmPath.');
end

% Add Osprey
addpath(genpath(ospreyPath));

% Find subject folders
subs = dir(fullfile(projectPath,'P*_S*'));
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
    fprintf('--- (%d/%d) %s ---\n', i, nSub, subj);

    % Find all the [.SDAT] files
    D = dir(fullfile(subjPath,'*PRESS_35_act.SDAT'));
    if isempty(D)
        warning(' No PRESS_act.SDAT in %s → skipping\n', subj);
        continue;
    end
    filename = D(1).name;
    fullSDAT = fullfile(subjPath, filename);
    fprintf(' Found: %s\n', filename);
    
    % Find corresponding [.SPAR] file
    fullSPAR = strrep(fullSDAT, '.SPAR', '.SPAR');
    if ~exist(fullSPAR, 'file')
        warning(' Corresponding SPAR file not found for %s → skipping\n', filename);
        continue;
    end

   % Create output folder
    outFld = fullfile(projectPath,'OspreyOutput_Philips_Gaba',subj);
    if ~exist(outFld,'dir'), mkdir(outFld); end

    % Create temporary job file
    jobFile = fullfile(subjPath,'job_temp.m');
    fid = fopen(jobFile,'w');
    fprintf(fid, 'files = {''%s''};\n', fullSDAT);
    
    % Add water reference if available
    refSDAT = strrep(fullSDAT, 'act.SDAT', 'ref.SDAT');
    if exist(refSDAT, 'file')
        fprintf(fid, 'files_ref = {''%s''};\n', refSDAT);
    else
        fprintf(fid, 'files_ref = {};\n');
    end
    
    fprintf(fid, 'files_w = {};\n');
    fprintf(fid, 'files_nii = {};\n\n');
    fprintf(fid, 'outputFolder = ''%s'';\n\n', outFld);
    fprintf(fid, 'seqType = ''unedited'';\n');
    fprintf(fid, 'editTarget = {''none''};\n');
    fprintf(fid, 'MultiVoxel = ''SVS'';\n\n');
    fprintf(fid, 'opts = struct();\n');
    
    % Fitting settings
    fprintf(fid, 'opts.fit = struct();\n');
    fprintf(fid, 'opts.fit.method = ''Osprey'';\n');
    fprintf(fid, 'opts.fit.style = ''Concatenated'';\n');
    fprintf(fid, 'opts.fit.range = [0.5 4.2];\n');
    fprintf(fid, 'opts.fit.rangeWater = [2.0 7.4];\n');
    fprintf(fid, 'opts.fit.bLineKnotSpace = 0.4;\n');
    fprintf(fid, 'opts.fit.fitMM = 1;\n');
    fprintf(fid, 'opts.fit.coMM3 = ''none'';\n\n');
    
    % Processing settings
    fprintf(fid, 'opts.SpecReg = ''RobSpecReg'';\n');
    fprintf(fid, 'opts.SubSpecAlignment = ''L2Norm'';\n\n');
    
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
        MRSCont = OspreyJob(jobFile, 0, '11');
        
        fprintf(' ▶ OspreyLoad\n');  % Load data
        MRSCont = OspreyLoad(MRSCont);

        fprintf(' ▶ OspreyProcess\n');  % Preprocessing
        MRSCont = OspreyProcess(MRSCont);

        fprintf(' ▶ OspreyFit\n');  % Fit spectra
        MRSCont = OspreyFit(MRSCont);

        % Save containers (manually to make sure)
        save(fullfile(outFld,'processed.mat'),'MRSCont','-v7.3');
        fitMat = MRSCont.fit; %#ok<NASGU>
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
        fprintf(2, 'FAILED: %s\n', subj, ME.message);
    end

    delete(jobFile);
end

% Save QC logs to CSV
logTbl = logTbl(1:row,:);
outCSV = fullfile(projectPath,'OspreyOutput_Philips_Gaba','processing_log.csv');
writetable(logTbl,outCSV);
fprintf('Finished. QC log → %s\n', outCSV);