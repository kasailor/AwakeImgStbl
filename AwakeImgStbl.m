%2022-11-25 AwakeImgStbl written by Kurt A. Sailor Institut Pasteur, Paris,
%France using Matlab 2018a. This code has input of multiple serial z-stack tiff image
%sequences for running a 2D cross correlation of the same image plane
%across the separate z-series. Initially individual images for trials are registered to the
%mean image for that trial. A 2D correlation is computed for each slice
%versus each other slice from a trial. The 2D correlation is normalized for
%each trial (important) to the maximum value. The user inputs a correlation
%value threshold and the z-stack is reconstructed by taking the mean of all
%the trlas above this threshold for each slice. A final registration
%is performed on the compiled stack.
clear all
%************Dependencies********************
%1. Matlab Image Processing Toolbox
%2. 'saveastiff' https://www.mathworks.com/matlabcentral/fileexchange/35684-multipage-tiff-stack 
%3. 'loadtiff' https://github.com/flatironinstitute/CaImAn-MATLAB/blob/master/loadtiff.m 
%4. Install Miji by going to Set Path, Add with Subfolders..., go
%to your Fiji directory and open scripts.
%5. Fiji plugin "multistackreg" https://github.com/miura/MultiStackRegistration 
%If using large images, in Matlab goto preferences->Java Heap Memory and increase the
%allocated memory.
%__________________________________________________________________________
%Data structure:
%rootdir folder containing multiple ZSeries from one experiment
%--ZSeries folders-PrairieView output ZSeries !Important! must have
%'ZSeries' as beginning of folder name otherwise change DirID input variable.

%__________________________________________________________________________
%Input variables:
%reRun set to 0 for first analysis. Set to 1 if want to re-do threshold.
reRun = 0;
%rootdir is directory containing 2 trials from Bruker imaging.
rootdir = 'D:\Awake structure imaging code\Data\Dataset 9\20221018_158_stretch_1';
%outputdir is where you want final processed stack saved.
outputdir = 'D:\Awake structure imaging code\Output\Dataset 9';
%DirID is unique search name for z-series folders
DirID = 'ZSeries*';
%PCR is percent to crop X-Y dimensions of images after registration for better 2d correlation.
PCR = 15;
%RegMethod set to 1 for rigid or 2 for affine registration of original data before thresholding.
RegMethod = 1;
%__________________________________________________________________________
if reRun == 0
    cd(rootdir)
    
    %Find folders starting with 'ZSeries', !!change this if data
    %structure is different.
    folders = dir(DirID);
    
    %For loading individual tiffs segregated by folder
    for k = 1:length(folders)
        cd([rootdir '\' folders(k).name]);
        files = dir('*.tif');
        for h = 1:length(files)
            stack(k,:,:,h) = imread(files(h).name);
        end
    end
    
    %Make 3D tiffs of individual slices across all trials
    cd(outputdir)
    mkdir('RawData')
    for h = 1:size(stack,4)
        for k = 1:size(stack,1)
            temp(:,:,k) = squeeze(stack(k,:,:,h));
        end
        temp = uint16(temp);
        saveastiff(temp,[outputdir '\RawData\rawSlice' num2str(h,'%03d') '.tif']);
    end
    meanstack = uint16(squeeze(mean(stack,1)));
    saveastiff(meanstack,[outputdir '\OriginalAveraged.tif']);
    
    %Make handshake math problem table for cross comparisons
    counter = 1;
    stepnumb = 0;
    for k = 1:size(stack,1)
        counter2 = size(stack,1) - k;
        for j = 1:counter2
            crosstable(1,counter) = k;
            crosstable(2,counter) = j+1+stepnumb;
            counter2 = counter2 + 1;
            counter = counter + 1;
        end
        stepnumb = stepnumb + 1;
        clear counter2
    end
    
    %Run stackreg across all trials in one stack with mean as template
    %performing rigid registration.
    Miji(false)
    cd([outputdir '\RawData\']);
    Ztrials = dir('raw*');
    mkdir('temp');
    mkdir('registered');
    imgname = 'tempstack.tif';
    FileLoadFIJI = strrep(['path=[' outputdir '\RawData\temp\tempstack.tif' ']'], '\', '\\');
    for k = 1:size(meanstack,3)
        tempstack = loadtiff([outputdir '\RawData\' Ztrials(k).name]);
        for b = 1:size(tempstack,3)
            Image2reg = cat(3,meanstack(:,:,k),tempstack(:,:,b));
            Image2reg = uint16(Image2reg);
            saveastiff(Image2reg,[outputdir '\RawData\temp\tempstack.tif']);
            MIJ.run('Open...', ['''' FileLoadFIJI '''']);
            if RegMethod == 1
                MIJ.run('MultiStackReg', ['stack_1=' imgname ' action_1=Align file_1=[] stack_2=None action_2=Ignore file_2=[]']);
            end
            if RegMethod == 2
                MIJ.run('MultiStackReg', ['stack_1=' imgname ' action_1=Align file_1=[] stack_2=None action_2=Ignore file_2=[] transformation=Affine']);
            end
            temp2reg = MIJ.getCurrentImage;
            temp2reg = uint16(temp2reg);
            temp2regtrials(:,:,b) = temp2reg(:,:,2);
            MIJ.run('Close');
        end
        saveastiff(temp2regtrials,[outputdir '\RawData\registered\' Ztrials(k).name(1:end-4) '_reg.tif']);
        disp(['Sice ' num2str(k) ' out of ' num2str(size(meanstack,3)) ' registered']);
    end
    
    %Load registered trial stacks and record 2D correlation for each slice vs
    %slice in each trial for each z-slice.
    cd([outputdir '\RawData\registered\']);
    regStacks = dir('raw*');
    Dim = ceil((PCR/100)*(size(meanstack(:,:,1))));
    xDim = Dim(1):size(meanstack,1)-Dim(1);
    yDim = Dim(2):size(meanstack,2)-Dim(2);
    for k = 1:size(meanstack,3)
        tempstackReg = loadtiff(regStacks(k).name);
        for j = 1:length(crosstable)
            crosscompareOrig(k,j) = corr2(squeeze(tempstackReg(xDim,yDim,crosstable(1,j))),squeeze(tempstackReg(xDim,yDim,crosstable(2,j))));
        end
    end
    
    %Organize each slice across z-stacks with corr values
    for slices = 1:size(stack,4)
        for trials = 1:size(stack,1)
            binmap = sum(crosstable == trials,1);
            counter = 1;
            for h = 1:length(binmap)
                if binmap(h) == 1
                    sortedcorrs(trials,counter) = crosscompareOrig(slices,h);
                    counter = counter + 1;
                end
            end
        end
        corrvalssmean(slices,:) = mean(sortedcorrs,2);
        clear sortedcorrs
    end
    
    %Normalized corr values across the trials for each slice
    for k = 1:size(corrvalssmean,1)
        m = max(corrvalssmean(k,:));
        corrvalssmeanNorm(k,:) = corrvalssmean(k,:)/m;
    end
    
    %Plot heatmap of correlations Z-stack trial versus Z-slices
    figure;imagesc(corrvalssmeanNorm);xlabel('Trial number');ylabel('Z-slice number');
    colorbar('eastoutside');title('2D crosscorrelation mean for each slice for each trial');
    cd(outputdir);
    saveas(gcf,'CorrValues.png');
    
    %Pop-up message for inputing correlation cutoff
    prompt = 'Input correlation cutoff (0-1): ';
    Corrcut = input(prompt);
    
    %Take only slices above Corrcut threshold, take mean insert into final stack.
    %If no slices above threshold at Z-level, fill slice with zeros.
    cd([outputdir '\RawData\registered\']);
    regNames = dir('raw*');
    for slices = 1:size(corrvalssmeanNorm,1)
        regstack = loadtiff(regNames(slices).name);
        for trials = 1:size(corrvalssmeanNorm,2)
            counter = 1;
            if corrvalssmeanNorm(slices,trials) > Corrcut == 1
                meanSlice(:,:,counter) = squeeze(regstack(:,:,trials));
                counter = counter + 1;
            end
        end
        if exist('meanSlice') == 1
            FinalStack(:,:,slices) = mean(meanSlice,3);
            clear meanSlice
        else FinalStack(:,:,slices) = zeros(size(regstack,2),size(regstack,3));
            clear meanSlice
        end
    end
    
    FinalStack = uint16(FinalStack);
    saveastiff(FinalStack,[outputdir '\CorrectedStack.tif']);
    figure;imagesc(corrvalssmeanNorm > Corrcut);xlabel('Trial number');ylabel('Z-slice number');
    colorbar('eastoutside');title('Slices used for mean reconstruction');
    colormap('gray');
    cd(outputdir);
    saveas(gcf,'Reconstruction map.png');
    
    FileLoadFIJI = strrep(['path=[' outputdir '\CorrectedStack.tif' ']'], '\', '\\');
    MIJ.run('Open...', ['''' FileLoadFIJI '''']);
    imgname = 'CorrectedStack.tif';
    MIJ.run('MultiStackReg', ['stack_1=' imgname ' action_1=Align file_1=[] stack_2=None action_2=Ignore file_2=[] transformation=Affine']);
    tempstackReg = MIJ.getCurrentImage;
    tempstackReg = uint16(tempstackReg);
    tempstackReg = tempstackReg(:,:,2:end);
    MIJ.run('Close');
    saveastiff(tempstackReg,[outputdir '\CorrectedStack_Registered.tif']);
    cd(outputdir);
    save('corrvalssmeanNorm.mat','corrvalssmeanNorm');
end

if reRun == 1
    cd(outputdir);
    load('corrvalssmeanNorm.mat');
    figure;imagesc(corrvalssmeanNorm);xlabel('Trial number');ylabel('Z-slice number');
    colorbar('eastoutside');title('2D crosscorrelation mean for each slice for each trial');
    
    prompt = 'Input correlation cutoff: ';
    Corrcut = input(prompt);
    
    %Take only slices above Corrcut threshold, take mean, if no slices above
    %threshold, fill with zeros
    cd([outputdir '\RawData\registered\']);
    regNames = dir('raw*');
    for slices = 1:size(corrvalssmeanNorm,1)
        regstack = loadtiff(regNames(slices).name);
        for trials = 1:size(corrvalssmeanNorm,2)
            counter = 1;
            if corrvalssmeanNorm(slices,trials) > Corrcut == 1
                meanSlice(:,:,counter) = squeeze(regstack(:,:,trials));
                counter = counter + 1;
            end
        end
        if exist('meanSlice') == 1
            FinalStack(:,:,slices) = mean(meanSlice,3);
            clear meanSlice
        else FinalStack(:,:,slices) = zeros(size(regstack,2),size(regstack,3));
            clear meanSlice
        end
    end
    
    FinalStack = uint16(FinalStack);
    saveastiff(FinalStack,[outputdir '\CorrectedStack.tif']);
    figure;imagesc(corrvalssmeanNorm > Corrcut);xlabel('Trial number');ylabel('Z-slice number');
    colorbar('eastoutside');title('Slices used for mean reconstruction');
    colormap('gray');
    saveas(gcf,'Reconstruction map.png');
    
    cd(outputdir);
    FileLoadFIJI = strrep(['path=[' outputdir '\CorrectedStack.tif' ']'], '\', '\\');
    MIJ.run('Open...', ['''' FileLoadFIJI '''']);
    imgname = 'CorrectedStack.tif';
    MIJ.run('MultiStackReg', ['stack_1=' imgname ' action_1=Align file_1=[] stack_2=None action_2=Ignore file_2=[] transformation=Affine']);
    tempstackReg = MIJ.getCurrentImage;
    tempstackReg = uint16(tempstackReg);
    tempstackReg = tempstackReg(:,:,2:end);
    MIJ.run('Close');
    saveastiff(tempstackReg,[outputdir '\CorrectedStack_Registered.tif']);
end