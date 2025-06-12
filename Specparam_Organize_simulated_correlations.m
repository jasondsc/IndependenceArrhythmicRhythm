% this code reads in the simulated time series csv values
% then imports this data into brainstorm for fooofing/ ms-specparam

% this script is specifically to organize the simulations associated with
% Figure 4

% % % % % % % % % % % % % % % % % % % % % % % % %
 % delete files in brainstorm db before running again!!!
% % % % % % % % % % % % % % % % % % % % % % % % % %
close all
clear 

expon_BIC_final= [];
offset_BIC_final= [];
mse_BIC_final= [];
rsq_BIC_final= [];
BF_model_final= [];
params= [];
corr_spec_task_final= [];
corr_spec_task_linear_final= [];
exponsim= [];
age = [];
BIC_peak_cf_final = [];
BIC_peak_amplitude_final = [];
BIC_peak_std_final = [];

for icorr=  [-0.0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9]
addpath('/Users/jason/Documents/SickKids/abagen_analysis/natsortfiles')

newfolder= ['/Users/jason/Downloads/TutorialOmega/data/Simulations/TriggerTest_SF031_20250204_effectBETAexpon_',num2str(icorr, '%.1f')];
templatefolder= '/Users/jason/Downloads/TutorialOmega/data/Simulations/TriggerTest_SF031_20250204_corrmedium/';
mkdir(newfolder)
copyfile(templatefolder, newfolder)
iprotocol= bst_get('Protocol', 'TutorialOmega');
db_reload_database(iprotocol)

file2fooof= [newfolder '/data_block001.mat'];
timeseries= load('~/Documents/McGill/simulations4Giu/Template_timeseries.mat');

timeseries.Time=  0:0.0020:29.9980;
cd(['/Users/jason/Documents/McGill/simulations4Giu/data/exponbetacorr/sim_', num2str(-1*icorr, '%.1f')])

simfiles= dir('sim_amp_*.csv');
simfiles=natsort({simfiles.name}); % sorted files à la mac

timeseries.F= zeros(308,length(timeseries.Time));

for i = 1:length(simfiles)
    
    temp=readmatrix(simfiles{i});
    timeseries.F(31,:)= temp;

save(file2fooof, '-struct', 'timeseries');

%Process: Power spectrum density (Welch)
sFiles = bst_process('CallProcess', 'process_psd', file2fooof, [], ...
    'timewindow',  [], ...
    'win_length',  3, ...
    'win_overlap', 50, ...
    'units',       'physical', ...  % Physical: U2/Hz
    'sensortypes', 'MLC11', ...
    'win_std',     0, ...
    'edit',        struct(...
         'Comment',         'Power', ...
         'TimeBands',       [], ...
         'Freqs',           [], ...
         'ClusterFuncTime', 'none', ...
         'Measure',         'power', ...
         'Output',          'all', ...
         'SaveKernel',      0));
     
     %Process: BIC-FOOOF
sFiles = bst_process('CallProcess', 'process_fooof_BIC', sFiles, [], ...
    'implementation', 'matlab', ...  % Matlab
    'freqrange',      [1, 50], ...
    'method',         'negloglike', ...  % Model selection
    'peaktype',       'gaussian', ...  % Gaussian
    'peakwidth',      [1, 24], ...
    'maxpeaks',       6, ...
    'minpeakheight',  1 / 10, ...
    'proxthresh',     2, ...
    'apermode',       'fixed', ...  % Fixed
    'guessweight',    'none', ...  % None
    'sorttype',       'param', ...  % Peak parameters
    'sortparam',      'frequency', ...  % Frequency
    'sortbands',      {'delta', '2, 4'; 'theta', '5, 7'; 'alpha', '8, 12'; 'beta', '15, 29'; 'gamma1', '30, 59'; 'gamma2', '60, 90'});

    
end


% % % % % % % % % % % % % % % % % % % % % % % % % %
% organize outputs of ms-specparam from brainstorm to csv for stats
% % % % % % % % % % % % % % % % % % % % % % % % % %
cd(newfolder);

colourmap=repmat([70:0.75:255],[3,1])'./255;
files_meg1=dir('./timefreq_psd_*fooof.mat');
numpar= length(files_meg1);
nchan=1;

PSD= zeros([numpar, 751]);
BIC= cell([numpar,1]);
mse_BIC= zeros([numpar,1]);
rsq_BIC=zeros([numpar,1]);
BF_model= ones(numpar, nchan);

count=1;
for i=1:numpar
    
    PsdMat=load(strcat(files_meg1(i).folder,'/',files_meg1(i).name));
    PSD(count,:,:)=squeeze(PsdMat.TF);
    BIC{count,1}=PsdMat.Options.FOOOF;
    count=count+1;
    
end 


% MODEL FITS 
for i=1:numpar
   
    mse_BIC(i)=  [BIC{i,1}.data.error];
    rsq_BIC(i)=  [BIC{i,1}.data.r_squared];
    fmse_BIC(i,:)=  squeeze(log10(PSD(i,4:151 )))-log10(BIC{i,1}.data(1).fooofed_spectrum);
    temp=[BIC{i,1}.data(1).models];
    [m,index]= min([temp.BIC]);
    BF_model(i,1)= temp(index).BF;

end


figure('Position', [100 100 800 400])
clear g

g(1,1)=gramm('x',[rsq_BIC],'color',reshape(repmat(repelem({'BIC'}, numpar), [1,nchan]), [numpar,nchan]));
g(1,2)=copy(g(1));
%Raw data as raster plot
g(1,1).geom_raster();g(1,1).set_color_options('map','lightness_range');
g(1,1).set_title('R^2');

%Histogram
g(1,2).stat_bin('nbins',8);
g(1,2).set_title('R^2');
g(1,2).set_color_options('map','lightness_range');
g.draw();


figure
clear g

g(1,1)=gramm('x',repmat(PsdMat.Freqs(4:151),[nchan,1]),'y', squeeze(mean(fmse_BIC,1)), 'color', categorical(repelem(1:nchan, 148)));
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Frequency-wise error Vanilla' );
g(1,1).no_legend();
g(1,1).set_color_options('map',colourmap);
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Error');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();
% saveas(gcf,'/Users/jason/Documents/FOOOF_BICvsVANILLA/figures/REVIEW_frequency_wise_error_vanilla.pdf')

% Foofed spectrum
for i=1:numpar
    
    temp1=BIC{i,1}.data;
    fpsd_BIC(i,:)=  temp1(1).fooofed_spectrum(1:148);
    temp=[BIC{i,1}.data.aperiodic_params]';
    offset_BIC(i,:)= temp(:,1);
    expon_BIC(i,:)= temp(:,2);
end

BIC_psd=squeeze(mean(log10(fpsd_BIC),1));

figure
clear g
g(1,1)=gramm('x',repmat(PsdMat.Freqs(4:151),[2,1]),'y', [mean(log10(PSD(:,4:151)),1); BIC_psd], 'color',{'PSD','FOOOF'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('FOOOFed spectrum');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
g.draw();
% saveas(gcf,'/Users/jason/Documents/FOOOF_BICvsVANILLA/figures/REVIEW_FOOOFed_spectrum.pdf')


% Collect fits of aperiodic and periodic 

for i=1:numpar
    
    temp12=reshape([BIC{i,1}.data.ap_fit], [size([BIC{i,1}.data.ap_fit],2)/nchan,nchan]);
    fitaper_BIC(i,:,:)=  temp12(1:148,:);

    temp12=reshape([BIC{i,1}.data.peak_fit], [size([BIC{i,1}.data.peak_fit],2)/nchan,nchan]);
    fitper_BIC(i,:,:)=  temp12(1:148,:);

end


% organize peak params 
% load matrix of peaks
% make the matrix 606 by 148, by 6 max num peaks 

BIC_peak_cf= NaN(numpar,6);
BIC_peak_amplitude= NaN(numpar,6);
BIC_peak_std= NaN(numpar,6);

for i=1:numpar
       
       tempp=  [BIC{i,1}.data.peak_params];
       BIC_peak_cf(i, 1:length(tempp(:,1)))= tempp(:,1);
       BIC_peak_amplitude(i,1:length(tempp(:,1)))= tempp(:,2);
       BIC_peak_std(i, 1:length(tempp(:,1)))= tempp(:,3);

end

% aperiodic corrected 
corr_spec_task=log10(PSD(:,4:151))-log10(fitaper_BIC);

figure
clear g

g(1,1)=gramm('x',PsdMat.Freqs(4:151),'y', squeeze(mean(corr_spec_task,1)), 'color', categorical(repelem(1:nchan, 148)));
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Aperiodic corrected spectra Vanilla');
g(1,1).no_legend();
g(1,1).set_color_options('map',colourmap);
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
g.draw();
%saveas(gcf,'/Users/jason/Documents/FOOOF_BICvsVANILLA/figures/REVIEW_Vanilla_aperiodic_corrected_spectra.pdf')
corr_spec_task_linear=PSD(:,4:151)-fitaper_BIC;

% append extracted values to matrix 
% will save these later outside for loop
expon_BIC_final= [expon_BIC_final; expon_BIC];
offset_BIC_final= [offset_BIC_final; offset_BIC];
mse_BIC_final= [mse_BIC_final; mse_BIC];
rsq_BIC_final= [rsq_BIC_final; rsq_BIC];
BF_model_final= [BF_model_final; BF_model];
params= [params; repelem( [icorr], [length(simfiles)])'];
corr_spec_task_final= [corr_spec_task_final; corr_spec_task];
corr_spec_task_linear_final= [corr_spec_task_linear_final; corr_spec_task_linear];

BIC_peak_cf_final = [BIC_peak_cf_final; BIC_peak_cf];
BIC_peak_amplitude_final = [BIC_peak_amplitude_final; BIC_peak_amplitude];
BIC_peak_std_final = [BIC_peak_std_final; BIC_peak_std];

% read in correlational data to make one nice big table
temp=readtable(['/Users/jason/Documents/McGill/simulations4Giu/data/exponbetacorr/simulated_correlations_' num2str(-1*icorr, '%.1f') '_sample.csv']);
exponsim= [exponsim; table2array(temp(:,2))];
age= [age; table2array(temp(:,3))];

% sanity check the correlations between alpha amplitudes
corr(mean(corr_spec_task(:,22:36),2), table2array(temp(:,2)))

corr(mean(corr_spec_task_linear(:,22:36),2), table2array(temp(:,2)))

WaitSecs(10)
close all
end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% % Write all data into a nice clean CSV for stats and modelling
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

cd('/Users/jason/Documents/McGill/simulations4Giu/')

output= table(params,expon_BIC_final, offset_BIC_final, mse_BIC_final, rsq_BIC_final, BF_model_final, age,exponsim );
writetable(output, './Simulations_of_correlations_exponBETAcorr_agerelationship_output.csv')

writematrix(BIC_peak_cf_final, './Simulations_of_correlations_exponBETAcorr_agerelationship_peaks_centerFreq.csv')
writematrix(BIC_peak_amplitude_final, './Simulations_of_correlations_exponBETAcorr_agerelationship_peaks_amplitude.csv')
writematrix(BIC_peak_std_final, './Simulations_of_correlations_exponBETAcorr_agerelationship_peaks_std.csv')

writematrix(corr_spec_task_final, 'Simulations_of_correlations_exponBETAcorr_agerelationship_corrspec.csv')
writematrix(corr_spec_task_linear_final, 'Simulations_of_correlations_exponBETAcorr_agerelationshipy_corrspec_linear.csv')
