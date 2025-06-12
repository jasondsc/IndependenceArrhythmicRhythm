% this code reads in the simulated time series csv values
% then imports this data into brainstorm for fooofing/ ms-specparam

% % % % % % % % % % % % % % % % % % % % % % % % %
 % delete files in brainstorm db before running again!!!
% % % % % % % % % % % % % % % % % % % % % % % % % %
close all
clear 

addpath('/Users/jason/Documents/SickKids/abagen_analysis/natsortfiles')
file2fooof= '/Users/jason/Downloads/TutorialOmega/data/Simulations/TriggerTest_SF031_20250204_01_freq/data_block001.mat';
timeseries= load('~/Documents/McGill/simulations4Giu/Template_timeseries.mat');

timeseries.Time=  0:0.0020:29.9980;
cd '/Users/jason/Documents/McGill/simulations4Giu/data/FreqEPHYS/'

simfiles= dir('sim_freq_*.csv');
simfiles=natsort({simfiles.name}); % sorted files à la mac

timeseries.F= zeros(308,length(timeseries.Time));

for i = 1001:length(simfiles)
    
    temp=readmatrix(simfiles{i});
    timeseries.F(31,:)= temp;

save(file2fooof, '-struct', 'timeseries');

% Process: Power spectrum density (Welch)
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
     
     % Process: BIC-FOOOF
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
%% organize outputs of ms-specparam from brainstorm to csv for stats
% % % % % % % % % % % % % % % % % % % % % % % % % %
cd('/Users/jason/Downloads/TutorialOmega/data/Simulations/TriggerTest_SF031_20250204_01_expon/');

colourmap=repmat([70:0.75:255],[3,1])'./255;
files_meg1=dir('./timefreq_psd_*fooof.mat');
numpar= length(files_meg1);
nchan=1;

params= repelem([0.6, 0.7, 0.8, 0.9, 1.0,1.1, 1.2, 1.3, 1.4, 1.5], 500); 
% amplitude repelem([0.3, 0.5, 0.7, 0.9, 1.0,1.2,1.4, 1.6, 1.8, 2], 500); 
% expon repelem([0.6, 0.7, 0.8, 0.9, 1.0,1.1, 1.2, 1.3, 1.4, 1.5], 500); 
% freqs repelem([6,7,8,9,10,11,12,13,14,], 500);

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
%% load matrix of peaks
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


close all
psdNo= 3818;
figure
clear g
g(1,1)=gramm('x',repmat(PsdMat.Freqs(4:151),[2,1]),'y', [log10(PSD(psdNo,4:151))-log10(fitaper_BIC(psdNo,:));log10(fitper_BIC(psdNo,:)) ], 'color',{'PSD','FOOOF'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('FOOOFed spectrum');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
g.draw();
%saveas(gcf,'/Users/jason/Documents/FOOOF_BICvsVANILLA/figures/SIMULATIONS_FOOOFed_correctedspectrum_888.pdf')

figure
clear g
g(1,1)=gramm('x',repmat(PsdMat.Freqs(4:151),[2,1]),'y', [log10(PSD(psdNo,4:151));log10(fitaper_BIC(psdNo,:)) ], 'color',{'PSD','FOOOF'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('FOOOFed spectrum');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
g.draw();
%saveas(gcf,'/Users/jason/Documents/FOOOF_BICvsVANILLA/figures/SIMULATIONS_FOOOFed_spectrum_888.pdf')


% Write all data into a nice clean CSV for stats and modelling

cd('/Users/jason/Documents/McGill/simulations4Giu/')
output= table(params',expon_BIC, offset_BIC, mse_BIC, rsq_BIC, BF_model);
writetable(output, './Simulations_frequency_output.csv')

writematrix(BIC_peak_cf, './Simulations_frequency_peaks_centerFreq.csv')
writematrix(BIC_peak_amplitude, './Simulations_frequency_peaks_amplitude.csv')
writematrix(BIC_peak_std, './Simulations_frequency_peaks_std.csv')

writematrix(fpsd_BIC, './Simulations_exponent_fitted_spectrum.csv')
writematrix(fitaper_BIC, './Simulations_exponent_fitted_aperiodicspectrum.csv')
writematrix(fitper_BIC, './Simulations_exponent_fitted_periodicspectrum.csv')
writematrix(log10(PSD(:,4:151)), './Simulations_exponent_PSD.csv')


y_norm2= 0.4*exp((-(PsdMat.Freqs(4:151) -19).^2)/(2*(5.^2)));
y_norm =0.7*exp((-(PsdMat.Freqs(4:151) -10).^2)/(2*(2.^2)));
peakdist_1= y_norm+y_norm2;

mean(peakdist_1(22:36)) % expected rhythmic alpha 0.6457

true_corrspec=[repmat(peakdist_1,[5000,1])];

corrspec_abs_error=(corr_spec_task-true_corrspec);

writematrix(corrspec_abs_error, 'Exponent_Simulations_error_signed.csv')
writematrix(corr_spec_task, 'Exponent_Simulations_corrspec.csv')

offset_error = log10(PSD(:,4)) - offset_BIC;
writematrix(offset_error, 'Exponent_Simulations__offset_error_signed.csv')

apert_1= log10(PSD(1:500,4)) - log10(PsdMat.Freqs(4:151).^0.6);
apert_2= log10(PSD(501:1000,4)) - log10(PsdMat.Freqs(4:151).^0.7);
apert_3= log10(PSD(1001:1500,4)) - log10(PsdMat.Freqs(4:151).^0.8);
apert_4= log10(PSD(1501:2000,4)) - log10(PsdMat.Freqs(4:151).^0.9);
apert_5= log10(PSD(2001:2500,4)) - log10(PsdMat.Freqs(4:151).^1.0);
apert_6= log10(PSD(2501:3000,4)) - log10(PsdMat.Freqs(4:151).^1.1);
apert_7= log10(PSD(3001:3500,4)) - log10(PsdMat.Freqs(4:151).^1.2);
apert_8= log10(PSD(3501:4000,4)) - log10(PsdMat.Freqs(4:151).^1.3);
apert_9= log10(PSD(4001:4500,4)) - log10(PsdMat.Freqs(4:151).^1.4);
apert_10= log10(PSD(4501:5000,4)) - log10(PsdMat.Freqs(4:151).^1.5);

true_aperiodic=[apert_1;apert_2;...
    apert_3;apert_4;...
    apert_5;apert_6;...
    apert_7;apert_8;...
    apert_9;apert_10];


aper_abs_error=(log10(fitaper_BIC)-true_aperiodic);
writematrix(aper_abs_error, 'Exponent_Simulations__aperiodic_error_signed.csv')

% Frequency simulations 
y_norm2= 0.4*exp((-(PsdMat.Freqs(4:151) -19).^2)/(2*(5.^2)));

y_norm =0.7*exp((-(PsdMat.Freqs(4:151) -6).^2)/(2*(2.^2)));
peakdist_1= y_norm+y_norm2;

y_norm =0.7*exp((-(PsdMat.Freqs(4:151) -7).^2)/(2*(2.^2)));
peakdist_2= y_norm+y_norm2;

y_norm =0.7*exp((-(PsdMat.Freqs(4:151) -8).^2)/(2*(2.^2)));
peakdist_3= y_norm+y_norm2;

y_norm =0.7*exp((-(PsdMat.Freqs(4:151) -9).^2)/(2*(2.^2)));
peakdist_4= y_norm+y_norm2;

y_norm =0.7*exp((-(PsdMat.Freqs(4:151) -10).^2)/(2*(2.^2)));
peakdist_5= y_norm+y_norm2;

y_norm =0.7*exp((-(PsdMat.Freqs(4:151) -11).^2)/(2*(2.^2)));
peakdist_6= y_norm+y_norm2;

y_norm =0.7*exp((-(PsdMat.Freqs(4:151) -12).^2)/(2*(2.^2)));
peakdist_7= y_norm+y_norm2;

y_norm =0.7*exp((-(PsdMat.Freqs(4:151) -13).^2)/(2*(2.^2)));
peakdist_8= y_norm+y_norm2;

y_norm =0.7*exp((-(PsdMat.Freqs(4:151) -14).^2)/(2*(2.^2)));
peakdist_9= y_norm+y_norm2;

true_corrspec=[repmat(peakdist_1,[500,1]);repmat(peakdist_2,[500,1]);...
    repmat(peakdist_3,[500,1]);repmat(peakdist_4,[500,1]);...
    repmat(peakdist_5,[500,1]);repmat(peakdist_6,[500,1]);repmat(peakdist_7,[500,1]);repmat(peakdist_8,[500,1]);repmat(peakdist_9,[500,1]);];

corrspec_abs_error=(corr_spec_task-true_corrspec);

writematrix(corrspec_abs_error, 'Frequency_Simulations_error_signed.csv')


% Amplitude simulations
y_norm2= 0.4*exp((-(PsdMat.Freqs(4:151) -19).^2)/(2*(5.^2)));

y_norm =0.3*exp((-(PsdMat.Freqs(4:151) -10).^2)/(2*(2.^2)));
peakdist_3= y_norm+y_norm2;

y_norm =0.5*exp((-(PsdMat.Freqs(4:151) -10).^2)/(2*(2.^2)));
peakdist_4= y_norm+y_norm2;

y_norm =0.7*exp((-(PsdMat.Freqs(4:151) -10).^2)/(2*(2.^2)));
peakdist_5= y_norm+y_norm2;

y_norm =0.9*exp((-(PsdMat.Freqs(4:151) -10).^2)/(2*(2.^2)));
peakdist_6= y_norm+y_norm2;

y_norm =1.0*exp((-(PsdMat.Freqs(4:151) -10).^2)/(2*(2.^2)));
peakdist_7= y_norm+y_norm2;

y_norm =1.2*exp((-(PsdMat.Freqs(4:151) -10).^2)/(2*(2.^2)));
peakdist_8= y_norm+y_norm2;

y_norm =1.4*exp((-(PsdMat.Freqs(4:151) -10).^2)/(2*(2.^2)));
peakdist_9= y_norm+y_norm2;

y_norm =1.6*exp((-(PsdMat.Freqs(4:151) -10).^2)/(2*(2.^2)));
peakdist_10= y_norm+y_norm2;

y_norm =1.8*exp((-(PsdMat.Freqs(4:151) -10).^2)/(2*(2.^2)));
peakdist_11= y_norm+y_norm2;

y_norm =2.0*exp((-(PsdMat.Freqs(4:151) -10).^2)/(2*(2.^2)));
peakdist_12= y_norm+y_norm2;

true_corrspec=[repmat(peakdist_3,[500,1]);repmat(peakdist_4,[500,1]);...
    repmat(peakdist_5,[500,1]);repmat(peakdist_6,[500,1]);repmat(peakdist_7,[500,1]);repmat(peakdist_8,[500,1]);...
    repmat(peakdist_9,[500,1]);repmat(peakdist_10,[500,1]);repmat(peakdist_11,[500,1]);repmat(peakdist_12,[500,1])];

corrspec_abs_error=(corr_spec_task-true_corrspec);

writematrix(corrspec_abs_error, 'Amplitude_Simulations_error_signed.csv')

