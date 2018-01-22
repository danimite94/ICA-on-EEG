% session7_1_2.m
% October 2015

clear; close all; clc; 
addpath ./robustica_package/

%electrode setup (channel names)
CHANnames = {'O2','O1', 'P4', 'Pz', 'P3', 'C4', 'Cz', 'C3', 'Fz', 'F4', 'F3', 'Fpz'};

%% Data
load EEG.mat; 

X=data;
figure('name','Plot of all EEG Sources')
title('Plot of All EEG Sources')
eegplot_simple(X,CHANnames)
title('Plot of All EEG Sources')


%% P300 grand average I
% Cut P300 trials from the Fpz (the last) channel
% Take 1 second (128 samples) of data starting at all event indices
% Plot all trials. Plot the grand average in bold.

fpz = data(end,:);

%fpz length is about 35*128 samples. So we need 35 fpz segments. We will
%create a matrix where each row represents a segment.
fpz_segments=[];

%% in here we will store the segments mentioned above

for i=1:35
    
    begin=i*128-127;
    final=i*128;
    
    fpz_segments=[fpz_segments;fpz(begin:final)]; 
    
end
%% in here we will compute the average vector

pfz_averaged=linspace(0,0,128);

for i=1:35 
     
    pfz_averaged=pfz_averaged+fpz_segments(i,:);
    
end

pfz_averaged=pfz_averaged/35;

figure('name','Synchronized Averaging')
for i=1:35
    plot(fpz_segments(i,:))
    hold on
end
plot(pfz_averaged,'r','LineWidth',3);
title('Grand Averaging on the Original FPZ source')
xlabel('Samples')
ylabel('Amplitude')
axis tight


%% ICA Algorithms
COMPnames = arrayfun(@(x) ['Cmp ' num2str(x)], 1:size(X,1),'UniformOutput',false); %array of component names

% ROBUST ICA
tol = 1e-3;     % termination threshold parameter
max_it = 1000;   % maximum number of iterations per independent component
[y_robustICA, A_robustICA] = robustica(X, [], tol, max_it, 1, 'r', 0, [], 0);   % regression-based deflation

figure('name','EEG ROBUST ICA'); eegplot_simple(y_robustICA, COMPnames);
title('Sources estimated using ROBUST ICA')
figure('name','TopoPlot ROBUST ICA'); topoPlot(A_robustICA, COMPnames);
title('TOPOGRAPHY estimated using ROBUST ICA')
% I_robustICA = input('Please enter the component number(s)...\n>');

% in here there are the components from which we think there are related to eye blink movement
I_robustICA=[1 3 7 12];
clc;


% SOBI
[A_SOBI,y_SOBI] = acsobiro(X,size(X,1),100);

figure('name','EEG SOBI'); eegplot_simple(y_SOBI, COMPnames);
title('Sources estimated using ROBUST ICA')
figure('name','TopoPlot SOBI'); topoPlot(A_SOBI, COMPnames);
title('TOPOGRAPHY estimated using ROBUST ICA')

% I_SOBI = input('Please enter the component number(s)...\n>');

% in here there are the components from which we think there are related to eye blink movement
I_SOBI=[1 6 8];


 

%% Artefact removal 
% based on the temporal and spatial information obtained by ICA, select
% which components to remove for each ICA algorithm. I.e. by using the I_robustICA 
% defined earlier. Adjust the mixing matrices accordingly and use them to 
% construct a version of the EEG without artifacts. 
% Plot the clean EEG for each algorithm and compare to the first plot in this
% exercise.

%Robust ICA

%we are making value 0 the components which are responsible for the eye
%blink movement. So we make modify the mixing matrix in order to those
%components to have weight 0.
for (i=1:length(I_robustICA))
    A_robustICA(:,I_robustICA(i))=0;
end

reconstructed_robustICA=A_robustICA*y_robustICA;

figure('name','EEG reconstructed ROBUST ICA'); eegplot_simple(reconstructed_robustICA, COMPnames);
figure('name','TopoPlot reconstructed ROBUST ICA'); topoPlot(A_robustICA, COMPnames);


%SOBI

%we are making value 0 the components which are responsible for the eye
%blink movement. So we make modify the mixing matrix in order to those
%components to have weight 0.

for (i=1:length(I_SOBI))
    A_SOBI(:,I_SOBI(i))=0;
end

reconstructed_SOBI=A_SOBI*y_SOBI;

figure('name','EEG reconstructed SOBI'); eegplot_simple(reconstructed_SOBI, COMPnames);
title('EEG reconstructed SOBI')
figure('name','TopoPlot reconstructed using SOBI'); topoPlot(A_SOBI, COMPnames);
title('Topography of the EEG reconstructed using SOBI')

%% P300 grand average II
% Repeat the P300 synchronized averaging, but this time on the Fpz channel 
% of the EEG from which artifacts have been removed (for all ICA algorithms). 
% Plot the trials and the grand average in bold. Compare the original P300 to 
% the ones obtained after artifact removal. 


%Grand average for the EEG reconstructed with Robust ICA

fpzICA = reconstructed_robustICA(end,:);

%fpz length is about 35*128 samples. So we need 35 fpz segments. We will
%create a matrix where each row represents a segment.
fpz_segmentsICA=[];


for i=1:35
    
    begin=i*128-127;
    final=i*128;
    
    fpz_segmentsICA=[fpz_segmentsICA;fpzICA(begin:final)]; 
    
end

pfz_averagedICA=linspace(0,0,128);

for i=1:35 
     
    pfz_averagedICA=pfz_averagedICA+fpz_segmentsICA(i,:);
    
end

pfz_averagedICA=pfz_averagedICA/35;

figure('name','Synchronized Averaging - After reconstruction with ROBUST ICA')
for i=1:35
    plot(fpz_segmentsICA(i,:))
    hold on
end
plot(pfz_averagedICA,'r','LineWidth',3);
title('Synchronized Averaging - After reconstruction with ROBUST ICA')
xlabel('Samples')
ylabel('Amplitude')
axis tight


%SOBI

%Grand average for the EEG reconstructed with SOBI

fpzSOBI = reconstructed_SOBI(end,:);

%fpz length is about 35*128 samples. So we need 35 fpz segments. We will
%create a matrix where each row represents a segment.
fpz_segmentsSOBI=[];


for i=1:35
    
    begin=i*128-127;
    final=i*128;
    
    fpz_segmentsSOBI=[fpz_segmentsSOBI;fpzSOBI(begin:final)]; 
    
end

pfz_averagedSOBI=linspace(0,0,128);

for i=1:35 
     
    pfz_averagedSOBI=pfz_averagedSOBI+fpz_segmentsSOBI(i,:);
    
end

pfz_averagedSOBI=pfz_averagedSOBI/35;

figure('name','Synchronized Averaging - After reconstruction with SOBI')
for i=1:35
    plot(fpz_segmentsSOBI(i,:))
    hold on
end
plot(pfz_averagedSOBI,'r','LineWidth',3);
title('Synchronized Averaging - After reconstruction with SOBI')
xlabel('Samples')
ylabel('Amplitude')
axis tight

%% A plot with the 3 grand averages computed in a subplot

figure('name','Grand Averaging - Comparations')
subplot(3,1,1)
for i=1:35
    plot(fpz_segments(i,:))
    hold on
end
plot(pfz_averaged,'r','LineWidth',3);
axis tight
title('Grand Averaging - Original Signal')
xlabel('Samples')
ylabel('Amplitude')

subplot(3,1,2)
for i=1:35
    plot(fpz_segmentsICA(i,:))
    hold on
end
plot(pfz_averagedICA,'r','LineWidth',3);
axis tight
title('Grand Averaging - Reconstructed Signal with robust ICA')
xlabel('Samples')
ylabel('Amplitude')

subplot(3,1,3)
for i=1:35
    plot(fpz_segmentsSOBI(i,:))
    hold on
end
plot(pfz_averagedSOBI,'r','LineWidth',3);
axis tight
title('Grand Averaging - Reconstructed Signal with SOBI')
xlabel('Samples')
ylabel('Amplitude')






