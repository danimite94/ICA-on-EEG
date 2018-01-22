% session_7_1_1.m
% October 2015

clear; close all; clc; 
addpath ./robustica_package/

%% Data
t=(0:9999 )*0.001; 
Sources=zeros(4, size(t,2)); 

Sources(1,:)= 2*sin(2*pi*7*t); 
Sources(2,:)= 0.1*randn(1,size(t,2));
Sources(3,:)= sawtooth(2*pi*3*t) ;

d=[0:1:10; 0.8.^ (0:10)]';
Sources(4,:) =2* pulstran(t,d,'gauspuls',10,0.3);
reference=Sources(3,:); 

% normalize the sources
Sources = zscore(Sources,0,2); 

% A == Mixing matrix generating 6 mixtures
A=randn(6,4); 

X=A*Sources; 

%display
figure(1); 
for i=1:size(Sources,1)
   plot(Sources(i,:)); 
   title(['Source s(t) #', num2str(i)]);
   pause; 
   clf;  
    
end


figure(1); 

for i=1:size(X,1)
   plot(X(i,:));
   title(['Mixture x(t) #', num2str(i)]);
   pause; 
   clf;  
    
end

%% Algorithms

% ROBUST ICA
tol = 1e-3;     % termination threshold parameter
max_it = 1000;   % maximum number of iterations per independent component
[y_robustICA, A_robustICA] = robustica(X, [], tol, max_it, 1, 'r', 0, [], 0);   % regression-based deflation

% SOBI
[A_SOBI,y_SOBI] = acsobiro(X,size(X,1),100);

y_SOBI_unormalized=y_SOBI;

% JADE
H_JADE=jader(X,size(X,1)); 
y_JADE=H_JADE*X; 
y_JADE=real(y_JADE); 

% CCA
[y_CCA,H_CCA,r] = ccaqr(X);
y_CCA = real(y_CCA);

        % comment: some algorithms are developed to process real and complex
        % signals, therefore use of the function "real"
        
        
 %% Display
 disp('robustICA'); 
for i=1:size(y_robustICA,1)
   plot(y_robustICA(i,:)); 
   title(['Source estimate robustICA #', num2str(i)]);
   pause; 
   clf; 
   
end


disp('SOBI')
for i=1:size(y_SOBI,1)
   plot(y_SOBI(i,:)); 
    title(['Source estimate SOBI #', num2str(i)]);
   pause; 
   clf; 
       
end

disp('JADE')
for i=1:size(y_JADE,1)
   plot(y_JADE(i,:)); 
   title(['Source estimate JADE #', num2str(i)]);
   pause; 
   clf; 
       
end

disp('CCA')
for i=1:size(y_CCA,1)
   plot(y_CCA(i,:)); 
   title(['Source estimate CCA #', num2str(i)]);
   pause; 
   clf; 
       
end

%% Exercise

% write the CODE here...

% 1)Normalize the y-estimates: set mean value to 0 and variance to 1. Do you know why?

for i=1:6
    y_robustICA(i,:)=zscore(y_robustICA(i,:));
    y_CCA(i,:)=zscore(y_CCA(i,:));
    y_JADE(i,:)=zscore(y_JADE(i,:));
    y_SOBI(i,:)=zscore(y_SOBI(i,:));
end


% 2) Automatically match each Source with (one of) its best estimate(s) (Source(1) with y_SOBI(#,:),y_robustICA(#,:))
%  Use squared or absolute correlation as matching criterion

% Robustica

Best_Estimates_Robustica=[];
positive_negative_corr_Robustica=[];

%we calculated the correlation for each sources with all the estimates.
%Since it is a correlation between 2 vectors, we will get a vector. So, in order
%to know which one was the best matched, we summed each squared correlation
%vector in order to have a comparable value between all correlations calculated.

for i=1:4
    source_evaluated=Sources(i,:);
    Estimations_correlations=[];
    Estimations_correlations_squared_value=[];
    for j=1:6
        correlation=xcorr(source_evaluated,y_robustICA(j,:));
        Estimations_correlations=[Estimations_correlations; correlation];
        Estimations_correlations_squared_value=[Estimations_correlations_squared_value; sum((correlation.^2)) ];
    end
    
    index=find(Estimations_correlations_squared_value==max(Estimations_correlations_squared_value));
    
    Best_Estimates_Robustica=[Best_Estimates_Robustica; y_robustICA(index,:) ];
    
    if i==1
        index_source1=index;   
    end

    positive=sum(Estimations_correlations(index))>=1;
    positive_negative_corr_Robustica=[positive_negative_corr_Robustica; positive];
    
end

% we have verified that sometimes the first and sources will be match paired with the same estimated
% component. And we have verified as well (by plotting), that the best pair
% corresponds to the first source. so we need to find a different estimate
% to the source 2 by eliminating the hypothesis of the component already
% chosen. We don't need to generalize this method and this verification to
% all the components since we have verified that will only occur on that
% situation for the given situation. on that way, we will save some
% computation

if Best_Estimates_Robustica(2,:)==Best_Estimates_Robustica(1,:)
    new_mixtures=y_robustICA;
    new_mixtures(index_source1,:)=[];
    
    source_evaluated=Sources(2,:);
    Estimations_correlations=[];
    Estimations_correlations_squared_value=[];
    
    for j=1:5
        correlation=xcorr(source_evaluated,new_mixtures(j,:));
        Estimations_correlations=[Estimations_correlations; correlation];
        Estimations_correlations_squared_value=[Estimations_correlations_squared_value; sum((correlation.^2)) ];
    end
    
    index=find(Estimations_correlations_squared_value==max(Estimations_correlations_squared_value));
    
    Best_Estimates_Robustica(2,:)=new_mixtures(index,:);
    
end

% the plot of each source with the best estimate found

% Display
 disp('robustICA estimates with Sources'); 
for i=1:4
    
   subplot(2,1,1)
   plot(Best_Estimates_Robustica(i,:)); 
   title(['Source estimate robustICA #', num2str(i)]);
   
   subplot(2,1,2)
   plot(Sources(i,:));
   title(['Original Source #', num2str(i)]);
   pause; 
   clf; 
   
end



%% CCA

Best_Estimates_CCA=[];
positive_negative_corr_CCA=[];

for i=1:4
    source_evaluated=Sources(i,:);
    Estimations_correlations=[];
    Estimations_correlations_squared_value=[];
    for j=1:6
        correlation=xcorr(source_evaluated,y_CCA(j,:));
        Estimations_correlations=[Estimations_correlations; correlation];
        Estimations_correlations_squared_value=[Estimations_correlations_squared_value; sum((correlation.^2)) ];
    end
    
    index=find(Estimations_correlations_squared_value==max(Estimations_correlations_squared_value));
    
    if i==1
        index_source1=index;   
    end
    
    Best_Estimates_CCA=[Best_Estimates_CCA; y_CCA(index,:) ];

    positive=sum(Estimations_correlations(index))>=1;
    positive_negative_corr_CCA=[positive_negative_corr_CCA; positive];
    
end

% correcting the 2nd source estimate if needed

if Best_Estimates_CCA(2,:)==Best_Estimates_CCA(1,:)
    new_mixtures=y_CCA;
    new_mixtures(index_source1,:)=[];
    
    source_evaluated=Sources(2,:);
    Estimations_correlations=[];
    Estimations_correlations_squared_value=[];
    
    for j=1:5
        correlation=xcorr(source_evaluated,new_mixtures(j,:));
        Estimations_correlations=[Estimations_correlations; correlation];
        Estimations_correlations_squared_value=[Estimations_correlations_squared_value; sum((correlation.^2)) ];
    end
    
    index=find(Estimations_correlations_squared_value==max(Estimations_correlations_squared_value));
    
    Best_Estimates_CCA(2,:)=new_mixtures(index,:);
    
end




% Display
 disp('CCA estimates with Sources'); 
for i=1:4
    
   subplot(2,1,1)
   plot(Best_Estimates_CCA(i,:)); 
   title(['Source estimate CCA #', num2str(i)]);
   
   subplot(2,1,2)
   plot(Sources(i,:));
   title(['Original Source #', num2str(i)]);
   pause; 
   clf; 
   
end

%% SOBI

Best_Estimates_SOBI=[];
positive_negative_corr_SOBI=[];
index_sawtooth=0;

for i=1:4
    source_evaluated=Sources(i,:);
    Estimations_correlations=[];
    Estimations_correlations_squared_value=[];
    for j=1:6
        correlation=xcorr(source_evaluated,y_SOBI(j,:));
        Estimations_correlations=[Estimations_correlations; correlation];
        Estimations_correlations_squared_value=[Estimations_correlations_squared_value; sum((correlation.^2)) ];
    end
    
    index=find(Estimations_correlations_squared_value==max(Estimations_correlations_squared_value));
    
    if i==3
        index_sawtooth=index;
    end
    
    if i==1
        index_source1=index;
    end
    Best_Estimates_SOBI=[Best_Estimates_SOBI; y_SOBI(index,:) ];

    positive=sum(Estimations_correlations(index))>=1;
    positive_negative_corr_SOBI=[positive_negative_corr_SOBI; positive];
    
end

%%Correcting the 2nd source estimate if needed

if Best_Estimates_SOBI(2,:)==Best_Estimates_SOBI(1,:)
    new_mixtures=y_SOBI;
    new_mixtures(index_source1,:)=[];
    
    source_evaluated=Sources(2,:);
    Estimations_correlations=[];
    Estimations_correlations_squared_value=[];
    
    for j=1:5
        correlation=xcorr(source_evaluated,new_mixtures(j,:));
        Estimations_correlations=[Estimations_correlations; correlation];
        Estimations_correlations_squared_value=[Estimations_correlations_squared_value; sum((correlation.^2)) ];
    end
    
    index=find(Estimations_correlations_squared_value==max(Estimations_correlations_squared_value));
    
    Best_Estimates_SOBI(2,:)=new_mixtures(index,:);
    
end

% Display
 disp('SOBI estimates with Sources'); 
for i=1:4
    
   subplot(2,1,1)
   plot(Best_Estimates_SOBI(i,:)); 
   title(['Source estimate SOBI #', num2str(i)]);
   
   subplot(2,1,2)
   plot(Sources(i,:));
   title(['Original Source #', num2str(i)]);
   pause; 
   clf; 
   
end


%% JADE

Best_Estimates_JADE=[];
positive_negative_corr_JADE=[];

for i=1:4
    source_evaluated=Sources(i,:);
    Estimations_correlations=[];
    Estimations_correlations_squared_value=[];
    for j=1:6
        correlation=xcorr(source_evaluated,y_JADE(j,:));
        Estimations_correlations=[Estimations_correlations; correlation];
        Estimations_correlations_squared_value=[Estimations_correlations_squared_value; sum((correlation.^2)) ];
    end
    
    index=find(Estimations_correlations_squared_value==max(Estimations_correlations_squared_value));
    
    Best_Estimates_JADE=[Best_Estimates_JADE; y_JADE(index,:) ];
    
    if i==1
        index_source1=index;
    end
    
    positive=sum(Estimations_correlations(index))>=1;
    positive_negative_corr_JADE=[positive_negative_corr_JADE; positive];
    
    
end

%correcting 2nd source estimate if needed

if Best_Estimates_JADE(2,:)==Best_Estimates_JADE(1,:)
    new_mixtures=y_JADE;
    new_mixtures(index_source1,:)=[];
    
    source_evaluated=Sources(2,:);
    Estimations_correlations=[];
    Estimations_correlations_squared_value=[];
    
    for j=1:5
        correlation=xcorr(source_evaluated,new_mixtures(j,:));
        Estimations_correlations=[Estimations_correlations; correlation];
        Estimations_correlations_squared_value=[Estimations_correlations_squared_value; sum((correlation.^2)) ];
    end
    
    index=find(Estimations_correlations_squared_value==max(Estimations_correlations_squared_value));
    
    Best_Estimates_JADE(2,:)=new_mixtures(index,:);
    
end

% Display
 disp('JADE estimates with Sources'); 
for i=1:4
    
   subplot(2,1,1)
   plot(Best_Estimates_JADE(i,:)); 
   title(['Source estimate JADE #', num2str(i)]);
   
   subplot(2,1,2)
   plot(Sources(i,:));
   title(['Original Source #', num2str(i)]);
   pause; 
   clf; 
   
end
 
% ...
 
% 3) calculate the Root Mean Squared Error (RMSE) between matched Source and estimate
%  y(t). Pay attention that the signals are zero mean,
%  standardized and that they sould have the same sign. That is, use lower
%  RMSE value (e.g. RMSE (Source, y_SOBI ) and RMSE (Source, - y_SOBI))

RMSE_robustICA=[];
RMSE_CCA=[];
RMSE_SOBI=[];
RMSE_JADE=[];

% RobustICA method

for i=1:4
    source=Sources(i,:);
    estimate=Best_Estimates_Robustica(i,:);
    inversed_estimate=fliplr(estimate);
    
    RMES1=sqrt(sum((source-estimate).^(2))/length(source));
    RMES2=sqrt(sum((source-inversed_estimate).^(2))/length(source));
   
    RMSE_val=min(RMES1,RMES2);
    
    RMSE_robustICA=[RMSE_robustICA RMSE_val];
end

% CCA method

for i=1:4
    source=Sources(i,:);
    estimate=Best_Estimates_CCA(i,:);
    inversed_estimate=fliplr(estimate);
    
    RMES1=sqrt(sum((source(1:end-1)-estimate).^(2))/length(source(1:end-1)));
    RMES2=sqrt(sum((source(1:end-1)-inversed_estimate).^(2))/length(source(1:end-1)));
   
    RMSE_val=min(RMES1,RMES2);
    
    RMSE_CCA=[RMSE_CCA RMSE_val];
end

% SOBI method

for i=1:4
    
    source=Sources(i,:);
    estimate=Best_Estimates_SOBI(i,:);
    inversed_estimate=fliplr(estimate);
    
    RMES1=sqrt(sum((source-estimate).^(2))/length(source));
    RMES2=sqrt(sum((source-inversed_estimate).^(2))/length(source));
   
    RMSE_val=min(RMES1,RMES2);
    
    RMSE_SOBI=[RMSE_SOBI RMSE_val];
    
end

% JADE method

for i=1:4
    
    source=Sources(i,:);
    estimate=Best_Estimates_JADE(i,:);
    inversed_estimate=fliplr(estimate);
    
    RMES1=sqrt(sum((source-estimate).^(2))/length(source));
    RMES2=sqrt(sum((source-inversed_estimate).^(2))/length(source));
   
    RMSE_val=min(RMES1,RMES2);
    
    RMSE_JADE=[RMSE_JADE RMSE_val];
    
end


% 4) Repeat the RMSE calculations 20 times for different A - automatically 
% for each Source and each algorithm. Find the mean values.
% Use the code of the previous parts to create a separate script RMSE.m for this

%in here, by using the function RMSE, we will obtain the mean of the RMSE
%value for each source for each algorithm. We will obtain a matrix 4X4 in
%which each row represents a algorithm: Row 1: Robust ICA, Row 2: CCA, Row
%3: SOBI, Row 4: JADE. Each column will represent a source. Column 1:
%Source 1, Column 2: Source 2, Column 3 : Source 3, Column 4: Source 4.

meansRMSE=RMSE(Sources,20);

% 5) Use SOBI to remove the "sawtooth" signal and reconstruct the mixtures signal. 
% The sawtooth signal is saved in the variable reference
% so USE this variable 'reference' to find the sawtooth component and
% change A_SOBI in order to remove it from the mixtures
% hint: x=A_SOBI*y_SOBI (y_SOBI should be unnormalized here, why?)

clear Sources  

% we need to know which row of y_SOBI_unormalized corresponds to the
% sawtooth. so, we can use our vector Best_Estimates_SOBI to know which
% signal corresponds to the sawtooth. Sawtooth is our source number 3! (the
% 3rd row)
%index_sawtooth will provide uns that information: which row of y_Sobi
%belongs to the sawtooth. If we know that index, we can put the respective
%column on A_SOBI with 0 values, in order to reconstruct the signal without
%the sawtooth estimate. 

A_SOBI(:,index_sawtooth)=0;

xx=A_SOBI*y_SOBI_unormalized; %reconstruction 

%make subplots of the original mixtures next to the mixtures from which the
%sawtooth has been removed

%%

figure('name','Original Mixtures and Reconstructed Mixtures without the sawtooth')

for i=1:6
    subplot(6,2,2*i-1)
    plot(X(i,3000:4000))
    title(['Original Mixture #', num2str(i)]);
    axis tight
    
    subplot(6,2,2*i)
    plot(xx(i,3000:4000))
    title(['Mixture without sawtooth #', num2str(i)]);
    axis tight
end
