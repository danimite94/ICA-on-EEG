% A simple demonstration illustrating the performance of RobustICA on synthetic mixtures.
%
%
% Please, report any bugs, comments or suggestions to <a href = "mailto:zarzoso@i3s.unice.fr">zarzoso(a)i3s.unice.fr</a>.
%
%
% HISTORY:
% 
%    <Please add modification date here>: - <please add modification details here>
%
%    2008/03/29: - added a bit of interactivity with the user
%                  (pause, "press a key to continue", etc.)
%
%    2008/03/27: - created by Vicente Zarzoso.


clc;

disp(' ')
disp('=== ROBUSTICA DEMO =====================================================')

%%% Generate different types of synthetic sources

disp(' ');
disp('>>> Generating random source signals');

T = 1e3;            % sample size
K = 20;             % number of sources and observations

disp(['- sample size: ', num2str(T)]);

j = sqrt(-1);

S_true = zeros(K, T);

disp(['> sources #1-', num2str(K/4), ': real-valued, sub-Gaussian']);
S_true(1:K/4,:) = (rand(K/4, T) > 0.5);     % a few sub-Gaussian real-valued sources

disp(['> sources #', num2str(K/4+1), '-', num2str(K/2), ': real-valued, super-Gaussian']);
S_true(K/4+1:K/2,:) = (rand(K/4, T) > 0.9);     % a few super-Gaussian real-valued sources 

disp(['> sources #', num2str(K/2+1), '-', num2str(3*K/4), ': complex-valued, sub-Gaussian']);
S_true(K/2+1:3*K/4,:) = (rand(K/4, T) > 0.5) + j*(rand(K/4, T) > 0.5);   % a few sub-Gaussian complex-valued sources

disp(['> sources #', num2str(3*K/4+1), '-', num2str(K), ': complex-valued, super-Gaussian']);
S_true(3*K/4+1:K,:) = (rand(K/4, T) > 0.9) + j*(rand(K/4, T) >= 0.9);    % a few super-Gaussian complex-valued sources

S_true = S_true - mean(S_true')'*ones(1, T);            % remove mean
S_true = diag(T./sqrt(diag(S_true*S_true')))*S_true;    % unit-variance normalization



disp(' '); disp('Please, press any key to continue...'); disp(' ');
pause



disp(' ');
disp('>>> Plotting original sources (real part only)');

Tplot = min(T, 100);  % plot just a few samples, for clarity

figure
set(gcf, 'Name', 'RobustICA demo: Original source signals');
subplot(K, 1, 1);
for k = 1:K
    subplot(K, 1, k);
    plot(0:Tplot-1, real(S_true(k, 1:Tplot)), 'b');   % plot real parts only
    set(gca, 'XTickLabel', []);
    ylabel(['s_{', num2str(k), '}']);
end
set(gca, 'XTickLabelMode', 'auto');
xlabel('sample number');
subplot(K, 1, 1);
title('Original source signals');



disp(' '); disp('Please, press any key to continue...'); disp(' ');
pause




disp(' ');
disp('>>> Normalized kurtosis of original sources:')

kurt = zeros(1, K);

for k = 1:K
    s = S_true(k, :);
    kurt(k) = (mean(abs(s).^4) - 2*mean(abs(s).^2)^2 - abs(mean(s.^2))^2)/mean(abs(s.^2))^2;
    disp(['- source #', num2str(k), ': ', num2str(kurt(k))]);
end



disp(' '); disp('Please, press any key to continue...'); disp(' ');
pause



%%% Generate mixing matix and observed signals

disp(' ');
disp('>>> Generating complex-valued random mixture');

H_true = randn(K, K) + j*randn(K, K);       % random complex-valued mixing matrix

X = H_true*S_true;      % observed mixture


% plot true souces and estimated sources

disp(' ');
disp('>>> Plotting observed mixture (real part only)');

figure
set(gcf, 'Name', 'RobustICA demo: Observed mixtures');
subplot(K, 1, 1);
for k = 1:K
    subplot(K, 1, k);
    plot(0:Tplot-1, real(X(k, 1:Tplot)), 'b');    % plot real parts only
    set(gca, 'XTickLabel', []);
    ylabel(['x_{', num2str(k), '}']);
end
set(gca, 'XTickLabelMode', 'auto');
xlabel('sample number');
subplot(K, 1, 1);
title('Observed mixtures');



disp(' '); disp('Please, press any key to continue...'); disp(' ');
pause



%%% Perform source separation

tol = 1e-3;     % termination threshold parameter
max_it = 1e3;   % maximum number of iterations per independent component

% please, comment out one of the following lines to test a different implementation of RobustICA
% (other modes of operation are possible; please type 'doc robustica' or 'help robustica' for details)

[S, H, iter, W] = robustica(X, [], tol, max_it, 1, 'r', 0, [], 1);   % regression-based deflation

%[S, H, iter, W] = robustica(X, [], tol, max_it, 1, 'o', 0, [], 1);   % deflationary orthogonalization

%[S, H, iter, W] = robustica(X, [1, 1, zeros(1, K-2)], -1, 10, 1, 'r', 0, [], 1);  % extract two super-Gaussian sources first, running a fixed number of
                                                                                  % iterations per independent component

%[S, H, iter, W]= robustica(X, [-1, -1, zeros(1, K-2)], tol, max_it, 1, 'r', 0, [], 1); % extract two sub-Gaussian sources first



disp('Please, press any key to continue...'); disp(' ');
pause



%%% Measure source separation performance

% test type of estimated sources (to see, e.g., if extraction ordering is as required)

disp(' ');
disp('>>> Normalized kurtosis of extracted sources:')

kurt = zeros(1, K);

for k = 1:K
    s = S(k, :);
    kurt(k) = (mean(abs(s).^4) - 2*mean(abs(s).^2)^2 - abs(mean(s.^2))^2)/mean(abs(s.^2))^2;
    disp(['- source #', num2str(k), ': ', num2str(kurt(k))]);
end



disp(' '); disp('Please, press any key to continue...'); disp(' ');
pause
    



disp(' ');
disp('>>> Sorting and scaling the estimated sources');

% re-arrange the sources and work out reconstruction error
% (the following ordering method is suboptimal, but works fine of sources are well separated)

P = zeros(K, K);        % permutation matrix
for k = 1:K
    s = S(k, :);
    r = S_true*s';      % cross-correlation
    [r_max, pos_max] = max(abs(r));
    P(pos_max, k) = r(pos_max)/(s*s');  % ordering with optimum scale
end

S = P*S;

% plot true souces and estimated sources

disp(' ');
disp('>>> Plotting actual and estimated sources (real part only)');

figure
set(gcf, 'Name', 'RobustICA demo: Source separation results');
subplot(K, 1, 1);
for k = 1:K
    subplot(K, 1, k);
    plot(0:Tplot-1, real(S_true(k, 1:Tplot)), 'b'); hold on;    % plot real parts only
    plot(0:Tplot-1, real(S(k, 1:Tplot)), 'r');
    set(gca, 'XTickLabel', []);
    ylabel(['s_{', num2str(k), '}']);
end
set(gca, 'XTickLabelMode', 'auto');
xlabel('sample number');
subplot(K, 1, 1);
title('RobustICA''s source separation results; blue lines: actual sources; red lines: estimated sources');



disp(' '); disp('Please, press any key to continue...'); disp(' ');
pause



disp(' ');
disp('>>> Normalized mean square error of estimated sources:')

SMSE = sum(sum(abs((S_true - S).^2)))/sum(sum(abs(S_true).^2));

disp(['SMSE = ', num2str(100*SMSE), '% = ', num2str(10*log10(SMSE)), ' dB']);

disp(' '); 

disp('========================================================================')
         
disp(' '); 




