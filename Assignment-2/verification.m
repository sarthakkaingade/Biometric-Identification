clc;
clear all;  %Clear Workspace
close all;  %Close all windows

[f,p] = uigetfile;      %Accept Similarity Matrix from user
S = load(cat(2,p,f));
[S_rows,S_cols] = size(S);      %Size of matrix

%% ROC Curve
S_genuine_scores = diag(S);                                 %Separate Genuine Scores
S_imposter_scores = zeros((S_rows * S_cols)- S_rows,1);
T = 1000;
smin = min(S(:));
smax = max(S(:));
GAR = zeros(1,T);
FAR = zeros(1,T);
FRR = zeros(1,T);

k = 1;
for i = 1:S_rows
    for j = 1:S_cols
        if (~(i == j)) 
        S_imposter_scores(k) = S(i,j);                      %Separate Imposter Scores
        k = k+1;
        end
    end
end

p = (smax - smin)/(T-1);
for j = 1:T
    nj = smin + ((j - 1) * p);
    a = find(S_genuine_scores <= nj);
    GAR(j) = length(a)/length(S_genuine_scores);            %Calculate Genuine Accept Rate
    FRR(j) = 1 - GAR(j);                                    %Calculate Flase Reject Rate
    
    b = find(S_imposter_scores <= nj);
    FAR(j) = length(b)/length(S_imposter_scores);           %Calculate False Accept Rate
end

figure
plot(FAR * 100,GAR * 100);
xlabel('False Accept Rate (%)');
ylabel('Genuine Accept Rate (%)');
title('Receiver Operating Charactristic (ROC) Curve');

%% Equal Error Rate

[a,b] = min(abs(FRR - FAR));
op = smin + ((b - 1) * p);                              %Operating Point
EER = FRR(b);
msgbox(sprintf('EER(at %f) = %f',op,EER * 100),'Equal Error Rate');

%% FRR computation at various FAR

a = find(FAR >= 0.01);
FRR1 = FRR(a(1) - 1);
h = msgbox(sprintf('FRR(%%) when FAR is 1%% = %f',FRR1 * 100),'FRR');
set(h, 'position', [100 400 200 50]);

a = find(FAR >= 0.05);
FRR1 = FRR(a(1) - 1);
h = msgbox(sprintf('FRR(%%) when FAR is 5%% = %f',FRR1 * 100),'FRR');
set(h, 'position', [100 300 200 50]);

a = find(FAR >= 0.1);
FRR1 = FRR(a(1) - 1);
h = msgbox(sprintf('FRR(%%) when FAR is 10%% = %f',FRR1 * 100),'FRR');
set(h, 'position', [100 200 200 50]);

a = find(FAR >= 0.2);
FRR1 = FRR(a(1) - 1);
h = msgbox(sprintf('FRR(%%) when FAR is 20%% = %f',FRR1 * 100),'FRR');
set(h, 'position', [100 100 200 50]);
