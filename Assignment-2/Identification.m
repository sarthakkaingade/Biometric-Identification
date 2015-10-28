clc;
clear all;      %Clear Workspace
close all;      %Close all windows

[f,p] = uigetfile;                  %Accept Similarity Matrix from user
S = load(cat(2,p,f));
[S_rows,S_cols] = size(S);          %Size of matrix

%% Score Distributions
S_genuine_scores = diag(S);                             %Separate Genuine Scores
S_imposter_scores = zeros((S_rows * S_cols)- S_rows,1);
P_genuine_scores_interval = zeros(200,1);
P_imposter_scores_interval = zeros(200,1);
mean_genuine = 0;
mean_imposter = 0;
std_genuine = 0;
std_imposter = 0;
k = 1;

for i = 1:S_rows
    for j = 1:S_cols
        if (~(i == j)) 
        S_imposter_scores(k) = S(i,j);                   %Separate Imposter Scores
        k = k+1;
        end
    end
end

for i = 1:length(P_genuine_scores_interval)
     a = find(S_genuine_scores > (1 * (i -1)) & S_genuine_scores <= (1*i));
     P_genuine_scores_interval(i) = length(a)/length(S_genuine_scores);          %Calculate Probability of Genuine Scores
     mean_genuine = mean_genuine + (P_genuine_scores_interval(i) * i);
end

for i = 1:length(P_imposter_scores_interval)
     a = find(S_imposter_scores > (1 * (i -1)) & S_imposter_scores <= (1*i));
     P_imposter_scores_interval(i) = length(a)/length(S_imposter_scores);        %Calculate Probability of Imposter Scores
     mean_imposter = mean_imposter + (P_imposter_scores_interval(i) * i);
end

figure
plot([1:1:200],P_genuine_scores_interval,'-b',[1:1:200],P_imposter_scores_interval,'-r');
legend('Genuine','Imposter');
xlabel('Distance Score');
ylabel('Probability(score)');
title('Score Distributions');

%% CMC Curve

user_scores = zeros(1,S_rows);
user_rank = zeros(1,S_rows);
RankT_percentage = zeros(1,S_rows);
CMC = zeros(1,S_rows);
temp = [];

for i = 1 : S_rows
    user_scores = S(i,:);
    user_scores = sort(user_scores);
    temp = find(user_scores == S(i,i));
    user_rank(i) = temp(1);                              %Calculate User Rank
end

for i = 1:S_rows
    RankT_percentage(i) = length(find(user_rank == i))/S_rows;
end

CMC = cumsum(RankT_percentage);
CMC = CMC*100;

figure;
plot([1:S_rows],CMC);
xlabel('Rank(t)');
ylabel('Rank-t Identification Rate(%)');
title('Cumulative Match Characteristic(CMC) Curve');

%% d-prime value

for i = 1:length(P_genuine_scores_interval)
    std_genuine = std_genuine + (P_genuine_scores_interval(i) * ((i - mean_genuine)^2));
end

for i = 1:length(P_imposter_scores_interval)
    std_imposter = std_imposter + (P_imposter_scores_interval(i) * ((i - mean_imposter)^2));
end

mean_genuine = mean(S_genuine_scores);
mean_imposter = mean(S_imposter_scores);
std_genuine = std(S_genuine_scores);
std_imposter = std(S_imposter_scores);

Dindex = (sqrt(2) * abs(mean_genuine - mean_imposter))/(sqrt(std_genuine^2 + std_imposter^2));

msgbox(sprintf('d'' = %f',Dindex),'d-prime value');

%% Lowest Rank at recognition performance of atleast 85%

Sum = 0;

for i = 1 : length(RankT_percentage)
    Sum = Sum + RankT_percentage(i);
    if Sum >= 0.85
        break;
    end
end

h = msgbox(sprintf('Rank(atleast 85%%) = %d',i),'Lowest Rank');
set(h, 'position', [100 300 200 50]);
