clc;
clear all;
close all;

SCORE_GENERATION = 1;

%% Input Data Loading
for user = 1 : 100 
    imgno = 2;
    name = ['gallery_set\subject', int2str(user),'_img', int2str(imgno),'.pgm'] ;
    image = double(imread(name));
    A(:,( (2*user) - 1 ) ) = image(:);
    imgno = 3;
    name = ['gallery_set\subject', int2str(user),'_img', int2str(imgno),'.pgm'] ;
    image = double(imread(name));
    A(:,(2*user)) = image(:);
    name = ['probe_set\subject', int2str(user),'_img1.pgm'] ;
    image = double(imread(name));
    probe(:,user) = image(:);
end

%% Calculate Covariance Matrix and Eigen Vectors
meanface = mean(A.');
meanface = meanface';
figure;
imshow(reshape(uint8(meanface),50,50));

for gallery = 1 : (2 * user)
    A(:,gallery) = A(:,gallery) - meanface;
end
for p = 1 : user
    probe(:,p) = probe(:,p) - meanface;
end

At = A';
C = At * A;
% C = C / 200;
[V,eigenvalue] = eig(C);
eigenvalue = sort(diag(eigenvalue),'descend');
U = A * V;

%% Selecting Eigen Faces
R1RR = [];
for nc = 30 : 10 : 100
    eigenvector = zeros(2500,nc);
    j = 1;
    if (nc == 30)
        figure;
    end
    for i = 200 : -1 : (200 - nc + 1)
        if ( (j <= 10) && (nc == 30) )
            subplot(2,5,j);
            imshow(reshape(uint8(U(:,i) + meanface),50,50));
        end
        eigenvector(:,j) = U(:,i);
        j = j + 1;
    end

    %% Weight Calculation
    gallery_weights = zeros(nc,200);
    probe_weights = zeros(nc,100);
    gallery_weights = eigenvector' * A;
    probe_weights = eigenvector' * probe;

    %% Score Generation

    scores = zeros(100,200);
    for i = 1 : user
        for j = 1 : user
            if (SCORE_GENERATION == 1)
                scores(i,((2*j) - 1)) = sqrt(sum((probe_weights(:,i) - gallery_weights(:,(2*j) - 1 )).^2)) / 100000;
                scores(i,(2*j)) = sqrt(sum((probe_weights(:,i) - gallery_weights(:,(2*j))).^2)) / 100000;
            end
            if (SCORE_GENERATION == 2)
                scores(i,((2*j) - 1)) = sum(((probe_weights(:,i) - gallery_weights(:,(2*j) - 1 )).^2) ./ eigenvalue(1:nc,1)) / 100000;
                scores(i,(2*j)) = sum(((probe_weights(:,i) - gallery_weights(:,(2*j))).^2) ./ eigenvalue(1:nc,1)) / 100000;
            end
            if (SCORE_GENERATION == 3)
                scores(i,((2*j) - 1)) = sum( abs( probe_weights(:,i) - gallery_weights(:,(2*j) - 1 ) ) ) / 100000;
                scores(i,(2*j)) = sum( abs( probe_weights(:,i) - gallery_weights(:,(2*j)) ) ) / 100000;
            end
            if (SCORE_GENERATION == 4)
                scores(i,((2*j) - 1)) = (sum(abs( probe_weights(:,i) - gallery_weights(:,(2*j) - 1) ) .^3) ^ (1/3)) / 100000;
                scores(i,(2*j)) = (sum(abs( probe_weights(:,i) - gallery_weights(:,(2*j)) ) .^3) ^ (1/3)) / 100000;
            end
        end
    end
    scores = round(scores);

    %% Score Distributions
    S_genuine_scores = [];
    S_imposter_scores = [];
    for i = 1 : user
        S_genuine_scores(((2*i) - 1) : (2*i),1) = scores(i,((2*i) - 1) : (2*i));
        if (~(i == 1))
            S_imposter_scores = [S_imposter_scores scores(i,1 : ((2*i) - 2))];
        end
        S_imposter_scores = [S_imposter_scores scores(i,((2*i) + 1): 200)];
    end

    S_imposter_scores = S_imposter_scores';
    P_genuine_scores_interval = zeros( ((max(S_imposter_scores) + (100 - mod(max(S_imposter_scores),100))) / 5), 1);
    P_imposter_scores_interval = zeros( ((max(S_imposter_scores) + (100 - mod(max(S_imposter_scores),100))) / 5), 1);
    mean_genuine = 0;
    mean_imposter = 0;
    std_genuine = 0;
    std_imposter = 0;


    for i = 1:length(P_genuine_scores_interval)
         a = find(S_genuine_scores > (5 * (i -1)) & S_genuine_scores <= (5*i));
         P_genuine_scores_interval(i) = length(a)/length(S_genuine_scores);
         mean_genuine = mean_genuine + (P_genuine_scores_interval(i) * i);
    end

    for i = 1:length(P_imposter_scores_interval)
         a = find(S_imposter_scores > (5 * (i -1)) & S_imposter_scores <= (5*i));
         P_imposter_scores_interval(i) = length(a)/length(S_imposter_scores);
         mean_imposter = mean_imposter + (P_imposter_scores_interval(i) * i);
    end
    if (nc == 30)
        figure
        plot([1:5:(max(S_imposter_scores) + (100 - mod(max(S_imposter_scores),100)))],P_genuine_scores_interval,'-b',[1:5:(max(S_imposter_scores) + (100 - mod(max(S_imposter_scores),100)))],P_imposter_scores_interval,'-r');
        legend('Genuine','Imposter');
        xlabel('Distance Score');
        ylabel('Probability(score)');
        title('Score Distributions');
    end

    %% CMC Curve

    [S_rows,S_cols] = size(scores);
    user_scores = zeros(1,S_rows);
    user_rank = zeros(1,S_rows);
    RankT_percentage = zeros(1,S_rows);
    CMC = zeros(1,S_rows);
    temp1 = [];
    temp2 = [];

    for i = 1 : S_rows
        user_scores = scores(i,:);
        user_scores = sort(user_scores);
        temp1 = find(user_scores == scores(i,(2 *i) - 1));
        temp2 = find(user_scores == scores(i,(2 * i)));
        if (temp1(1) < temp2(1))
            user_rank(i) = temp1(1);
        else
            user_rank(i) = temp2(1);
        end
    end

    for i = 1:S_rows
        RankT_percentage(i) = length(find(user_rank == i))/S_rows;
    end
    
    R1RR = [R1RR RankT_percentage(1)];

    CMC = cumsum(RankT_percentage);
    CMC = CMC*100;

    if (nc == 30)
        figure;
        plot([1:S_rows],CMC);
        xlabel('Rank(t)');
        ylabel('Rank-t Identification Rate(%)');
        title('Cumulative Match Characteristic(CMC) Curve');
    end

    %% ROC Curve
    T = 1000;
    smin = min(scores(:));
    smax = max(scores(:));
    GAR = zeros(1,T);
    FAR = zeros(1,T);
    FRR = zeros(1,T);

    p = (smax - smin)/(T-1);
    for j = 1:T
        nj = smin + ((j - 1) * p);
        a = find(S_genuine_scores <= nj);
        GAR(j) = length(a)/length(S_genuine_scores);
        FRR(j) = 1 - GAR(j);

        b = find(S_imposter_scores <= nj);
        FAR(j) = length(b)/length(S_imposter_scores);
    end

    if (nc == 30)
        figure
        plot(FAR * 100,GAR * 100);
        xlabel('False Accept Rate (%)');
        ylabel('Genuine Accept Rate (%)');
        title('Receiver Operating Charactristic (ROC) Curve');
    end
end

%% Plot of Recognition Rate vs No of coefficients
for i = length(R1RR) : -1 : 2
     if (R1RR(i) ~= R1RR(i - 1))
         break;
     end
end
R1RR(2,:) = [30:10:100];
figure
plot(R1RR(2,:),R1RR(1,:) * 100);
xlabel('No of Coefficients');
ylabel('Rank - 1 Recognition Rate');
title('Rank - 1 Recognition Rate vs No of coefficients');
        
h = msgbox(sprintf('No of coefficients = %d',(i + 2) * 10),'Recognition Rate Constant');
set(h, 'position', [100 300 200 50]);