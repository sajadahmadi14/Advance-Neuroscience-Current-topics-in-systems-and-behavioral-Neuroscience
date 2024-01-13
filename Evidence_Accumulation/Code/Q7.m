%% Question 7:
clear; clc; close all;

x0 = 0;
B = 0.001;
sigma = 1;
m = 0;
dt = 0.1;
thrU = 5;
thrD = -5;

NumTrl = 500;

figure;
for itrl = 1: NumTrl
    if mod(itrl, 100) == 0, disp(itrl); end
    [choice(itrl), RT(itrl), ~] = two_choice_trial(x0, B, dt, m, sigma, thrU, thrD);
    if choice(itrl) == 1, scatter(RT(itrl), choice(itrl), 'filled', 'r'); hold on; end
    if choice(itrl) == -1, scatter(RT(itrl), choice(itrl), 'filled', 'b'); hold on; end
end
ylabel('choice'); xlabel('Reaction Time'); axis tight


figure;
crctInd = find(choice == 1);
wrInd = find(choice == -1);
nbin = 50;
hist(RT,nbin);
ylabel('Histogram of all choices'); xlabel('RT'); title('all choices')

figure;
hist(RT(crctInd),nbin);
ylabel('Histogram of choice 1'); xlabel('RT'); title('choice 1')
figure;
hist(RT(wrInd),nbin);
ylabel('Histogram of choice -1'); xlabel('RT'); title('choice 2')


function [choice, RT, x] = two_choice_trial(x0, B, dt, m, sigma, thrU, thrD)
    x(1) = x0;
    i = 2;
    while 1
        x(i) = x(i-1) + B*dt + sigma * normrnd(m, sqrt(dt), 1);
        if x(i) > thrU
            RT = i * dt;
            choice = 1;
            x(i) = thrU;
            break;
        end
        if x(i) < thrD
            RT = i * dt;
            choice = -1;
            x(i) = thrD;
            break;
        end
        i = i + 1;
    end
end

