currFolders = ["All_1uM", "All_10nM"];

numOfExp = numel(currFolders);

s = zeros(864,numOfExp);
lcc = zeros(864,numOfExp);
minLengthT = 5;
usedLengthT = 2;
tx0 = zeros(minLengthT,numOfExp);
tx1 = zeros(minLengthT,numOfExp);

for i=1:numOfExp
    if contains(currFolders(i),'_1uM')
        s(:,i) = 1000;
    else
        s(:,i) = 10;
    end
    
    tzones = [1,72; 73,264; 265,336; 337,528; 529,792];

    tx0(1:length(tzones(:,1)),i) = tzones(1:end,1)';
    tx1(1:length(tzones(:,2)),i) = tzones(1:end,2)';
    if length(tzones(:,1)) < minLengthT
        minLengthT = length(tzones(:,2));
    end
end

tx0 = tx0(1:usedLengthT,:);
tx1 = tx1(1:usedLengthT,:);
tx1(end,:) = min(tx1(end,:));
lastFrame = tx1(end,1);
tmax = lastFrame/12;

fr = (1:lastFrame)';

s_set = cat(1,ones(tx1(1),1),zeros(tx1(2)-tx1(1),1));

s(1:length(s_set),:) = s_set.*s(1:length(s_set),:);

y = zeros(tx1(end,1),numOfExp);

%% get the experiments' data
xdata = [];

for i=1:numOfExp
    
    IntAvg = readmatrix(fullfile(currFolders(i),strcat(currFolders(i),"_s5_Int_Avg.csv")));
    DivTime = readmatrix(fullfile(currFolders(i),strcat(currFolders(i),"_avgDivTime_Roll50.csv")));
        
    lcc(1:lastFrame,i) = log(2)./DivTime(1:lastFrame,2);
    
    ts = IntAvg(fr,1);
    y(:,i) = IntAvg(fr,2);
    xdata = cat(1,xdata,ts+tmax*(i-1));
end

%% set the coefficients' boundaries
llim = 1;
ulim = 20;
startVal = 10;

xmin = [0,   4, 0.1, 1/16, 1000, 0.1,  20,    0,    1, 0.1,    1,    0]; % lower limits
xmax = [4,  40, 2.5,  0.7, 5000, 2.5, 500,    1,  100,  10,  500, 1e-2]; % upper limits

x0 = load(fullfile('loop','memory','1uMCoeffs_v2.mat'), 'meanCoeffVals1uM').meanCoeffVals1uM;

x0 = cat(2,x0(1:length(xmin)),startVal,1); % take only the first 12 params

xmin = cat(2,xmin,llim,1); % set st and fix rt
xmax = cat(2,xmax,ulim,1);

ydata = y(:);
    
%% plot experiment data
figure
%%
hold on
for i = y
    plot(ts,i)
end

%% get the previus fit's ID number
numberinfile = struct2cell(dir(fullfile('loop','memory')));
numberinfile = numberinfile(1,:);
numberinfile = regexp(numberinfile(contains(numberinfile,'multiple') & ~contains(numberinfile,'lsq') & ~contains(numberinfile,'rmse')),'\d*','Match');

maxnumber = 0;

for i = 1:length(numberinfile)
    if str2double(numberinfile{i}(1)) > maxnumber
        maxnumber = str2double(numberinfile{i}(1));
    end
end

maxnumber = maxnumber + 1;

%% fit
fitfun = fittype( @(c,l1,ld,lg,k2n,k2p,k3p,ksn,Ks,kg,a,A,st,rt,x)...
    multiple_exp_fit(x,tmax,lcc,s,c,l1,ld,lg,k2n,k2p,k3p,ksn,Ks,kg,a,A,st,rt,numOfExp));

[fitted_curve,gof] = fit(xdata,ydata,fitfun,'StartPoint',x0,'Lower',xmin,'Upper',xmax); % fit
coeffvals = coeffvalues(fitted_curve);
yfit = fitted_curve(xdata);
    
save(fullfile('loop','memory',strcat('multiple',string(maxnumber),'.mat')));
for j = 1:numOfExp
    plot(ts,yfit(fr+fr(end)*(j-1)),'DisplayName',strcat('multiple',string(maxnumber),'.mat',' - ',currFolders(j)))
end
