close all
clear all

degM = [14];

ndeg = length(degM);

fpavM = [3,5,7,9];

fpav = sprintf('Cheb%d_bet%d',degM(1),fpavM(1));
load(fpav)

nbet = length(fpavM);
nn = length(PredictTimeM);

vid = zeros(ndeg,nbet,nn);
vidlog = zeros(ndeg,nbet,nn);
minVer = zeros(ndeg,nbet,nn);
maxVer = zeros(ndeg,nbet,nn);
stdVer = zeros(ndeg,nbet,nn);
countN = zeros(ndeg,nbet,nn);
upper = zeros(ndeg,nbet,nn);
down = zeros(ndeg,nbet,nn);



for ideg = 1:length(degM)
    %for every beta
    for ii1 = 1:length(fpavM)
        fpav = sprintf('Cheb%d_bet%d',degM(ideg),fpavM(ii1));
        load(fpav)    
        vid(ideg,ii1,:) = median(ErrorM);
        vidlog(ideg,ii1,:) = median(log(ErrorM));
        minVer(ideg,ii1,:) = min(ErrorM);
        maxVer(ideg,ii1,:) = max(ErrorM);
        stdVer(ideg,ii1,:)= std(ErrorM);
        countN(ideg,ii1,:) = median(countNaNM);
        for ii2=1:length(PredictTimeM)
            upper(ideg,ii1,ii2) = median(ErrorM(ErrorM(:,ii2)>vid(ideg,ii1,ii2),ii2));
            down(ideg,ii1,ii2) = median(ErrorM(ErrorM(:,ii2)<vid(ideg,ii1,ii2),ii2));
        end
    end
end

clrs  =[0	0.447	0.741;
0.85	0.325	0.098;
0.929	0.694	0.125;
0.494	0.184	0.556;
0.466	0.674	0.188;
0.301	0.745	0.933;
0.635	0.078	0.184];

%legend
lg={};
figure(1), hold on
%for every degree
for ideg = 1:length(degM)
    %for every beta
    for ii1 = 1:length(fpavM)
        tmpvid = reshape(vid(ideg,ii1,:),[],1,1);
        tmpminver = reshape(minVer(ideg,ii1,:),[],1,1);
        tmpmaxver = reshape(maxVer(ideg,ii1,:),[],1,1);
        tmpstdver = reshape(stdVer(ideg,ii1,:),[],1,1);
        tmpup = reshape(upper(ideg,ii1,:),[],1,1);
        tmpdown = reshape(down(ideg,ii1,:),[],1,1);
        plot(PredictTimeM/CC.lyapTime, tmpvid,'Color',clrs(ii1,:))

        lg{ii1}=sprintf('$\\log(\\beta)=-%d$',fpavM(ii1));
    end
end


legend(lg,'Interpreter','latex')
xlabel('$\Lambda t_{pred}$','Interpreter','latex')
ylabel('RMSE')
set(gca,'yScale','log')




