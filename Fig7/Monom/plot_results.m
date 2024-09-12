close all
clear all

degM = [21];

ndeg = length(degM);

fpavM = [9];

reiksm = {'1proc.'};

fpav = sprintf('Monom%d_bet%d',degM(1),fpavM(1));
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

    
%for every degree


for ideg = 1:length(degM)
    %for every beta
    for ii1 = 1:length(fpavM)
        fpav = sprintf('Monom%d_bet%d',degM(ideg),fpavM(ii1));
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

%         for i3=1:size(ErrorM,1)
%             figure(1), hold on
%             scatter(PredictTimeM/CC.lyapTime,(ErrorM(i3,:)),'o','MarkerFaceColor','r','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none')
%         end
%         for i3=1:size(ErrorM,1)
%             figure(5), hold on
%             scatter(PredictTimeM/CC.lyapTime,ErrorM(i3,:),'o','MarkerFaceColor','b','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none')
%         end
%         for i3=1:size(ErrorM,2)
%             figure(i3+10), hold on
%             histogram(ErrorM(:,i3))
%             pav= sprintf('Lambda*t_{pred} = %.3f',PredictTimeM(i3)/CC.lyapTime);
%             title(pav)
%         end

    end
end

% lg={};
% figure(4), hold on
% %for every degree
% for ideg = 1:length(deg)
%     %for every beta
%     for ii1 = 1:length(fpavM)
%         tmpvid = reshape(vid(ideg,ii1,:),[],1,1);
%         tmpvidlog = reshape(vidlog(ideg,ii1,:),[],1,1);
%         tmpminver = reshape(minVer(ideg,ii1,:),[],1,1);
%         tmpmaxver = reshape(maxVer(ideg,ii1,:),[],1,1);
%         tmpstdver = reshape(stdVer(ideg,ii1,:),[],1,1);
%         errorbar(PredictTimeM/CC.lyapTime, log(tmpvid),tmpstdver,'gs-')
%         plot(PredictTimeM/CC.lyapTime, tmpvidlog,'bd-')
%         lg{ii1}=sprintf('$\\log(\\beta)=-%d$',fpavM(ii1));
%     end
% end

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
%         plot(PredictTimeM/CC.lyapTime, tmpvid)
        %plot(PredictTimeM/CC.lyapTime, tmpup)
        %plot(PredictTimeM/CC.lyapTime, tmpdown)
        shade(PredictTimeM/CC.lyapTime, tmpvid,'b--',PredictTimeM/CC.lyapTime, tmpup,'b-',PredictTimeM/CC.lyapTime, tmpdown,'b-','FillType',[2 3])

        %errorbar(PredictTimeM/CC.lyapTime, tmpvid,tmpstdver)
        lg{ii1}=sprintf('$\\log(\\beta)=-%d$',fpavM(ii1));
    end
end
% set(gca,'yScale','log')


figure(2), hold on
%for every degree
for ideg = 1:length(degM)
    %for every beta
    for ii1 = 1:length(fpavM)
        tmpvid = reshape(countN(ideg,ii1,:),[],1,1);
        plot(PredictTimeM/CC.lyapTime, tmpvid)
        lg{ii1}=sprintf('$\\log(\\beta)=-%d$',fpavM(ii1));
    end
end
legend(lg,'Interpreter','latex')
xlabel('$\Lambda t_{pred}$')
ylabel('NaNs')



% set(gca,'yScale','log')
xlabel('$\Lambda t_{pred}$','Interpreter','latex')
ylabel('Error')

legend({'Monom','Monom'})


figure(2), hold on
%for every degree
for ideg = 1:length(degM)
    %for every beta
    for ii1 = 1:length(fpavM)
        tmpvid = reshape(countN(ideg,ii1,:),[],1,1);
        plot(PredictTimeM/CC.lyapTime, tmpvid)
        lg{ii1}=sprintf('$\\log(\\beta)=-%d$',fpavM(ii1));
    end
end
legend(lg,'Interpreter','latex')
xlabel('$\Lambda t_{pred}$')
ylabel('NaNs')
legend({'Monom','Monom'})





