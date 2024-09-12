clear all
close all

i2 = 3
deg = 21

fpav = sprintf('Cheb%d_bet%d.mat',deg,i2);
load(fpav)

nn = length(hTM);

histM = histM1;

dh =histSteps(2)-histSteps(1);
hh= histSteps(2:end)-dh;

t1 = CC.lyapTime;
t2 = 2*CC.lyapTime;

nt1 = round(t1/dh)+1;
nt2 = round(t2/dh)+1;


figure, hold on
for i1 = 1:nn
    plot(hh,histM(i1,:))
    vid = sum(hh.*histM(i1,:))/sum(histM(i1,:));
    fprintf('h = %f, vidurkis = %f\n\n',hTM(i1),vid)
end

legend(num2str(hTM))
pav=sprintf('Legend:step [bet=%.2e]',bet);
title(pav)


disp('Vidurkiai:')
for i1 = 1:nn
    vid = sum(hh.*histM(i1,:))/sum(histM(i1,:));
    fprintf('%f\t',vid)
end
fprintf('\n')



disp('1x Lyap time:')
for i1 = 1:nn
    vid = sum(histM(i1,nt1:end))/sum(histM(i1,:));
    fprintf('%f\t',vid)
end
fprintf('\n')

disp('2x Lyap time:')
for i1 = 1:nn
    vid = sum(histM(i1,nt2:end))/sum(histM(i1,:));
    fprintf('%f\t',vid)
end
fprintf('\n')

