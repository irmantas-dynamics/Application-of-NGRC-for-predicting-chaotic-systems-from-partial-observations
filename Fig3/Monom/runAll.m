%array of polynomial degrees
degM=[14];

%array of regression parameter degrees
betM = [3,5,7,9];

for i1 = 1:length(degM)
    for i2 =1:length(betM)
        deg = 0:degM(i1);
        betdeg = betM(i2);
        tic
        run_program_function(deg,betdeg);
        toc
    end
end
