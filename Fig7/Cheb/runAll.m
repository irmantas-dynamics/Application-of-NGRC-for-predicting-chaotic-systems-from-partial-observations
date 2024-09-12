degM=[21];

betM = [9,7,5,3];

for i1 = 1:length(degM)
    for i2 =1:length(betM)
        deg = 1:2:degM(i1);
        betdeg = betM(i2);
        run_program_function(deg,betdeg);
    end
end
