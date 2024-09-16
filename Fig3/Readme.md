'*draw_quartiles.m*'  - draws Fig.3 from the paper. The script uses data from  files 'Cheb14_bet9.mat' and 'Monom14_bet9.mat'

Function '*runBoth.m*' calls scripts '*Cheb/runAll.m*' and '*Monom/runAll.m*', which generates data files with naming scheme '*ChebX_betY.mat*'. Here *X* marks the degree of polynomial and *Y* marks the degree of regression parameter *beta*. The degree of the polynomial is defined in the file '*runAll.m*' as variable *degM* and regression parameters are defined in the same script as *betM*. The comparison between RMSE of different calculations can be done with script '*plot)intermediate_results.m*'.

Warning: execution of "Cheb/runAll.m" and "Monom/runAll.m" is very time-consuming (i.e. may need a few hours) and also requires "Parallel Computing Toolbox".
