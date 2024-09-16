draw_quartiles.m  - draws Fig.3 from paper. Function uses results from  files 'Cheb14_bet9.mat' and 'Monom14_bet9.mat'

Function '*runBoth.m*' calls scripts '*Cheb/runAll.m*' and '*Monom/runAll.m*', which generates data files with naming scheme '*ChebX_betY.mat*'. Here *X* marks degree of polynome and *Y* marks degree of regression parameter *beta*. The degree of polynome is defined in the file '*runAll.m*' as variable *degM* and regression parameters are defined in the same script as *betM*. The comparisson between RMSE of different calculations can be done with script '*plot)intermediate_results.m*'.

Warning: execution of "Cheb/runAll.m" and "Monom/runAll.m" is very time consuming (i.e. may need few hours) also requires "Parallel Computing Toolbox".
