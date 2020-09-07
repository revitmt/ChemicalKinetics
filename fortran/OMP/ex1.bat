@echo off

set Npaths=1e6
set Tfin=10

call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=10,0,0,0 rates=1e2,1e2,1e2,1e2,1e2,1e2 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=10,0,0,0 rates=1e4,1e4,1e4,1e4,1e4,1e4 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=10,0,0,0 rates=1e6,1e6,1e6,1e6,1e6,1e6 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc

call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=100,0,0,0 rates=1e2,1e2,1e2,1e2,1e2,1e2 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=100,0,0,0 rates=1e4,1e4,1e4,1e4,1e4,1e4 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=100,0,0,0 rates=1e6,1e6,1e6,1e6,1e6,1e6 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc

call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e2,1e2,1e2,1e2,1e2,1e2 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e4,1e4,1e4,1e4,1e4,1e4 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e6,1e6,1e6,1e6,1e6,1e6 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc

call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e0,1e-1,1e0,1e0,1e-1,5e-1 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e4,1e3,1e4,1e4,1e3,5e3 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc

rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e2,1e2,1e2,1e2,1e2,1e2 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e4,1e4,1e4,1e4,1e4,1e4 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e6,1e6,1e6,1e6,1e6,1e6 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e4,1e3,1e4,1e4,1e3,5e3 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc

rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e2,1e2,1e2,1e2,1e2,1e2 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e4,1e4,1e4,1e4,1e4,1e4 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e6,1e6,1e6,1e6,1e6,1e6 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=1000,0,0,0 rates=1e4,1e3,1e4,1e4,1e3,5e3 prefix=ex_1 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc


call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e-4,1e-4,1e-4,1e-4,1e-4,1e-4 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e-3,1e-3,1e-3,1e-3,1e-3,1e-3 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e-2,1e-2,1e-2,1e-2,1e-2,1e-2 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e-1,1e-1,1e-1,1e-1,1e-1,1e-1 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e0,1e0,1e0,1e0,1e0,1e0 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e1,1e1,1e1,1e1,1e1,1e1 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e2,1e2,1e2,1e2,1e2,1e2 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e3,1e3,1e3,1e3,1e3,1e3 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
call result.exe method=ssTauLeap Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e4,1e4,1e4,1e4,1e4,1e4 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc

rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e-4,1e-4,1e-4,1e-4,1e-4,1e-4 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e-3,1e-3,1e-3,1e-3,1e-3,1e-3 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e-2,1e-2,1e-2,1e-2,1e-2,1e-2 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e-1,1e-1,1e-1,1e-1,1e-1,1e-1 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e0,1e0,1e0,1e0,1e0,1e0 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e1,1e1,1e1,1e1,1e1,1e1 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e2,1e2,1e2,1e2,1e2,1e2 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e3,1e3,1e3,1e3,1e3,1e3 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=1 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e4,1e4,1e4,1e4,1e4,1e4 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc

rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e-4,1e-4,1e-4,1e-4,1e-4,1e-4 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e-3,1e-3,1e-3,1e-3,1e-3,1e-3 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e-2,1e-2,1e-2,1e-2,1e-2,1e-2 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e-1,1e-1,1e-1,1e-1,1e-1,1e-1 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e0,1e0,1e0,1e0,1e0,1e0 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e1,1e1,1e1,1e1,1e1,1e1 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e2,1e2,1e2,1e2,1e2,1e2 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e3,1e3,1e3,1e3,1e3,1e3 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=1e0 X0=250,250,250,250 rates=1e4,1e4,1e4,1e4,1e4,1e4 prefix=ex_1 writeFinalTime=T NPaths=%NPaths% Nthreads=6 optim_mode=unc