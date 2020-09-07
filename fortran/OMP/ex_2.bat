@echo off

set Tfin=1e-2
set dt=1e-3
set NPaths=1e5
set Nthreads=6

call result.exe method=ssTauLeap              Tfin=%Tfin% dt=%dt% X0=1e3,1e3,1e6 rates=1e3,1e3,1e-5,10,1,1e6 prefix=ex_2 writeFinalTime=F NPaths=%NPaths% Nthreads=%Nthreads%
call result.exe method=ThetaTauLeap theta=1.0 Tfin=%Tfin% dt=%dt% X0=1e3,1e3,1e6 rates=1e3,1e3,1e-5,10,1,1e6 prefix=ex_2 writeFinalTime=F NPaths=%NPaths% Nthreads=%Nthreads%
call result.exe method=ThetaTauLeap theta=0.5 Tfin=%Tfin% dt=%dt% X0=1e3,1e3,1e6 rates=1e3,1e3,1e-5,10,1,1e6 prefix=ex_2 writeFinalTime=F NPaths=%NPaths% Nthreads=%Nthreads%
call result.exe method=SSA                    Tfin=%Tfin% dt=%dt% X0=1e3,1e3,1e6 rates=1e3,1e3,1e-5,10,1,1e6 prefix=ex_2 writeFinalTime=F NPaths=1e5 Nthreads=%Nthreads%	