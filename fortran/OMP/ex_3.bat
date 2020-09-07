@echo off

set dt=1e-1
set NPaths=1e6
set rates=3e-7,1e-4,1e-3,3.5


rem call result.exe method=ssTauLeap Tfin=10 dt=%dt% X0=1e5,100,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=SSA       Tfin=10 td=%dt% X0=1e5,100,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6

rem call result.exe method=ssTauLeap Tfin=10 dt=%dt% X0=1e5,200,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=SSA       Tfin=10 dt=%dt% X0=1e5,200,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6

rem call result.exe method=ssTauLeap Tfin=10 dt=%dt% X0=1e5,240,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=con
rem call result.exe method=SSA       Tfin=10 dt=%dt% X0=1e5,240,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6

rem call result.exe method=ssTauLeap Tfin=10 dt=%dt% X0=1e5,250,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=con
rem call result.exe method=SSA       Tfin=10 dt=%dt% X0=1e5,250,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6

rem call result.exe method=ssTauLeap Tfin=10 dt=%dt% X0=1e5,300,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=SSA       Tfin=10 dt=%dt% X0=1e5,300,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6

rem call result.exe method=ssTauLeap Tfin=10 dt=%dt% X0=1e5,400,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc
rem call result.exe method=SSA       Tfin=10 dt=%dt% X0=1e5,400,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6

rem call result.exe method=ssTauLeap_old Tfin=10 dt=%dt% X0=1e5,240,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3_old writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ssTauLeap_old Tfin=10 dt=%dt% X0=1e5,250,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3_old writeFinalTime=F NPaths=%NPaths% Nthreads=6



rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=10 dt=%dt% X0=1e5,100,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=10 dt=%dt% X0=1e5,200,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=10 dt=%dt% X0=1e5,240,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=10 dt=%dt% X0=1e5,250,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=10 dt=%dt% X0=1e5,300,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=10 dt=%dt% X0=1e5,400,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6

rem call result.exe method=ThetaTauLeap theta=1 Tfin=10 dt=%dt% X0=1e5,100,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ThetaTauLeap theta=1 Tfin=10 dt=%dt% X0=1e5,200,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ThetaTauLeap theta=1 Tfin=10 dt=%dt% X0=1e5,240,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ThetaTauLeap theta=1 Tfin=10 dt=%dt% X0=1e5,250,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ThetaTauLeap theta=1 Tfin=10 dt=%dt% X0=1e5,300,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ThetaTauLeap theta=1 Tfin=10 dt=%dt% X0=1e5,400,2e5 rates=3e-7,1e-4,1e-3,3.5 prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6


rem set dt=0.5
rem call result.exe method=ssTauLeapBistable      Tfin=10 dt=%dt% X0=1e5,250,2e5 rates=%rates% prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc lin_mode=end
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=10 dt=%dt% X0=1e5,250,2e5 rates=%rates% prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ThetaTauLeap theta=1.0 Tfin=10 dt=%dt% X0=1e5,250,2e5 rates=%rates% prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6

rem set dt=0.1
rem call result.exe method=ssTauLeapBistable      Tfin=10 dt=%dt% X0=1e5,250,2e5 rates=%rates% prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6 optim_mode=unc lin_mode=end
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=10 dt=%dt% X0=1e5,250,2e5 rates=%rates% prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6
rem call result.exe method=ThetaTauLeap theta=1.0 Tfin=10 dt=%dt% X0=1e5,250,2e5 rates=%rates% prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6

set dt=1.0
rem call result.exe method=ssTauLeapBistable      Tfin=20 dt=%dt% X0=1e5,250,2e5 rates=%rates% prefix=ex_3 writeFinalTime=F NPaths=1e4 Nthreads=6 optim_mode=unc lin_mode=end
call result.exe method=ssTauLeap_old      Tfin=10 dt=%dt% X0=1e5,250,2e5 rates=%rates% prefix=ex_3 writeFinalTime=F NPaths=1e4 Nthreads=6 optim_mode=unc lin_mode=end
rem call result.exe method=ThetaTauLeap theta=0.5 Tfin=100 dt=%dt% X0=1e5,250,2e5 rates=%rates% prefix=ex_3 writeFinalTime=F NPaths=1e4 Nthreads=6

rem call result.exe method=SSA                    Tfin=10 dt=1e-1 X0=1e5,250,2e5 rates=%rates% prefix=ex_3 writeFinalTime=F NPaths=%NPaths% Nthreads=6