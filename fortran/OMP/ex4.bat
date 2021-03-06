@echo off

rem xT = 200

rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=1,2,3e0,4e0 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6
rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=1,2,3e1,4e1 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6
rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=1,2,3e2,4e2 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6
rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=1,2,3e3,4e3 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6
rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=1,2,3e4,4e4 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6
rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=1,2,3e5,4e5 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6

rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=2,2,2e0,2e0 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6
rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=2,2,2e1,2e1 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6
rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=2,2,2e2,2e2 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6
rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=2,2,2e3,2e3 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6
rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=2,2,2e4,2e4 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6
rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=2,2,2e5,2e5 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6


rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=2e3,2e3,2e3,2e3 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6
rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=2e3,2e3,2e4,2e4 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6
rem call result.exe method=ssTauLeap theta=-1 Tfin=1000 dt=1 X0=10,0,0 rates=2e3,2e3,2e5,2e5 prefix=ex_3 writeFinalTime=T NPaths=100000 Nthreads=6

call result.exe method=ssTauLeap theta=-1 Tfin=100 dt=1 X0=100,100,0 rates=10,100 prefix=ex_4 writeFinalTime=T NPaths=100000 Nthreads=6
call result.exe method=ssTauLeap theta=-1 Tfin=100 dt=1 X0=100,100,0 rates=100,1000 prefix=ex_4 writeFinalTime=T NPaths=100000 Nthreads=6

rem call result.exe method=SSA Tfin=10 dt=1 X0=100,100,0 rates=10,100 prefix=ex_4 writeFinalTime=T NPaths=10000 Nthreads=6
rem call result.exe method=SSA Tfin=10 dt=1 X0=100,100,0 rates=100,1000 prefix=ex_4 writeFinalTime=T NPaths=10000 Nthreads=6

rem ################################

call result.exe method=ssTauLeap Tfin=1e0 dt=1e-1 X0=1e3,1e3,1e6 rates=1e3,1e3,1e-5,10,1,1e6 prefix=ex_3 writeFinalTime=F NPaths=1e4 Nthreads=6