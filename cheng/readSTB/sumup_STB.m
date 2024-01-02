clear 
clc

fin = 'D:\Chronos_Footage\1.1g\output_nnn';
fnum = 30;

tracer = [];
f = waitbar(0,'Please wait...');
for i = 1:fnum
    waitbar(i/fnum, f, 'Please wait...');
    t = load([fin '\tracer_STB_' num2str(i) '.mat']);
    tracer = [tracer,t.tracer];
end
waitbar(1, f, 'Saving...');
save([fin '\tracer_STB.mat'],'tracer', '-v7.3')
waitbar(1, f, 'Done!');
pause(0.5)
close(f)