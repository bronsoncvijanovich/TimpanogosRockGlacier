clear; close; clc

fids = "rho"+[string([2325,2670])]+".txt";

for i = 1:length(fids)
    GTeC('inpgrav_xyz.csv',fids(i),'topomaxi_95m.dat','topomidi_47m.dat','topomini_19m.dat','topovemini_9m.dat','topomicro_elcorr.dat')
    load crtp.dat
    tail = split(fids(i),'/');
    tail = tail(2);
    outfile = "CBA_"+tail;
    writematrix(crtp,outfile);
    close all
end