
data = struct;
data.time = D.time;
data.elec = leftlrpidx;
data.elecnames = {'C3', 'P3'};
data.raw = ft_lrp_allsubj(1:2);
data.label = {'Condition 1', 'Condition 2'};
data.ylim = {[-2 2] [-2 2]};
data.xTick = [-3200 -2400 -1600 -800 0 800 1600 2400 3200 4000];
data.colormap = zmap;
% data.xTick = [-1600 -800 0 800 1600 2400 3200 4000 4800 5600];
save('data', 'data');
