
%
% 07-04-2014
%
% Implemented by Leonardo S. Barbosa : lsbarbosa@gmail.com
%
% Calculate statistics per sample/electrode and then perform a simple
% clustering / permutation statistics folowing
%
% Maris, E., & Oostenveld, R. (2007). 
%     Nonparametric statistical testing of EEG- and MEG-data. Journal of neuroscience methods, 164(1), 177–90. doi:10.1016/j.jneumeth.2007.03.024
% 

%% load data 

% data = struct;
% data.time = D.time;
% data.elec = leftlrpidx;
% data.elecnames = {'C3', 'P3'};
% data.raw = ft_lrp_allsubj(1:2);
% data.label = {'Condition 1', 'Condition 2'};
% data.ylim = {[-2 2] [-2 2]};
% data.xTick = [-3200 -2400 -1600 -800 0 800 1600 2400 3200 4000];
% data.colormap = zmap;
% % data.xTick = [-1600 -800 0 800 1600 2400 3200 4000 4800 5600];
% save('data', 'data');

load('data');

% data.time : Array with time (in seconds) for each sample
% data.elec : electrodes to select
% data.raw  : one cell where each entry conatin one condition represented by a matrix of [electrodes samples subjects]

%% set parameters

% Cluster / Permutation parameters
montecarloalpha=0.05;
% montecarloalpha=0.1;

% clusteralpha=0.15;
clusteralpha=0.1;
% clusteralpha=0.05;
% tinv(clusteralpha, nsuj-1)

% BE CAREFUL TO NOT FULL YOUR MEMORY
% permutation monte carlo statistics 
%         nsim = 2000; 
%         nsim = 1500; 
nsim = 1000; 
%         nsim = 500; 
% NO MEMORY
%         nsim = 200; 
%         nsim = 40; 

% Plot parameters

if isfield(data, 'xTick'), xTick = data.xTick; else xTick = 1000*linspace(1,data.time(end),10); end
%         c={[1 0 1], [0 1 0], [1 1 0], [0 1 1], [1 0 0], [0 0 1], [0 0 0], [.5 0 .5], [0 .5 0], [.5 .5 0], [0 .5 .5], [.5 0 0], [0 0 .5]};
c={[0 0 0], [.5 .5 .5], [0 1 1], [1 0 0],  [1 1 1], [.8 .8 .8], [.5 .5 .5]};
linetype={'-', '--', '-.', ':', '-', '--', '-.'};
smoothing_method= 'moving';
flag_fastsmooth =  false; 
% smoothwindow = 500;
smoothwindow = 250;


%% Create raw data images

nelec = length(data.elec);
ncond = length(data.raw);
ntime = size(data.raw{1},2);
nsuj = size(data.raw{1},3);

bdata = zeros(nelec, ntime, nsuj, nsim, ncond);
rd = cell(1,ncond);

f=zeros(1,ncond);
h=cell(1,ncond);
tic

milliTime = data.time*1000; 

for iCond=1:ncond

    rdm = data.raw{iCond}(data.elec,:,:); % region of interest 
    [~,~,~,STATS] = ttest(rdm, 0, clusteralpha, 'both', 3);
    rd{iCond} = STATS.tstat;

    % Bootstrap
    fprintf('Calculating permutations...');
    tic;
    for ib =1:nsim
        % permute random number of subjects
        p = 2*(rand(1,nsuj) > .5)-1;
        bp = permute(repmat(p, [nelec 1 ntime]), [1 3 2]);

        bdata(:,:,:,ib,iCond) = bp .* rdm;
    end
    toc;

    % plot selected electrodes
    f(iCond) = figure('Color',[1 1 1]);
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) pos(3) * 1.2105 pos(4) * 1.2105]);
%             set(gcf, 'Position', [pos(1) pos(2) pos(3) * 1.8157 pos(4) * 1.8157]);
    hElec = zeros(1,nelec);
    for ielec=1:nelec

        % permutation mean
%         me = mean(bdata(ielec, :, :, :, iCond),3);
%         bmean = squeeze(mean(me,4));

        % plot raw data and t-values
        xlimVal = [min(milliTime)+100 max(milliTime)]; plot([xlimVal(1) xlimVal(2)], [0 0], '--k','linewidth', 0.5); 
        hold on;
        plot([0 0], data.ylim{iCond}, '--k','linewidth', 0.5);
        axis([xlimVal data.ylim{iCond}]); 
        V=axis ; 
        set(gca, 'FontSize',12); 
        set(gca, 'FontWeight', 'bold');
        legend off; 
        set(gcf,'color',[1 1 1]); 
        set(gca, 'Box', 'off'); 
        %                 set(gca, 'YTick', (-0.8:0.2:0.8)*iCond, 'XTick' , xTick);
        if iCond < 3
            spacing = 0.4;
        else
            spacing = 0.8;
        end
        set(gca, 'YTick', data.ylim{iCond}(1):spacing:data.ylim{iCond}(2), 'XTick' , xTick);


        % Smooth for visualization (not sure this should be done anyway)        
        data4plot = squeeze(mean(rdm(ielec,:,:),3));
        if smoothwindow
            if flag_fastsmooth
                splotdata = fastsmooth(data4plot,smoothwindow,1,1);
%                 sbmean = fastsmooth(bmean,smoothwindow,1,1);
            else
                splotdata = smooth(data4plot,smoothwindow,smoothing_method)';
%                 sbmean = smooth(bmean,smoothwindow,smoothing_method)';
            end
        else
            splotdata = data4plot;
%             sbmean = bmean;
        end

        % Plot scatter plot
%         if nelec < 3
%             plot(milliTime, splotdata, linetype{ielec},  'Color', [0 0 0], 'LineWidth', 10);
%         else
%             plot(milliTime, splotdata, 'Color', [0 0 0], 'LineWidth', 10);
%         end

        hElec(ielec) = scatter(milliTime(1:end-50), splotdata(1:end-50), 105, c{ielec}, 'filled');

        hold on;
        % Plot Sample Statistics
        scatter(milliTime(1:end-50), splotdata(1:end-50), 25, rd{iCond}(ielec, 1:end-50), 'filled');

%         % plot permutation mean
%         scatter(milliTime, sbmean, 5, [.5 .5 .5], 'filled');

%         bstd = squeeze(std(me,[],4));
%         scatter(milliTime, 2*bstd, 10, 'k', 'filled')
%         scatter(milliTime, -2*bstd, 10, 'k', 'filled')


    end

    set(gca,'ydir', 'reverse');
    th = title(data.label{iCond});
    set(th, 'FontSize',14); 
    set(th, 'FontWeight', 'bold');

    colormap(data.colormap);
    hc = colorbar;
    caxis([-4 4])

    xlh = xlabel('ms');
    set(xlh, 'FontSize',12); 
    set(xlh, 'FontWeight', 'bold');

    ylh = ylabel('\muV');
    set(ylh, 'FontSize',12); 
    set(ylh, 'FontWeight', 'bold');

    clh = get(hc,'ylabel');
    set(clh,'String', 't-value');
    set(clh, 'FontSize',12); 
    set(clh, 'FontWeight', 'bold');

    h{iCond} = hElec;
end

toc;


%% single pixel statistics of permutations

tic;

nsim = size(bdata,4);

df = size(bdata,3)-1;

pd = cell(1,ncond);
for iCond = 1:ncond

    if ~f(iCond)
        fprintf('Empty condition %d\n', iCond);
        continue;
    end  
    
    fprintf('Calculating single pixel statistics for condition %s...', data.label{iCond});

    tic;
    pdm = bdata(:,:,:,:,iCond);
    [~,~,~,STATS] = ttest(pdm, 0, clusteralpha, 'both', 3);
    pd{iCond} = STATS.tstat;

%     pd{iCond} = squeeze(STATS.tstat);
%     % avoid problems when there is only one channel (average)
%     pd{iCond} = reshape(STATS.tstat,[size(STATS.tstat,1),size(STATS.tstat,2),size(STATS.tstat,4)]);
    toc;
end
    

%% cluster statistics

for iCond = 1:ncond

    if ~f(iCond)
        fprintf('Empty condition %d\n', iCond);
        continue;
    end    
    
    fprintf( 'Computing significance for %s\n', data.label{iCond});
    [realpos, realneg] = findcluster(rd{iCond}, df, clusteralpha);

    nelec = size(rd{iCond}, 1);

%     if nelec < 7 % number of pure colors
%         cont = contrast(i);
%     else
        cont = 1;
%     end  
            
    for ie = 1:nelec
        realpos{ie}.pmonte = zeros(size(realpos{ie}.tclusters));
        realneg{ie}.pmonte = zeros(size(realneg{ie}.tclusters));
    end

    for isim = 1:nsim
        [simpos, simneg] =  findcluster(pd{iCond}(:,:,1,isim), df, clusteralpha);
        for ie = 1:nelec

            maxval = max(simpos{ie}.tclusters);
            if ~isempty(maxval)
                realpos{ie}.pmonte = realpos{ie}.pmonte + (realpos{ie}.tclusters < maxval)./nsim;
            end

            minval = min(simneg{ie}.tclusters);
            if ~isempty(minval)
                realneg{ie}.pmonte = realneg{ie}.pmonte + (realneg{ie}.tclusters > minval)./nsim;
            end
        end
    end

    figure(f(iCond));
    hold on;

    mylim = ylim;
    
    hElec = [];
    lElec = {};
    barpos = 1.05;
    barsep = 15;
    
    for ie = 1:nelec
        pmonte = realpos{ie}.pmonte;
        goodc = find(pmonte < montecarloalpha);
        contrast = linspace(.5, 1, length(goodc));
        ht = 0;
        for i = 1:length(goodc)
            ic = goodc(i);
            samples = realpos{ie}.clusters == ic;
            cint = [min(milliTime(samples)) max(milliTime(samples))];
            [peakv,peaki] = max(rd{iCond}(ielec,samples));
            cintsamples = find(samples);
            peakt = milliTime(cintsamples(peaki));
            fprintf('\telec %d | pos | p-value : %0.4f | time :  %1.0f %1.0f [peak : %1.0f]; ... \n', ie, pmonte(ic), cint, peakt);
            goodtime = find(realpos{ie}.clusters==ic);
            x = milliTime(goodtime);
            s = mylim(2)/barsep;
            p = mylim(2)/barpos;
            y = (p + (1-ie)*s) * ones(1,length(goodtime));
            
            ht = scatter(x,y,25,cont*c{ie});
            if cont < 1
                hElec = [hElec ht];
                lElec = [lElec sprintf('%s pos : %0.4f', elecnames{ie}, pmonte(ic))];
            end
            
        end
        
    end

    fprintf('\n');
    
    for ie = 1:nelec
        pmonte = realneg{ie}.pmonte;
        goodc = find(pmonte < montecarloalpha);
        contrast = linspace(.5, 1, length(goodc));
        ht = 0;
        for i = 1:length(goodc)
            ic = goodc(i);
            samples = realneg{ie}.clusters == ic;
            cint = [min(milliTime(samples)) max(milliTime(samples))];
            [peakv,peaki] = min(rd{iCond}(ielec,samples));
            cintsamples = find(samples);
            peakt = milliTime(cintsamples(peaki));
            fprintf('\telec %d | neg | p-value : %0.4f | time :  %1.0f %1.0f [peak : %1.0f]; ... \n', ie, pmonte(ic), cint, peakt);
            goodtime = find(realneg{ie}.clusters==ic);
            x = milliTime(goodtime);
            s = mylim(1)/15;
            p = mylim(1)/1.1;
            y = (p + (1-ie)*s) * ones(1,length(goodtime));

            ht = scatter(x,y,25,cont*c{ie});
            if cont < 1
                hElec = [hElec ht];
                lElec = [lElec sprintf('%s pos : %0.4f', elecnames{ie}, pmonte(ic))];
            end
            
        end
        
    end

    if isempty(hElec)
        hElec = h{iCond};
        lElec = data.elecnames;
    end
        
    hLegend = legend(hElec,lElec, 'Location', 'Best');
    
    % fill markers in the legend with the correct colors
    hMarkers = findobj(hLegend,'type','patch');
    for im = 1:length(hMarkers)
        hc = get(hMarkers(im), 'MarkerFaceColor');
        set(hMarkers(im), 'MarkerEdgeColor', hc, 'MarkerFaceColor', hc, 'MarkerSize', 9);
    end

end

toc;


