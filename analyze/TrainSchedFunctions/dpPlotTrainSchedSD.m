function dpPlotTrainSchedSD(data, varToPlot, axisLim)

% Plots std and sem data for experiment 'TrainSched' (Kristy's honors thesis project 2018-2019)
% data - data in a table format
% varToPlot - enter the dependent variable to plot as a string e.g. 'RT'
% axisLim - enter the axis limits in the format [xmin xmax ymin ymax]

rotations = [30 45 60];
color = {'b','r','g','m'};

% CREATE GROUP MEAN AND SEM
stdData = varfun(@nanstd, data, 'GroupingVariables', {'Group','abs_tgt_rot','TN_tgt'},'OutputFormat','table');
semData = varfun(@sem, data, 'GroupingVariables', {'Group','abs_tgt_rot','TN_tgt'},'OutputFormat','table');

% PLOT DATA
figure;
for ri = 1:length(rotations) ;
    for gi = unique(data.Group)';
        subplot(3,1,ri); hold on;
        
        idx = stdData.Group == gi & abs(stdData.abs_tgt_rot) == rotations(ri);
        
        nanstdVar = {strcat('nanstd_',varToPlot)};
        semVar = {strcat('sem_',varToPlot)};
        
        % Plot data
        x = stdData.TN_tgt(idx);
        y = stdData.(nanstdVar{1})(idx);
        err = semData.(semVar{1})(idx);
        shadedErrorBar(x, y, err, color{gi}, 1);
        
        % Title and axes labels
        axis(axisLim);
        str = sprintf('%iº Rotation', rotations(ri));
        title(str, 'interpreter','none');
        xlabel('Cycle Number'); % x-axis label
        
        ystr = sprintf('%s', varToPlot);
        ylabel(ystr, 'interpreter','none'); % y-axis label
        set(gca,'FontSize',13);
        
        % Draw lines
        drawline([5.5 15.5 105.5 115.5 145.5 155.5 165.5], 'dir', 'vert', 'linestyle', ':'); %between blocks
        drawline(0, 'dir', 'hor', 'linestyle', '-'); %0º hand angle
        
    end
end



end