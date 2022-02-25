% Unlike corrplot, plotmatrix cannot receive a table. 
% Form a matrix of variables in each column.
data = rset(:,1:5); 
nVars = size(data,2); 
% Remove any rows that contain NaN values. Otherwise corr() will 
% return NaN. 
data(any(isnan(data),2), :) = []; 
% Create plotmatrix
figure('Name', 'Corrplot')
[sh, ax, ~, hh] = plotmatrix(data);
% Add axis labels
label = labels(B(:,idxLambda1SE)~=0);
label = label(1:5);
arrayfun(@(h,lab)ylabel(h,lab),ax(:,1), label')
arrayfun(@(h,lab)xlabel(h,lab),ax(end,:), labels)
% Compute correlation for each scatter plot axis
[r,p] = arrayfun(@(h)corr(h.Children.XData(:),h.Children.YData(:)),ax(~eye(nVars)));
% Label the correlation and p values
arrayfun(@(h,r,p)text(h,min(xlim(h))+range(xlim(h))*.05,max(ylim(h)),...
    sprintf('r=%.2f,  p=%.3f',r,p),'Horiz','Left','Vert','top','FontSize',8,'Color','r'),...
    ax(~eye(nVars)),r,p)
% Change marker appearance
set(sh, 'Marker', 'o','MarkerSize', 2, 'MarkerEdgeColor', ax(1).ColorOrder(1,:)) 
lsh = arrayfun(@(h)lsline(h),ax(~eye(nVars)));
% Add least square regression line. 
set(lsh,'color', 'm')