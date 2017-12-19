function patt = lag_patterns(bold_traces,nuis_traces,tr,ranges,exclude_cutoffs,xlab,ylab,ntr_post,nmin,doplot)

%lag_patterns  Visualize lagged correlates of nuisance traces in BOLD signal, binned into the ranges specified, as described in Byrge & Kennedy, Identifying and characterizing systematic temporally-lagged BOLD artifacts.
%
% SYNTAX: patt = lag_patterns(bold_traces,nuis_traces,tr,ranges,exclude_cutoffs,xlab,ylab,ntr_post,nmin,doplot)
% 
% ARGS:
% bold_traces 
%             - 2-dimensional array containing the BOLD traces to be
%               analyzed. Must be organized with different scans in rows
%               and TRs/time in columns. BOLD traces can be global BOLD
%               signal or from other spatial scale (ROI, voxel). Can be a
%               single BOLD trace from one scan, multiple traces from the
%               same subject, or multiple traces from many subjects. 
%               (Can also be an array of other time series, e.g. non-BOLD, 
%               in order to examine interrelationships among time series.)
%		Input traces are z-scored (if not already z-scored).
%
% nuis_traces 
%             - 2-dimensional array containing nuisance traces to be
%               analyzed. Must be organized exactly the same as bold_epochs 
%               (e.g. rows must correspond to the same scans). Any nuisance 
%               measure sampled at (or interpolated to) the same rate as 
%               the BOLD traces can be used. Framewise displacement was 
%               used primarily in Byrge & Kennedy. Input traces
%		are assumed to have meaningful scale and not z-scored.
%
% tr      (optional)   
%            - TR for BOLD traces, in seconds. Used for visualization
%               only. If specified, plots are converted into seconds instead 
%               of TRs. (default is to plot in TR space.)
%
% ranges  (optional)
%            - How to subdivide the nuisance trace into ranges. 
%              Either:
%                  * number of percentiles to use for data-driven ranges
%              Or:
%                  * n X 2 matrix - hard-coded min and max for each of n 
%                    ranges. (default is 14-hard-coded framewise 
%                    displacement (FD) ranges from Byrge & Kennedy;
%                    default is not appropriate for nuisance 
%                    traces other than FD.)
%  
% exclude_cutoffs (optional)
%           - [ min_cutoff max_cutoff ] - epochs beginning with nuisance 
%             values <= min_cutoff or >= max_cutoff will be excluded from 
%             analysis. Appropriate for excluding e.g. excessive motions.
%             Recommended to exclude 0 if using framewise displacements.
%             Use e.g. [ NaN max_cutoff ] to specify one but not both
%             cutoffs.  (default: no cutoffs). 
%
% xlab (optional)
%           - string to use as x label text. 
%             (default: 'TRs post-displacement')
%
% ylab (optional)
%           - string to use as y label text. 
%             (default: 'BOLD (z)')
%
% ntr_post (optional)
%           - # of TRs after initial displacement (or other nuisance event)
%             to plot (default: 65)
%
% nmin (optional)
%           - number of instances of each range to require before omitting 
%             range. Not used for data-driven ranges. Most relevant when 
%             analyzing individual scans. (default: 20)
%
% doplot (optional)
%           - whether to generate the plot (default: true)
%
% note: use [] for defaults (for optional args). 
%           
%
%
% RETURNS:
% patt  
%       - ranges X ntr_post array containing mean lagged pattern for each
%         range. This is the array plotted in the figure (before optional
%         conversion from TRs to seconds).
%
%
% AUTHOR : Lisa Byrge, Ph.D.
% PLACE  : Indiana University, Social Brain Lab (PI: Dan Kennedy, Ph.D.)
% DATES  : 7/11/2017 LB initial version.
%
% Copyright 2017 Indiana University.
% All rights reserved. 
%
%
% EXAMPLES:
%
% lag_patterns(bold_traces, nuis_traces);
%       -- generates plot using default settings in TR space
%          ( 14 hard-coded ranges, no max values excluded, 65 TRs
%          post-displacement)
%
% lag_patterns(bold_traces, nuis_traces, .720);
%       -- generates plot using default settings as above, but converted
%          into seconds
%
% lag_patterns(bold_traces, nuis_traces, .720, 10);
%      -- generates plot using 10 data-driven ranges instead of hard-coded
%         ranges, and all other default values
%
% lag_patterns(bold_traces, nuis_traces, .720, [], [0 1.5]);
%       -- generates plot using default 14 hard-coded ranges but 
%          excludes any nuisance values exceeding 1.5. 
%
% lag_patterns(nuis1_traces, nuis2_traces, .720, 3, [],'TRs post-nuisance1','nuisance2 (z)');
%       -- generates plot showing relationships between nuis1 and nuis2 traces,
%          with epochs of the nuis1 trace plotted, and the nuis2 trace used
%          for ranges (3, percentile-based), and custom axis labels indicated.
%
% patt = lag_patterns(bold_traces, nuis_traces, [], [], [], [], [], [], [], false);
%       -- suppress plot and return the pattern that is normally plotted instead, 
%          with all default options.



if nargin<2
    error('Not enough input arguments');
end
if nargin<3 || length(tr)==0
    tr=1;convert_to_s=false;
else
    convert_to_s=true;
end
if nargin<4 || length(ranges)==0
    ranges = [ [ 0:.05:.55 .7 .9  ] ; [ .05:.05:.55 .7 .9 1.5 ] ]' ;
end
if nargin<5 || length(exclude_cutoffs)==0
    exclude_cutoffs=[NaN NaN];
end
if nargin<6 || length(xlab)==0
    xlab=[]; % set below
end
if nargin<7 || length(ylab)==0
    ylab='BOLD (z)';
end
if nargin<8 || length(ntr_post)==0
    ntr_post=65;
end
if nargin<9 || length(nmin)==0
    nmin=20;
end
if nargin<10 || length(doplot)==0
    doplot=true;
end

if tr>5
    error('Enter TR in seconds (not milliseconds)');
end

if not(length(exclude_cutoffs)==2)
    error('Invalid exclude_cutoffs spec.');
end

% number of TRs *prior* to the event of alignment (e.g. initial displacement)
% to visualize
ntr_pre=0; 

epoch_length=ntr_pre+ntr_post;

n_scans=size(bold_traces,1);
scan_length=size(bold_traces,2);

if size(bold_traces)~=size(nuis_traces)
    error('bold_traces and nuis_traces must have same size');
end

if size(nuis_traces,1)*(size(nuis_traces,2)-epoch_length) <= 0
    error('ERROR: make sure input matrices are organized as scans (rows) X TRs (columns). Must have more columns than epoch_length (default: 65).');
end

if size(bold_traces,1)>size(bold_traces,2)
    disp('WARNING: make sure input matrices are organized as scans (rows) X TRs (columns).');
    disp('Press any key to continue... (or ctrl-c to quit and transpose input matrices)');
    pause;
end



    
% exclude extreme values if indicated. omit final epoch_length TRs as they
% do not lead a complete epoch.
nuis_traces_flat=reshape(nuis_traces(:,1:(size(nuis_traces,2)-epoch_length)),size(nuis_traces,1)*(size(nuis_traces,2)-epoch_length),1);
nuis_traces_flat(nuis_traces_flat<=exclude_cutoffs(1))=NaN;
nuis_traces_flat(nuis_traces_flat>=exclude_cutoffs(2))=NaN;

counts=nan(size(ranges,1),1);

if length(ranges)==1 
    % generate percentile-based bins 
    
    n_prc=ranges;
    ranges=nan(n_prc,2);
    
    for iprc=1:n_prc
        if iprc==1
            ranges(iprc,1)=nanmin( nuis_traces_flat ) ;
        else
            ranges(iprc,1)=ranges(iprc-1,2);
        end
        ranges(iprc,2)=prctile( nuis_traces_flat, (100/n_prc)*iprc );

	counts(iprc)=length(find(nuis_traces_flat >=ranges(iprc,1) & nuis_traces_flat < ranges(iprc,2) ));
    end

else if length(ranges)>1 
        % if hard-coded ranges specified, make sure enough data to use them
        
        if size(ranges,2)~=2
            error('Invalid ranges spec.');
        end
                
        if ranges(1,1)<=exclude_cutoffs(1) || ranges(end,2) >= exclude_cutoffs(2)
           
            ranges(1,1)=nanmin(nuis_traces_flat);
            ranges(end,2)=nanmax(nuis_traces_flat);
                
            %disp('Ranges adjusted to exclude specified endpoints...');
            %disp(ranges');      
        end       

        for irange=1:size(ranges,1)
            
            if irange<size(ranges,1)
                counts(irange)=length ( find ( nuis_traces(:) >= ranges(irange,1) & nuis_traces(:) < ranges(irange,2)  ) );
            else
                counts(irange)=length ( find ( nuis_traces(:) >= ranges(irange,1) & nuis_traces(:) <= ranges(irange,2)  ) );
            end
            if sum(counts<nmin)>0
                disp('Warning: insufficient data for ranges spec. Results for some ranges will be missing. (To fix, widen range, decrease nmin, or try percentile-based bins.)');
            end
        end        
    else
        error('Invalid ranges spec.');
    end
end


ind_mat = nan(scan_length-epoch_length,epoch_length+1);
for irow=1:size(ind_mat,1)
    ind_mat(irow,:)=irow:(epoch_length+irow);
end


gb_needed=n_scans*size(ind_mat,1)*(epoch_length+1)*8*(9.313225746e-10);
if gb_needed>.5
    fprintf('WARNING: this will require at least %0.2f GB memory ... ',gb_needed);
    disp('Press any key to continue... (or ctrl-c to quit and re-run on a subset of input data)');
    pause;
end


if abs(nanmean(nanmean(bold_traces,2)))>.005 || abs(nanmean(nanstd(bold_traces'))-1)>.005
    %disp('z-scoring input BOLD traces....');
    bold_traces=nanzscore(bold_traces')';
end

bold_epochs=nan(n_scans*size(ind_mat,1),epoch_length+1);
%nuis_epochs=nan(n_scans*size(ind_mat,1),epoch_length+1); 
initial_nuis_val=nan(n_scans*size(ind_mat,1),1); 

scan_idx=nan(n_scans*size(ind_mat,1),1); % for reconstructing if needed
tr_idx=nan(n_scans*size(ind_mat,1),1); % for reconstructing if needed

% generate epochs from BOLD signal and save initial nuisance value for each
% epoch. note that bold_epochs(:,1) corresponds to the TR of initial
% nuisance value (corresponding to t=0 on the plots), and then 
% bold_epochs(:,2:end) corresponds to the TRs *following* the initial 
% nuisance value. 

for irow=1:size(ind_mat,1)
    ind=ind_mat(irow,:);
    bold_epochs(irow:size(ind_mat,1):n_scans*size(ind_mat,1),:) = bold_traces(:,ind);
    %nuis_epochs(irow:size(ind_mat,1):n_scans*size(ind_mat,1),:) = nuis_traces(:,ind);
    initial_nuis_val(irow:size(ind_mat,1):n_scans*size(ind_mat,1),1) = nuis_traces(:,ind(ntr_pre+1));

    scan_idx (irow:size(ind_mat,1):n_scans*size(ind_mat,1) ) = 1:n_scans;
    tr_idx (irow:size(ind_mat,1):n_scans*size(ind_mat,1) ) = irow;
end

% mean pattern following nuisance values within each range
patt=nan(length(ranges),epoch_length+1);
for irange=1:length(ranges);
   if counts(irange)>nmin
      patt(irange,:)=nanmean(bold_epochs( find(  ( initial_nuis_val >=ranges(irange,1) & initial_nuis_val < ranges(irange,2))), :),1);
   end
end


if doplot    
    f=figure;
    set(f,'defaultAxesFontSize',24);
    subplot(2,3,[1 2 4 5]);
    cmap=parula(length(ranges));
    hold on;
    
    % shift x-axis so bold_epochs(:,1) which corresponds to the initial
    % nuisance value is plotted as t=0. 
    x_axis=(1:(epoch_length+1))-1-ntr_pre; 
    if convert_to_s
        x_axis=x_axis*tr;
        if length(xlab)==0
            xlab='seconds post-displacement';
        end
    else
        if length(xlab)==0
            xlab='TRs post-displacement';
        end
    end
    
    for irange=1:length(ranges);
        plot(x_axis,patt(irange,:),'LineWidth',3,'col',cmap(irange,:));
    end
    xlim([nanmin(x_axis) nanmax(x_axis)]);
    set(gca,'XGrid','on');
    set(gca,'YGrid','on');
    set(gca,'GridLineStyle',':');
    xlabel(xlab);
    ylabel(ylab);
     
    subplot(2,3,[3 6]);
    hold on;
    for irange=1:length(ranges);
        plot(NaN,NaN,'LineWidth',10,'Color',cmap(irange,:));
        if irange==1
            leg{irange}=sprintf('    < % 0.2f',ranges(irange,2));
        else if irange==length(ranges)
                leg{irange}=sprintf('    > % 0.2f ',ranges(irange,1));
            else
                leg{irange}=sprintf('    % 0.2f  to % 0.2f ',ranges(irange,1),ranges(irange,2));
            end
        end
    end
    axis off;
    legend(leg);
    
    set(f,'Position',[446   605   851   493]);
end

end


function zx = nanzscore(x)
    zx = bsxfun(@rdivide, bsxfun(@minus, x, nanmean(x)), nanstd(x));
end
