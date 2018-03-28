function [ subset ] = generateSubset( m, xnodes, unodes )
% subset = generateSubset( m, nodes )
%   m must contain fields S, Species, nx,  
%   Output arguments
%       subset
%           .plotDiagram(t,m,con,sim,opts)
%               Input arguments
%                   opts
%                       .expandNodes
%                       .eliminateInactiveNodes

if nargin < 3
    unodes = [];
end
if nargin < 2
    xnodes = [];
end

nx = m.nx;
nu = m.nu;
Sx = full(m.S);
statenames = {m.States.Name}';
inputnames = {m.Inputs.Name}';

% Determine the participation matrix for the inputs. This indicates which
% reactions (the rows) are affected by which inputs (the columns).
% If an input affects a reaction rate, then it is a "reactant" to that
% reaction.
ru = m.drdu(rand, rand(m.nx,1), rand(m.nu,1));
ru = ru ~= 0;
S = [Sx; -ru'];

% By default, include all states and inputs, plotting the entire network
if isempty(xnodes)
    xnodes = statenames;
end
if isempty(unodes)
    unodes = inputnames;
end

subset.states = xnodes;
subset.nx = length(subset.states);
[~,state_i] = ismember(subset.states, {m.States.Name});
subset.C1 = zeros(1,nx); subset.C1(state_i) = 1;

subset.inputs = unodes;
subset.nu = length(subset.inputs);
[~,input_i] = ismember(subset.inputs, {m.Inputs.Name});
subset.U1 = zeros(1,nu); subset.U1(input_i) = 1;

subset.N1 = [subset.C1 subset.U1];

% Reorder nodes so that they are in the same order as in the model
subset.states = statenames(subset.C1 ~= 0);
subset.inputs = inputnames(subset.U1 ~= 0);
subset.nodes = [subset.states; subset.inputs];

%   Find which reactions flow into and out of these nodes and the
%   reactions' stoichiometries
subset.S = S; % Add both states and inputs as nodes
isnodeinsubset = logical([subset.C1 subset.U1]);
subset.S(~isnodeinsubset,:) = 0;
isreactioninsubset = any(subset.S ~= 0,1);
% Reassign reactions from S to include states that are reactants or
% products of subset reactions
subset.Sexpanded = subset.S;
subset.Sexpanded(:,isreactioninsubset) = S(:,isreactioninsubset);

% Construct reduced stoichiometric matrix
[subset.Si,subset.Sj] = find(subset.S ~= 0);
subset.Sr = S(unique(subset.Si),unique(subset.Sj));
[isx,isu] = deal(false(nx+nu, 1));
isx(1:nx) = true;
isu(nx+1:end) = true;
subset.nxr = sum(isx(unique(subset.Si)));
subset.nur = sum(isu(unique(subset.Sj)));
subset.nr = size(subset.Sr,2);

[subset.Siexpanded,subset.Sjexpanded] = find(subset.Sexpanded ~= 0);
subset.Srexpanded = S(unique(subset.Siexpanded),unique(subset.Sjexpanded));
subset.nxrexpanded = sum(isx(unique(subset.Siexpanded)));
subset.nurexpanded = sum(isu(unique(subset.Sjexpanded)));

% Set reaction names
rnames = {m.Reactions(unique(subset.Sj)).Name}';
subset.rnames = rnames;

% Set state names to include reactants and products of reactions, in
% addition to user-specified nodes
isnodeinexpandedsubset = any(subset.Sexpanded ~= 0,2);
nodenames = [statenames; inputnames];
subset.nodesexpanded = nodenames(isnodeinexpandedsubset);
subset.N1expanded = 1*isnodeinexpandedsubset;

subset.plotReactions = @plotReactions;
subset.plotDiagram = @plotDiagram;

    function bg = plotDiagram(t,m,con,sim,opts)
        % bg = plotDiagram(t,m,con,sim,opts)
        
        if nargin < 5
            opts = struct('expandNodes',false,'eliminateInactiveNodes',false);
            if nargin < 4
                sim = getSim(t,m,con);
            end
        end
        
        if ~isfield(opts, 'expandNodes')
            opts.expandNodes = false;
        end
        if ~isfield(opts, 'eliminateInactiveNodes')
            opts.eliminateInactiveNodes = false;
        end
        
        % Get simulation data
        [~,r] = getReactionRates(t,sim,m,con);
        x = getStateValues(t,sim,opts.expandNodes);
        u = getInputValues(t,sim,opts.expandNodes);
        xu = [x;u];
        
        
        nsims = size(r,2);
        
        % Add in null nodes or expanded nodes
        isnull = false(length(xu),1);
        if opts.expandNodes
            Srnull = subset.Srexpanded;
            nodesnullbase = subset.nodesexpanded;
        else
            Srnull = subset.Sr;
            nodesnullbase = subset.nodes;
            
            nrr = size(subset.Sr,2);
            %{
            nxur = size(subset.Sr,1);
            nullcounter = 1;
            for xui = 1:nxur
                isrxntoaddnull = all(Srnull([1:xui-1 xui+1:end],:) == 0,1);
                if any(isrxntoaddnull)
                    Srnull = [Srnull; zeros(1,nrr)];
                    Srnull(end,isrxntoaddnull) = -Srnull(xui,isrxntoaddnull);
                    nodesnullbase = [nodesnullbase; {sprintf('null%1d',nullcounter)}];
                    isnull = [isnull; true];
                    nullcounter = nullcounter + 1;
                end
            end
            %}
            hasonlyreactants = all(Srnull <= 0, 1);
            hasonlyinputs = all(Srnull(1:subset.nxr,:) == 0, 1);
            isnullrxn = hasonlyreactants & ~hasonlyinputs;
            nnullrxns = sum(isnullrxn);
            Srnull_cat = zeros(nnullrxns, nrr);
            indonlyreactants = find(isnullrxn);
            for i = 1:nnullrxns
                Srnull_cat(i,indonlyreactants(i)) = 1;
            end
            Srnull = [Srnull; Srnull_cat];
            isnull = [isnull; true(nnullrxns,1)];
            nodesnullbase = [nodesnullbase; arrayfun(@(i){sprintf('null%1d',i)}, (1:nnullrxns)')];
        end
        
        % Remove reactions that only involve inputs
        if opts.expandNodes
            reactionhasinputsonly = all(subset.Srexpanded(1:subset.nxrexpanded,:)==0, 1);
        else
            reactionhasinputsonly = all(subset.Sr(1:subset.nxr,:)==0, 1);
        end
        r(reactionhasinputsonly,:) = [];
        Srnull(:,reactionhasinputsonly) = [];
        rnames_ = subset.rnames;
        rnames_(reactionhasinputsonly) = [];
        
        nr = size(r,1);
        
        % Add in reaction nodes
        reactionsnullbase = textscan(sprintf('r%d\n',1:nr),'%s\n');
        reactionsnullbase = reactionsnullbase{1};
        
        %reactionsnullbase = repmat({''},nr,1);
        
        % Convert stoichiometry matrix to connectivity matrix
        Cbase = stoichiometric2connectivity(Srnull);
        
        % Initialize biographs
        temp = biograph(Cbase,[nodesnullbase;reactionsnullbase]);
        rindices_base = [NaN(size(nodesnullbase));(1:nr).'];
        bg = repmat(temp,nsims,1);
        
        % Set minimum line width and minimum line color
        linewidthscale = 30; % multiplier by which to multiply normalized fractions to determine line width. This is necessary because line widths less than 1 all look the same.
        minimumlinewidth = 0.01;
        minimumlinecolor = [0.7 0.7 0.7];
        aboveminimumlinecolor = [1 0 0];
        
        for si = 1:nsims

            nxunull = size(Srnull,1);
            nrnull = size(Srnull,2);
            
            C = Cbase; % reset C and states list in case some inactive nodes were eliminated
            nodesnull = nodesnullbase;
            reactionsnull = reactionsnullbase;
            rindices = rindices_base;
            
            % Get node values for all nodes for this simulation
            xu_si = zeros(nxunull,1);
            xu_si(~isnull) = xu(:,si);
            
            % Get reaction rates for all reactions for this simulation
            r_si = full(r(:,si));
            
            % Generate node labels as (1) names of species with concentration
            % values in square brackets for states and inputs and (2) reaction rate
            % values for reactions
            getlabel = @(name,value)[name ' [' sprintf('%1.2e',value) ']'];
            getlabel_r = @(value)sprintf('%1.2e',value);
            nodelabels = [
                cellfun(getlabel,nodesnullbase,num2cell(xu_si),'UniformOutput',false);
                cellfun(getlabel_r,num2cell(r_si),'UniformOutput',false)
                ];
            
            % Don't include concentrations for null nodes, for obvious
            % reasons
            nodelabels(isnull) = nodesnull(isnull);
            
            % Calculate edge weights and line widths

            weights_C = zeros(nxunull);
            linewidths_C = zeros(nxunull);
            
            for ri = 1:nrnull
                nodeisreactant = C(:,nxunull+ri) ~= 0;
                nodeisproduct = C(nxunull+ri,:) ~= 0;
                weights_C(nodeisreactant,nxunull+ri) = r_si(ri);
                weights_C(nxunull+ri,nodeisproduct) = r_si(ri);
                linewidths_C(nodeisreactant,nxunull+ri) = r_si(ri);
                linewidths_C(nxunull+ri,nodeisproduct) = r_si(ri);
            end
            
            %{
          risadded = false(nrnull,1);
            for ri = 1:nrnull
                % Skip this reaction if it was counted in an earlier
                % iteration
                if risadded(ri)
                    continue
                end
                sources = find(Srnull(:,ri) < 0);
                sinks = find(Srnull(:,ri) > 0);
                for soi = 1:length(sources)
                    thissource = sources(soi);
                    for sii = 1:length(sinks)
                        thissink = sinks(sii);
                        % Find all r's with the indicated source and sink
                        hasthesamesourcesandsinks = cellfun(@isequal,...
                            mat2cell(Srnull,nxnull,ones(nrnull,1)),...
                            repmat({Srnull(:,ri)},1,nrnull)...
                            );
                        risadded(hasthesamesourcesandsinks) = true;
                        thisr = sum(r(hasthesamesourcesandsinks,si));
                        %thisr_normalized = sum(rscaled_normalized(hasthesamesourcesandsinks,si));
                        % Switch source and sink if r is negative
                        if thisr < 0
                            thisr = -thisr;
                            thissinktemp = thissource;
                            thissourcetemp = thissink;
                            thisCelement = C(thissource,thissink);
                            % Add elements to C to ensure reverse direction
                            % reaction exists
                            C(thissourcetemp,thissinktemp) = C(thissourcetemp,thissinktemp) + thisCelement;
                        else
                            thissinktemp = thissink;
                            thissourcetemp = thissource;
                        end
                        weights_C(thissourcetemp,thissinktemp) = weights_C(thissourcetemp,thissinktemp) + thisr;
                        linewidths_C(thissourcetemp,thissinktemp) = linewidths_C(thissourcetemp,thissinktemp) + thisr;
%                         if linewidths_C(thissourcetemp,thissinktemp) <= minimumlinewidth
%                             linewidths_C(thissourcetemp,thissinktemp) = minimumlinewidth;
%                             linecolors_C{thissourcetemp,thissinktemp} = minimumlinecolor;
%                         else
%                             linecolors_C{thissourcetemp,thissinktemp} = aboveminimumlinecolor;
%                         end
                        % Remove edge for current direction
                        % If it doesn't exist, add edge for reverse direction

                    end
                end
            end
%}
            
            % Normalize line widths
            newtotalwidthamount = sum(sum(abs(linewidths_C)));
            if si == 1
                totalwidthamount = newtotalwidthamount;
            end
            linewidths_C_new = linewidths_C/newtotalwidthamount;
            linewidths_C = linewidths_C/totalwidthamount;
            
            
            % Set small line widths to minimum value
            isminlinewidth = linewidths_C <= minimumlinewidth;
            isminlinewidth_new = linewidths_C_new <= minimumlinewidth;
            linewidths_C(isminlinewidth) = minimumlinewidth;
            
            % Set line colors
            linecolors_C = repmat({aboveminimumlinecolor},nxunull+nrnull,nxunull+nrnull);
            linecolors_C(isminlinewidth) = {minimumlinecolor};
            linecolors_C = reshape(linecolors_C,nxunull+nrnull,nxunull+nrnull);
            
            % Eliminate inactive nodes, if requested
            if opts.eliminateInactiveNodes
                % A node is considered inactive if all reactions entering
                % and exiting the node are at the minimum line width,
                % normalizing to the total width in the current simulation
                % instead of to the base simulation to avoid global changes
                % from preventing reactions from being shown
                nodeisinactive = all(isminlinewidth_new,1)' & all(isminlinewidth_new,2);
                if si == 1
                    nodeisinactive_sim1 = nodeisinactive;
                end
                nodeisinactive = nodeisinactive & nodeisinactive_sim1; % only eliminate nodes that are inactive in both the original simulation and the current simulation
                C(nodeisinactive,:) = [];
                C(:,nodeisinactive) = [];
                nodesnull(nodeisinactive(1:nxunull)) = [];
                reactionsnull(nodeisinactive(nxunull+1:nxunull+nrnull)) = [];
                rindices(nodeisinactive) = [];
                nxunull = length(nodesnull);
                nrnull = length(reactionsnull);
                weights_C(nodeisinactive,:) = [];
                weights_C(:,nodeisinactive) = [];
                linewidths_C(nodeisinactive,:) = [];
                linewidths_C(:,nodeisinactive) = [];
                linecolors_C(nodeisinactive,:) = [];
                linecolors_C(:,nodeisinactive) = [];
                nodelabels(nodeisinactive) = [];
            end
            
            % biographs are handle objects, so I need to create a new one
            % each time to make copies
            bg(si) = biograph(C,[nodesnull;reactionsnull]);
            
            isreactionnode = false(nxunull+nrnull,1);
            isreactionnode(nxunull+1:end) = true;
            
            % Set node labels
            nodeids = get(bg(si).Nodes,'ID'); % Get node IDs
            [~,statestonodes] = ismember([nodesnull;reactionsnull],nodeids); % Map states & reactions list to nodes list, in case they are in different orders
            reactionnodecolor = [1 1 1];
            reactionnodeshape = 'ellipse';
            for ni = 1:length(nodeids)
                set(bg(si).Nodes(ni),'Label',nodelabels{statestonodes(ni)});
                if isreactionnode(statestonodes(ni))
                    ri = rindices(ni);
                    reactionnodedescription = rnames_{ri};
                    set(bg(si).Nodes(ni),'Color',reactionnodecolor,'Shape',reactionnodeshape,'Description',reactionnodedescription)
                end
            end
            
            set(bg(si).Edges,'Weight',0); % Set weights to 0 to start
            set(bg(si).Edges,'LineWidth',0.1); % Set line widths to 0 to start

            edges = get(bg(si),'Edges');
            nEdges = length(edges);
            edgesourcesandsinks = regexp(get(edges,'ID'),'(?<source>.*)\s->\s(?<sink>.*)','names');
            edgesourcesandsinks = vertcat(edgesourcesandsinks{:});
            [~,edgesources] = ismember({edgesourcesandsinks.source}',[nodesnull(:);reactionsnull(:)]);
            [~,edgesinks] = ismember({edgesourcesandsinks.sink}',[nodesnull(:);reactionsnull(:)]);
            for ei = 1:nEdges
                thisedge = edges(ei);
                thissource = edgesources(ei);
                thissink = edgesinks(ei);
                set(thisedge,'Weight',round(weights_C(thissource,thissink),2,'significant'));
                set(thisedge,'LineWidth',linewidthscale*linewidths_C(thissource,thissink));
                set(thisedge,'LineColor',linecolors_C{thissource,thissink});
            end
            
            
            % Plot network diagram
            bgview = bg(si).view;
            %set(bgview,'ShowWeights','on','ShowTextInNodes','label')
            set(bgview,'ShowTextInNodes','label')
        
        end
        
    end

    function rscaled = plotReactions(t,m,con,sim)
        
        if nargin < 4
            sim = getSim(t,m,con);
        end
        rscaled = getReactionRates(t,sim,m,con);
        
        % Plot the reaction rates
        figure;
        %plotfun = @(rscaled,p) plotbar(rscaled(sum(subset.Sr,1) ~= 0,:)',size(m,1),size(con,2),p);
        %hb = cellfun(plotfun,rscaled,reshape(num2cell(1:numel(con)),size(m,1),size(con,2)),'UniformOutput',false);
        rxnstoplot = any(abs(rscaled) >= repmat(0.01*max(abs(rscaled),[],1),size(rscaled,1),1), 2);
        
        % If there are multiple model/experiment combinations to plot...
        if size(rscaled,2) > 1
            hb = bar(rscaled(rxnstoplot,:)');
            colormap(jet(sum(rxnstoplot)));
            legend(regexprep(subset.rnames(rxnstoplot),'_','\\_'),'Location','EastOutside')
        % If there is only one model/experiment combination...
        else
            hb = bar(rscaled(rxnstoplot));
            set(gca,'XTickLabel',subset.rnames(rxnstoplot),'XTickLabelRotation',45);
        end
        
        hf = gcf; pos = hf.Position;
        hf.Position = [pos(1:2)      1086         474];
        ylabel('Reaction''s contribution to d[species]/dt')
    end

    function [rscaled,r] = getReactionRates(t,sim,m,con)
        
        assert(isscalar(t),'t must be a scalar')
        assert(isempty(con) || (isscalar(con) && isstruct(con)), 'con must be a scalar experiment struct')
        t = repmat(t,size(m));
        
        % Update m to include parameter resets in con
        if ~isempty(con)
            for ii = 1:numel(m)
                if isfield(con(ii), 'k') && ~isempty(con(ii).k)
                    m(ii).k(con(ii).k(:,1)) = con(ii).k(:,2);
                    m(ii) = m(ii).Update(m(ii).k,m(ii).x0,m(ii).q);
                end
            end
        end
        
        % Set up functions to convert state and input values to reaction
        % rates to stoichiometry-and-C1-scaled reaction rates
        if isa(sim.x, 'function_handle')
            xfun = @(t)sim.x(t);
        else
            xfun = @(t)sim.x(:,t == sim.t);
        end
        if isa(sim.u, 'function_handle')
            ufun = @(t)sim.u(t);
        else
            ufun = @(t)sim.u(:,t == sim.t);
        end
        getfullr =  @(m,t) m.r(t,xfun(t),ufun(t));
        uniqueSj = unique(subset.Sj);
        getreducedr = @(r) r(uniqueSj);
        scaler = @(r) sum(bsxfun(@times,bsxfun(@times,r,subset.Sr'),subset.N1(subset.N1 ~= 0)),2);
        getr = @(m,t) getreducedr(getfullr(m,t));
        getrscaled = @(m,t) scaler(getreducedr(getfullr(m,t)));
        
        % Perform the conversions
        r = arrayfun(getr,m,t,'UniformOutput',false);
        rscaled = arrayfun(getrscaled,m,t,'UniformOutput',false);
        r = [r{:}];
        rscaled = [rscaled{:}];
        
    end

    function x = getStateValues(t,sim,isexpanded)
        
        if nargin < 3
            isexpanded = false;
        end
        
        
        nsims = numel(sim);
        
        if isexpanded
            xstoget = subset.C1expanded ~= 0;
        else
            xstoget = subset.C1 ~= 0;
        end
        
        nxr = sum(xstoget);
        
        x = NaN(nxr,nsims);
        for si = 1:numel(sim)
            thissim = sim(si);
            if isa(thissim.x, 'function_handle')
                xtemp = thissim.x(t);
            else
                xtemp = thissim.x(:,t == thissim.t);
            end
            x(:,si) = xtemp(xstoget,:);
        end
        
    end

    function u = getInputValues(t,sim,isexpanded)
        
        if nargin < 3
            isexpanded = false;
        end
        
        nsims = numel(sim);
        
        if isexpanded
            ustoget = subset.U1expanded ~= 0;
        else
            ustoget = subset.U1 ~= 0;
        end
        
        nur = sum(ustoget);
        
        u = NaN(nur,nsims);
        for si = 1:numel(sim)
            thissim = sim(si);
            if isa(thissim.u, 'function_handle')
                utemp = thissim.u(t);
            else
                utemp = thissim.u(:,t == thissim.t);
            end
            u(:,si) = utemp(ustoget,:);
        end
        
    end

    function sim = getSim(t,m,con)
        
        % Check inputs
        assert(isscalar(t),'t must be a scalar')
        
        % Ensure m, con, and t are the same size, and error if they are
        % incompatible
        [m,con,t] = sxrepmat(m,con,t);
        sim = arrayfun(@SimulateSystem,m,con,t);
        
%         % Simulate or reassign simulation to sim
%         if strncmp(con.Type,'Experiment',10)
% 
%         elseif strncmp(con.Type,'Simulation',10)
%             sim = con;
%         else
%             error('Unrecognized input for con.')
%         end
        
    end


%     function hb = plotbar(data,nm,ncon,p)
%         subplot(nm,ncon,p)
%         hb = bar(data);
%         [icon,~] = ind2sub([nm ncon],p);
% %         if icon == ncon
% %             set(gca,'XTick',1:length(data));
% %             set(gca,'XTickLabel',subset.rnames);
% %             set(gca,'XTickLabelRotation',90)
% %         end
%     end

end

