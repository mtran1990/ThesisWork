function createMainFig(figNum,sim)

f = figure(figNum);
clf(f);

%% initialize figure
xwidth = 600;
ywidth = 400;
set(f,'Visible','on','Position',[360,500,xwidth,ywidth]);

data = struct('nodeMode',true,'lastPopVal',1);
guidata(f,data);

%% create plot axes
ax = axes('units','normalized','Position',[0.1 0.1 0.6 0.8]);

sim.initMap(ax);
%% create info panel
[p,infoTog] = createInfoPanel(f,ax,sim);

%% create option panel
p2 = uipanel(f,'title','Options','units','normalized',...
    'position',[0.75 0.05 0.2 0.45]);

% create popup menu
str = getPopUpString(p2,sim);
nodePop = uicontrol(p2,'Style','popupmenu','units','normalized',...
    'HorizontalAlignment','Left','string',str,...
    'position',[0.1 0.9 0.8 0.1]);

% create buttons
nodeButtons = zeros(1,2);

nodeButtons(1) = uicontrol(p2,'Style','pushbutton','units','normalized',...
    'HorizontalAlignment','Left','string','Show Tracks',...
    'position',[0.05 0.6 0.9 0.2],...
    'callback',@(~,~)showTracks(sim,nodePop,ax));

nodeButtons(2) = uicontrol(p2,'Style','pushbutton','units','normalized',...
    'HorizontalAlignment','Left','string','Show Point View',...
    'position',[0.05 0.4 0.9 0.2],...
    'callback',@(~,~)toggleMode(nodePop,sim,infoTog));

end

function str = getPopUpString(nodePop,sim)

data = guidata(nodePop);

if(data.nodeMode)
    N = sim.getNumNodes;

    loc = sim.getNodeLoc;

    str = cell(1,N);

    for k = 1:N

        x = loc(1,k);
        y = loc(2,k);

        str{k} = sprintf('Node %2d: (%.1f,%.1f)',k,x,y);

    end
else
    
    idx = get(nodePop,'Value');
    N = sim.getNumTracks(idx);
    
    str = cell(1,N);
    
    for k = 1:N
        str{k} = sprintf('Track %2d',k);
    end
    
end
    
end

function togglePopup(nodePop,sim)

str = getPopUpString(nodePop,sim);
set(nodePop,'value',1);
set(nodePop,'string',str);

end

function toggleMode(nodePop,sim,infoTog)

sim.hidePlots;

data = guidata(nodePop);

if(data.nodeMode)
    data.nodeMode = false;
else
    data.nodeMode = true;
end
data.lastPopVal = get(nodePop,'value');

guidata(nodePop,data);

infoTog(nodePop);
togglePopup(nodePop,sim);

end

function showTracks(sim,nodePop,ax)

data = guidata(nodePop);

if(data.nodeMode)
    
    node = get(nodePop,'value');
    sim.showTracks(ax,node);
    
else
    
    node = data.lastPopVal;
    track = get(nodePop,'value');
    sim.showTracks(ax,node,track,0);
    
end

end
