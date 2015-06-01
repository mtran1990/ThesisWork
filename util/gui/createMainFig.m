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
    'position',[0.05 0.4 0.9 0.2]);

% create panel to show scroll through point-view
timeBox(1) = uicontrol(p2,'Style','text','string','Current Time:',...
    'units','normalized','HorizontalAlignment','Left','visible','off',...
    'position',[0.05 0.3 0.9 0.1]);

timeBox(2) = uicontrol(p2,'Style','text','string','Time Range (s):',...
    'units','normalized','HorizontalAlignment','Left','visible','off',...
    'position',[0.05 0.2 0.9 0.1]);

timeBox(3) = uicontrol(p2,'Style','text','string','[0 0]',...
    'units','normalized','HorizontalAlignment','Left','visible','off',...
    'position',[0.05 0.1 0.9 0.1]);

sBar = uicontrol(p2,'Style','slider','units','normalized',...
    'position',[0.05 0 0.9 0.1],'visible','off');

scrollHandles = [timeBox,sBar];

% need to set callback for point-view button
set(nodeButtons(2),...
    'callback',@(nodeBut,~)toggleMode(nodeBut,nodePop,scrollHandles,sim,infoTog));

% need to set callback for popup menu
set(nodePop,...
    'callback',@(~,~)popUpCallback(nodePop,scrollHandles));

% need to set callback for slider
set(sBar,...
    'callback',@(h,eData)sliderCallback(h,ax,nodePop,scrollHandles,sim,eData));

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

function popUpCallback(nodePop,scrollHandles)

data = guidata(nodePop);

if(~data.nodeMode)
    toggleScroll(nodePop,scrollHandles);    
end

end

function toggleScroll(nodePop,scrollHandles)

data = guidata(nodePop);

if(data.nodeMode)
    set(scrollHandles,'visible','off');
else
    % get the slider callback function and run it to update the map
    initPlot = get(scrollHandles(end),'callback');
    initPlot(scrollHandles(end),true);
    
    
    
    set(scrollHandles,'visible','on');
end


end

function setTimeBoxes(scrollHandles,minT,maxT)

t = get(scrollHandles(end),'value');
str = sprintf('Current Time: %.2fs',t);
set(scrollHandles(1),'string',str);
str = sprintf('[%.2f %.2f]',minT,maxT);
set(scrollHandles(3),'string',str);

end

function sliderCallback(h,ax,nodePop,scrollHandles,sim,bool)

data = guidata(h);

t = get(h,'value');
node = data.lastPopVal;
track = get(nodePop,'value');

info = sim.nodeList(node).getTrackInfo(track);

dt = sim.tParams.dt;
maxT = max(info.t);
minT = min(info.t);
sliderStep(2) = 0.2;
sliderStep(1) = dt/(maxT-minT);

set(scrollHandles(end),'min',minT,'max',maxT,...
    'sliderStep',sliderStep);

if(bool)
    set(scrollHandles(end),'value',minT);
end

setTimeBoxes(scrollHandles,minT,maxT);

sim.showTracks(ax,node,track,t);

end

function toggleMode(nodeBut,nodePop,scrollHandles,sim,infoTog)

sim.hidePlots;

data = guidata(nodePop);

if(data.nodeMode)
    data.nodeMode = false;
    str = 'Show Node View';
else
    data.nodeMode = true;
    str = 'Show Point View';
end
set(nodeBut,'string',str);
data.lastPopVal = get(nodePop,'value');

guidata(nodePop,data);

infoTog(nodePop);
togglePopup(nodePop,sim);
toggleScroll(nodePop,scrollHandles);

end

function showTracks(sim,nodePop,ax)

data = guidata(nodePop);

if(data.nodeMode)
    
    node = get(nodePop,'value');
    sim.showTracks(ax,node);
    
else
    
%     node = data.lastPopVal;
%     track = get(nodePop,'value');
%     sim.showTracks(ax,node,track,0);
    
end

end