function [p,infoTog] = createInfoPanel(f,ax,sim)

p = uipanel(f,'title','Information','units','normalized',...
    'position',[0.75 0.55 0.2 0.4]);

%% create text boxes inside panel
n = 2;
tbox = zeros(1,n);
for k = 1:n
    tbox(k) = uicontrol(p,'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'position',[0.1 0.85-(k-1)*0.15 0.8 0.1]);    
end

%% add labels
addTextLabels(tbox,sim);

infoTog = @(popup)(addTextLabels(tbox,sim,popup));

%% create buttons
n = 2;
infoButtons = zeros(1,n);

infoButtons(1) = uicontrol(p,'Style','pushbutton','units','normalized',...
    'HorizontalAlignment','Left','string','Show Nodes',...
    'position',[0.05 0.45 0.9 0.2],'callback',@(~,~)sim.showMap(ax));

infoButtons(2) = uicontrol(p,'Style','pushbutton','units','normalized',...
    'HorizontalAlignment','Left','string','Show Connections',...
    'position',[0.05 0.25 0.9 0.2],...
    'callback',@(~,~)sim.showConnections(ax));

% show the connections and nodes by default
sim.showMap(ax);
sim.showConnections(ax);

end

function addTextLabels(tbox,sim,popup)

data = guidata(tbox(1));

if(data.nodeMode)

    num = sim.getNumNodes;
    setTBoxTxt(tbox(1),'# of Nodes',num);
    setTBoxTxt(tbox(2),'# of Targets',0);
    
else
    
    idx = get(popup,'value');
    loc = sim.getNodeLoc(idx);
    
    str = sprintf('Node %2d: (%.1f,%.1f)',idx,loc(1),loc(2));        
    set(tbox(1),'string',str);
    
    num = sim.getNumTracks(idx);
    setTBoxTxt(tbox(2),'# of Tracks',num);
    
end

end

function setTBoxTxt(h,txt,num)

txt = sprintf('%12s: %3d',txt,num);
set(h,'string',txt);

end