function createMainFig(figNum,sim)

f = figure(figNum);
clf(f);

xwidth = 600;
ywidth = 400;
set(f,'Visible','on','Position',[360,500,xwidth,ywidth]);

%% create plot axes
ax = axes('units','normalized','Position',[0.1 0.1 0.6 0.8]);

sim.initMap(ax);
%% create info panel
p = createInfoPanel(f,ax,sim);

%% create option panel
p2 = uipanel(f,'title','Options','units','normalized',...
    'position',[0.75 0.05 0.2 0.45]);

% create popup menu
str = getPopUpString(sim);
nodePop = uicontrol(p2,'Style','popupmenu','units','normalized',...
    'HorizontalAlignment','Left','string',str,...
    'position',[0.1 0.9 0.8 0.1]);

% create buttons
nodeButtons = zeros(1,2);

nodeButtons(1) = uicontrol(p2,'Style','pushbutton','units','normalized',...
    'HorizontalAlignment','Left','string','Show Tracks',...
    'position',[0.05 0.6 0.9 0.2],...
    'callback',@(~,~)sim.showTracks(nodePop,ax));

nodeButtons(2) = uicontrol(p2,'Style','pushbutton','units','normalized',...
    'HorizontalAlignment','Left','string','Show Point View',...
    'position',[0.05 0.4 0.9 0.2],...
    'callback',@(~,~)sim.hidePlots);

end

function str = getPopUpString(sim)

N = sim.getNumNodes;

loc = sim.getNodeLoc;

str = cell(1,N);

for k = 1:N
    
    x = loc(1,k);
    y = loc(2,k);
    
    str{k} = sprintf('Node %2d: (%.1f,%.1f)',k,x,y);
    
end

end


