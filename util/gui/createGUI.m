function createGUI(sim)

f = figure(1);
clf(f);

xwidth = 600;
ywidth = 400;
set(f,'Visible','on','Position',[360,500,xwidth,ywidth]);

%% create plot axes
ax = axes('units','normalized','Position',[0.1 0.1 0.6 0.8]);

%% create info panel
p = uipanel(f,'title','Information','units','normalized',...
    'position',[0.75 0.55 0.2 0.4]);

% create text boxes inside panel
n = 2;
tbox = zeros(1,n);
for k = 1:n
    tbox(k) = uicontrol(p,'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'position',[0.1 0.85-(k-1)*0.15 0.8 0.1]);    
end

% add labels
setTBoxTxt(tbox(1),'# of Nodes',0);
setTBoxTxt(tbox(2),'# of Targets',0);

% create buttons
n = 2;
infoButtons = zeros(1,n);

infoButtons(1) = uicontrol(p,'Style','pushbutton','units','normalized',...
    'HorizontalAlignment','Left','string','Show Nodes',...
    'position',[0.05 0.45 0.9 0.2]);

infoButtons(2) = uicontrol(p,'Style','pushbutton','units','normalized',...
    'HorizontalAlignment','Left','string','Show Connections',...
    'position',[0.05 0.25 0.9 0.2]);

%% create option panel
p2 = uipanel(f,'title','Options','units','normalized',...
    'position',[0.75 0.05 0.2 0.45]);

% create popup menu
nodeBox = uicontrol(p2,'Style','popupmenu','units','normalized',...
    'HorizontalAlignment','Left','string','Node 1: (x,y)',...
    'position',[0.1 0.9 0.8 0.1]);

% create buttons
nodeButtons = uicontrol(p2,'Style','pushbutton','units','normalized',...
    'HorizontalAlignment','Left','string','Show Tracks',...
    'position',[0.05 0.6 0.9 0.2]);

end

function setTBoxTxt(h,txt,num)

txt = sprintf('%12s: %3d',txt,num);
set(h,'string',txt);

end

function showNodes(h,sim)

    

end