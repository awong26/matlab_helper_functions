function varargout = markdataGUI(varargin)
% markdataGUI    visualization program for movement data
%      markdataGUI(dat) creates a user interface for visualizing
%      and/or marking the data in dat from a single trial (movement). dat
%      can have one of two formats; it may be an array of structures where
%      each structure contains the x, y, and z coordinates of the movement
%      in the format m(i).x, m(i).y, and m(i).z where i is the index of the
%      marker (up to 8 trackers can be provided). alternatively, dat may be
%      a 3D array with dimensions [samples, axes, markers]; this should
%      therefore be a [: by 3 by 8] array. If a structure or array with
%      fewer markers is provided, the remaining dimensions are assumed to
%      be NaN; any additional markers beyond 8 will be ignored.
%
%      markdataGUI assumes that the markers are as follows:
%          1: thumb*
%          2: index finger*
%          3: middle finger
%          4: back of the hand
%          5: wrist*
%          6: elbow*
%          7: shoulder*
%          8: opposite shoulder
%      where the markers indicated with an asterisk are actually displayed
%      in the GUI. Note, markers 7 and 8 are used to align/center the other
%      markers as a reference point, so using these indices for other data
%      may cause confusion in the visualization.
%
%      [inds] = markdataGUI(dat) provides as output the array of indices
%      that were marked during the call to markdataGUI. This output is
%      only populated when the GUI window is closed.
%
%      Clicking on the slider or on any of the axes on the right side of
%      the figure will position a dashed line on the right-hand axes
%      indicating the current index, and will update the 3D plot. When an
%      index of interest has been identified, it can be saved by using the
%      "Mark" button. All currently saved indices will be denoted by a
%      vertical red line on the right-hand axes. To remove a mark, choose
%      a position close to the the mark to be removed and click the
%      "Unmark" button, which will remove the closest marked index. 
%
%      The "play" and "stop" features will automatically step through the
%      data set like a movie. Note that in some cases (i.e., when the
%      processor or Matlab are very slow) this can cause unintended
%      blocking, preventing additional updates to be made. In this case, it
%      may become necessary to force-quit out of the GUI.



% MARKDATAGUI MATLAB code for markdataGUI.fig
%      MARKDATAGUI, by itself, creates a new MARKDATAGUI or raises the existing
%      singleton*.
%
%      H = MARKDATAGUI returns the handle to a new MARKDATAGUI or the handle to
%      the existing singleton*.
%
%      MARKDATAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARKDATAGUI.M with the given input arguments.
%
%      MARKDATAGUI('Property','Value',...) creates a new MARKDATAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before markdataGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to markdataGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help markdataGUI

% Last Modified by GUIDE v2.5 02-Feb-2018 11:05:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @markdataGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @markdataGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before markdataGUI is made visible.
function markdataGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to markdataGUI (see VARARGIN)

%create the data plot
m = varargin{1};

inds = [];
nmarkers = 6;

i = 2;
while i <= length(varargin)
    if ~ischar(varargin{i}) && ~isscalar(varargin{i}) %catch cases where invalid input into switch
        continue;
    end
    switch(varargin{i})
        case 'rot'
            if isstruct(m)
                for a = 1:length(m)
                    m(a).x = m(a).rotx;
                    m(a).y = m(a).roty;
                    m(a).z = m(a).rotz;
                end
            end
            i = i+1;
        case 'title'
            ttlstr = varargin{i+1};
            set(gcf,'Name',sprintf('markdataGUI: %s',ttlstr));
            i = i+2;
        case 'mark'
            inds = varargin{i+1};
            i = i+2;
        case 'nmarkers'
            nmarkers = varargin{i+1};
            i = i+2;
        otherwise
            i = i+1;
    end
end
            

if ~isstruct(m)
    for a = 1:size(m,3)
        n(a).x = m(:,1,a);
        n(a).y = m(:,2,a);
        n(a).z = m(:,3,a);
    end
    m = n;
end

for a = length(m)+1:7
    m(a).x = zeros(size(m(1).x));
    m(a).y = zeros(size(m(1).y));
    m(a).z = zeros(size(m(1).z));
end

if length(m) < 8
    m(8).x = m(7).x-30;
    m(8).y = m(7).y;
    m(8).z = m(7).z;
end


% for a = 1:size(m,3)
%     m(a).y = -m(a).y;
% end

handles.torso = [(m(7).x+m(8).x)/2 (m(7).y+m(8).y)/2 (m(7).z+m(8).z)/2];
%handles.torso = [zeros(size(m(7).x)) zeros(size(m(7).y)) zeros(size(m(7).z))];

handles.joint.refshoulder = [m(8).x m(8).y m(8).z] - handles.torso;
handles.joint.shoulder = [m(7).x m(7).y m(7).z] - handles.torso;
handles.joint.elbow =    [m(6).x m(6).y m(6).z] - handles.torso;
handles.joint.wrist =    [m(5).x m(5).y m(5).z] - handles.torso;
handles.joint.hand =     [m(4).x m(4).y m(4).z] - handles.torso;
handles.joint.middlefinger = [m(3).x m(3).y m(3).z] - handles.torso;
handles.joint.indexfinger = [m(2).x m(2).y m(2).z] - handles.torso;
handles.joint.thumbfinger = [m(1).x m(1).y m(1).z] - handles.torso;
handles.nmarkers = nmarkers;

xmin = min([handles.joint.thumbfinger(:,1); handles.joint.indexfinger(:,1); handles.joint.wrist(:,1); handles.joint.elbow(:,1); handles.joint.shoulder(:,1); handles.joint.refshoulder(:,1)]);
xmax = max([handles.joint.thumbfinger(:,1); handles.joint.indexfinger(:,1); handles.joint.wrist(:,1); handles.joint.elbow(:,1); handles.joint.shoulder(:,1); handles.joint.refshoulder(:,1)]);
ymin = min([handles.joint.thumbfinger(:,2); handles.joint.indexfinger(:,2); handles.joint.wrist(:,2); handles.joint.elbow(:,2); handles.joint.shoulder(:,2); handles.joint.refshoulder(:,2)]);
ymax = max([handles.joint.thumbfinger(:,2); handles.joint.indexfinger(:,2); handles.joint.wrist(:,2); handles.joint.elbow(:,2); handles.joint.shoulder(:,2); handles.joint.refshoulder(:,2)]);
zmin = min([handles.joint.thumbfinger(:,3); handles.joint.indexfinger(:,3); handles.joint.wrist(:,3); handles.joint.elbow(:,3); handles.joint.shoulder(:,3); handles.joint.refshoulder(:,3)]);
zmax = max([handles.joint.thumbfinger(:,3); handles.joint.indexfinger(:,3); handles.joint.wrist(:,3); handles.joint.elbow(:,3); handles.joint.shoulder(:,3); handles.joint.refshoulder(:,3)]);

handles.plot.xmin = floor(xmin/5)*5;
handles.plot.xmax = ceil(xmax/5)*5;
handles.plot.ymin = floor(ymin/5)*5;
handles.plot.ymax = ceil(ymax/5)*5;
handles.plot.zmin = floor(zmin/5)*5;
handles.plot.zmax = ceil(zmax/5)*5;
if( handles.plot.zmax < (max([handles.joint.shoulder(:,3); handles.joint.refshoulder(:,3)]) + abs(handles.joint.refshoulder(1,1)-handles.joint.shoulder(1,1))/4) )
    handles.plot.zmax = (max([handles.joint.shoulder(:,3); handles.joint.refshoulder(:,3)]) + abs(handles.joint.refshoulder(1,1)-handles.joint.shoulder(1,1))/2*1.5);
    %handles.plot.zmax = ceil(zmax/5)*5;
end

% handles.plot.xmin = -80;
% handles.plot.xmax = 80;
% handles.plot.ymin = -70;
% handles.plot.ymax = 70;
% handles.plot.zmin = -60;
% handles.plot.zmax = 60;

%handles.viewaz = -25;
%handles.viewel = 25;
handles.viewaz = 195;
handles.viewel = 20;

set(handles.slider1,'Min',1);
set(handles.slider1,'Max',size(handles.joint.shoulder,1));
set(handles.slider1,'SliderStep',[1/size(handles.joint.shoulder,1) , 10/size(handles.joint.shoulder,1)]);

%handles.c = 1;
c = 1;
setappdata(handles.figure1,'c',c);

%handles.inds = [];
%inds = [];
setappdata(handles.figure1,'inds',inds);

playtimer = timer('Name','DataPlayer','Period',0.25,'StartDelay',0.1,'ExecutionMode','fixedRate','TimerFcn',{@timerCallbackFcn,handles});
stop(playtimer);
setappdata(handles.figure1,'playtimer',playtimer);

updateplot(handles);

% Choose default command line output for markdataGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes markdataGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% function [] = timerCallback(~,~,guiHandle)
% 
% if ~isempty(guiHandle)
%     
%     % get the handles for the GUI/figure
%     handles = guidata(guiHandle);
% 
% end


% --- Outputs from this function are returned to the command line.
function varargout = markdataGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
inds = getappdata(handles.figure1,'inds');
varargout{1} = inds;

delete(handles.figure1);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton4,'enable','off');
set(handles.pushbutton5,'enable','off');

playtimer = getappdata(handles.figure1,'playtimer');
start(playtimer);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
playtimer = getappdata(handles.figure1,'playtimer');
stop(playtimer);

set(handles.pushbutton4,'enable','on');
set(handles.pushbutton5,'enable','on');



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
inds = getappdata(handles.figure1,'inds');
c = getappdata(handles.figure1,'c');

inds = [inds; c];
inds = unique(inds); %get rid of duplicate marks
inds = sort(inds);   %resort the indices in ascending order
setappdata(handles.figure1,'inds',inds);

updateplot(handles);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

inds = getappdata(handles.figure1,'inds');
c = getappdata(handles.figure1,'c');
if ~isempty(inds)
    diffinds = inds-c;
    
    [~,idiff] = min(abs(diffinds));
    inds(idiff) = [];
end

setappdata(handles.figure1,'inds',inds);

updateplot(handles);



% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

playtimer = getappdata(handles.figure1,'playtimer');
stop(playtimer);

set(handles.pushbutton4,'enable','on');
set(handles.pushbutton5,'enable','on');

c = get(handles.slider1,'Value');
c = round(c);
%cmin = get(handles.slider1,'Min');
%cmax = get(handles.slider1,'Max');
%c = round((c-cmin)/(cmax-cmin) * size(handles.joint.shoulder,1));

if c < 1
    c = 1;
end
if c > size(handles.joint.shoulder,1)
    c = size(handles.joint.shoulder,1);
end

setappdata(handles.figure1,'c',c);

updateplot(handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Min',1);
set(hObject,'Max',2000);
set(hObject,'Value',1);
set(hObject,'SliderStep',[1/2000 , 10/2000]);




% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doClickCallback(handles,handles.axes2)


% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doClickCallback(handles,handles.axes3)


% --- Executes on mouse press over axes background.
function axes4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doClickCallback(handles,handles.axes4)


% --- Executes on mouse press over axes background.
function axes5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doClickCallback(handles,handles.axes5)


% --- Executes on mouse press over axes background.
function axes6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doClickCallback(handles,handles.axes6)




% --- create new function to update plot
function updateplot(handles)

c = getappdata(handles.figure1,'c');
inds = getappdata(handles.figure1,'inds');

if handles.nmarkers == 6
    %do simplified arm without wrist/hand
    [xSE,ySE,zSE] = cylinder2P(.45,8,handles.joint.shoulder(c,:),handles.joint.elbow(c,:));
    [xEW,yEW,zEW] = cylinder2P(.375,8,handles.joint.elbow(c,:),handles.joint.wrist(c,:));
    [xWF,yWF,zWF] = cylinder2P([.25 .25 .25 .1 .1 .1],8,handles.joint.wrist(c,:),handles.joint.indexfinger(c,:));
    [xWT,yWT,zWT] = cylinder2P([.25 .25 .25 .1 .1 .1],8,handles.joint.wrist(c,:),handles.joint.thumbfinger(c,:));
elseif handles.nmarkers == 8
    [xSE,ySE,zSE] = cylinder2P(.45,8,handles.joint.shoulder(c,:),handles.joint.elbow(c,:));
    [xEW,yEW,zEW] = cylinder2P(.375,8,handles.joint.elbow(c,:),handles.joint.wrist(c,:));
    [xWH,yWH,zWH] = cylinder2P(.375,8,handles.joint.wrist(c,:),handles.joint.hand(c,:));
    [xHM,yHM,zHM] = cylinder2P([.25 .25 .25 .1 .1 .1],8,handles.joint.hand(c,:),handles.joint.middlefinger(c,:));
    [xWF,yWF,zWF] = cylinder2P([.25 .25 .25 .1 .1 .1],8,handles.joint.hand(c,:),handles.joint.indexfinger(c,:));
    [xWT,yWT,zWT] = cylinder2P([.25 .25 .25 .1 .1 .1],8,handles.joint.hand(c,:),handles.joint.thumbfinger(c,:));
end

%rectangle representing the body, assuming a vertical plane
torsolength = abs(handles.joint.shoulder(1,1)-handles.joint.refshoulder(1,1))*1.5;
xtorso = [handles.joint.shoulder(c,1); handles.joint.refshoulder(c,1); handles.joint.refshoulder(c,1); handles.joint.shoulder(c,1)];
ytorso = [handles.joint.shoulder(c,2); handles.joint.refshoulder(c,2); handles.joint.refshoulder(c,2); handles.joint.shoulder(c,2)];
ztorso = [handles.joint.shoulder(c,3); handles.joint.refshoulder(c,3); handles.joint.refshoulder(c,3)-torsolength; handles.joint.shoulder(c,3)-torsolength];

%headrad = 9;
headrad = abs(handles.joint.refshoulder(1,1)-handles.joint.shoulder(1,1))/4;
headc = handles.torso(c,:) + [0 0 headrad];
theta = [0:pi/10:2*pi];
headx = headrad*cos(theta);
heady = zeros(size(headx));
headz = headrad*sin(theta)+headrad;

axes(handles.axes1)
cla(handles.axes1);
h = patch('XData',xtorso,'YData',ytorso,'ZData',ztorso);
set(h,'FaceColor',[.5 .5 .5]);
hold on;
h = patch('XData',headx,'YData',heady,'ZData',headz);
set(h,'FaceColor',[.5 .5 .5]);
h = surf(xSE, ySE, zSE);
set(h,'FaceColor',[1 0 0]);
h = surf(xEW, yEW, zEW);
set(h,'FaceColor',[1 1 0]);
if handles.nmarkers == 8
    h = surf(xWH, yWH, zWH);
    set(h,'FaceColor',[1 1 0]);
    h = surf(xHM, yHM, zHM);
    set(h,'FaceColor',[0 0 1]);
end
h = surf(xWF, yWF, zWF);
set(h,'FaceColor',[0 0 1]);
h = surf(xWT, yWT, zWT);
set(h,'FaceColor',[0 0 1]);
hold off;
grid on;
set(gca,'xlim',[handles.plot.xmin handles.plot.xmax],'ylim',[handles.plot.ymin handles.plot.ymax],'zlim',[handles.plot.zmin handles.plot.zmax]);
view(handles.viewaz,handles.viewel);
xlabel('x');
ylabel('y');
zlabel('z');

axes(handles.axes2)
cla(handles.axes2);
set(handles.axes2,'NextPlot','add');
h = plot(handles.joint.shoulder,'-');
set(h,'HitTest','off')
hold on
h = plot([c c],get(gca,'ylim'),'k:');
set(h,'HitTest','off')
for a = 1:length(inds)
    h = plot([inds(a) inds(a)],get(gca,'ylim'),'r-');
    set(h,'HitTest','off')
end
hold off;
ylabel('shoulder')


axes(handles.axes3)
cla(handles.axes3);
set(handles.axes3,'NextPlot','add');
h = plot(handles.joint.elbow,'-');
set(h,'HitTest','off');
hold on
h = plot([c c],get(gca,'ylim'),'k:');
set(h,'HitTest','off');
for a = 1:length(inds)
    h = plot([inds(a) inds(a)],get(gca,'ylim'),'r-');
    set(h,'HitTest','off')
end
hold off;
ylabel('elbow')

axes(handles.axes4)
cla(handles.axes4);
set(handles.axes4,'NextPlot','add');
h = plot(handles.joint.wrist,'-');
set(h,'HitTest','off');
hold on
h = plot([c c],get(gca,'ylim'),'k:');
set(h,'HitTest','off');
for a = 1:length(inds)
    h = plot([inds(a) inds(a)],get(gca,'ylim'),'r-');
    set(h,'HitTest','off')
end
hold off;
ylabel('wrist')

axes(handles.axes5)
cla(handles.axes5);
set(handles.axes5,'NextPlot','add');
h = plot(handles.joint.indexfinger,'-');
set(h,'HitTest','off');
hold on
h = plot([c c],get(gca,'ylim'),'k:');
set(h,'HitTest','off');
for a = 1:length(inds)
    h = plot([inds(a) inds(a)],get(gca,'ylim'),'r-');
    set(h,'HitTest','off')
end
hold off;
ylabel('index')

axes(handles.axes6)
cla(handles.axes6);
set(handles.axes6,'NextPlot','add');
h = plot(handles.joint.thumbfinger,'-');
set(h,'HitTest','off');
hold on
h = plot([c c],get(gca,'ylim'),'k:');
set(h,'HitTest','off');
for a = 1:length(inds)
    h = plot([inds(a) inds(a)],get(gca,'ylim'),'r-');
    set(h,'HitTest','off')
end
hold off;
ylabel('thumb')



function doClickCallback(handles,haxes)

playtimer = getappdata(handles.figure1,'playtimer');
stop(playtimer);

set(handles.pushbutton4,'enable','on');
set(handles.pushbutton5,'enable','on');

coordinates = get(haxes,'CurrentPoint');
coords = coordinates(1,1:2);
c = round(coords(1));

%cmin = get(handles.slider1,'Min');
%cmax = get(handles.slider1,'Max');
%set(handles.slider1,'Value',(c/size(handles.joint.shoulder,1)*(cmax-cmin)+cmin));
set(handles.slider1,'Value',c);

setappdata(handles.figure1,'c',c);

updateplot(handles)


function timerCallbackFcn(obj, event, handles)
c = getappdata(handles.figure1,'c');
c = c+10;

if c > size(handles.joint.shoulder,1)
    c = 1;
    playtimer = getappdata(handles.figure1,'playtimer');
    stop(playtimer);
    
    set(handles.pushbutton4,'enable','on');
    set(handles.pushbutton5,'enable','on');
end

setappdata(handles.figure1,'c',c);

%cmin = get(handles.slider1,'Min');
%cmax = get(handles.slider1,'Max');
%set(handles.slider1,'Value',(c/size(handles.joint.shoulder,1)*(cmax-cmin)+cmin));
set(handles.slider1,'Value',c);

updateplot(handles);




% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

playtimer = getappdata(handles.figure1,'playtimer');
stop(playtimer);

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

% Hint: delete(hObject) closes the figure
%delete(hObject);
