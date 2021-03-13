function minecraft
%% Initialization
clear; clc; close all; opengl hardware; import java.awt.Robot;
% global variables
UPDATE_LIMIT = 0;       % minimum time before game state is updated (performance)
DRAW_LIMIT = 0;         % minimum time before frame is redrawn (performance)
CHUNK_SIZE = 17;        % chunk width in block lengths (requires int)
SEED = tic;             % unique game seed
RANGE = 100;            % camera range in block lengths
FOV = 45;               % camera field of view in degrees
MOVE_SPEED = 15;        % move speed multiplier
ROTATION_SPEED = 25e-4; % rotation speed multiplier
KEY_SHIFT = 'shift';    % shift key id
KEY_W = 'w';            % w key id
KEY_A = 'a';            % a key id
KEY_S = 's';            % s key id
KEY_D = 'd';            % d key id
KEY_SPACE = 'space';    % space key id
KEY_ESCAPE = 'escape';	% escape key id
KEY_ALT = 'alt';        % alt key id

% helper functions
rotate_full = @(cp,sp,cy,sy) [cy cp*sy -sp*sy; -sy cp*cy -cy*sp; 0 sp cp];
rotate = @(r) rotate_full(cos(r(1)),sin(r(1)),cos(r(2)),sin(r(2)));

% global variables
running = true;     % whether the program is running
lock_mouse = false; % whether to lock the mouse to the center of the figure
pos = [0 0 2.5];      % camera position (x,y,z) in block lengths
vel = [0 0 0];      % camera velocity (x,y,z) in block lengths / second
rot = [0 0];        % camera rotation (pitch, yaw) in degrees
forward = [0 0 -1]; % camera forward direction
up = [0 1 0];       % camera up direction
                    % camera right direction

% initialize key map
keys = {KEY_SHIFT,KEY_W,KEY_A,KEY_S,KEY_D,KEY_SPACE};
keymap = containers.Map(keys,zeros(1,length(keys)));

% initialize chunk map
chunk_map = containers.Map("KeyType","int64","ValueType","any");

% initialize window
screen_h = get(0,'screensize');
screen_h = screen_h(4);
window_center = [0 0];
window = figure("CloseRequestFcn",@on_close,"KeyPressFcn",@on_key_down,...
    "KeyReleaseFcn",@on_key_up,"WindowButtonDownFcn",@on_mouse_down,...
    "ResizeFcn",@on_resize,"WindowButtonUpFcn",@on_mouse_up);
on_resize(window);
datacursormode off;
set(window,"MenuBar","none");
set(window,"ToolBar","none");

% initialize axis
axis vis3d;
ax = window.CurrentAxes;
set(ax,'visible','off');
disableDefaultInteractivity(ax);
daspect([1 1 1]);
xlim([-RANGE,RANGE]+pos(1));
ylim([-RANGE,RANGE]+pos(2));
zlim([-RANGE,RANGE]+pos(3));

% initialize camera directions
T = rotate(rot);
forward = (T*(forward'))';
up = (T*(up'))';
right = cross(forward,up);

% initialize camera
camproj("perspective");
camva(FOV);
set(ax,'cameraposition',pos);
set(ax,'cameratarget',pos+forward);
set(ax,'cameraupvector',up);

% initialize mouse
mouse = Robot;
mouse_offset = [0 0];

%% Game Loop
generate;
drawnow;
draw_tick = 0;
update_tick = 0;
while running
    dt = toc-update_tick;
    if dt>=UPDATE_LIMIT
        % update game state
        update_tick = toc;
        update(dt);
        if lock_mouse
            mouse_offset = get(0,"PointerLocation") - [window_center(1)+1 screen_h-window_center(2)];
            if mouse_offset(1)~=0 || mouse_offset(2)~=0
                mouse.mouseMove(window_center(1),window_center(2));
            end
        end
    end
    dt = toc-draw_tick;
    if dt>=DRAW_LIMIT
        % redraw game
        draw_tick = toc;
        target = pos+forward;
        if ~isequal(get(ax,"cameraposition"),pos)
            set(ax,"cameraposition",pos);
            xlim([-RANGE,RANGE]+pos(1));
            ylim([-RANGE,RANGE]+pos(2));
            zlim([-RANGE,RANGE]+pos(3));
        end
        if ~isequal(get(ax,"cameratarget"),target)
            set(ax,"cameratarget",target);
        end
        if ~isequal(get(ax,"cameraupvector"),up)
            set(ax,"cameraupvector",up);
        end
        %fprintf("dt:\t%d\n",dt*1000);
    end
    drawnow;
end
close all;

%% Terrain Generation
function blocks = generate_chunk(cx,cy,cz)
    %fprintf("generating chunk (%d,%d,%d)\n",cx,cy,cz);
    blocks = zeros(CHUNK_SIZE,CHUNK_SIZE,CHUNK_SIZE);
    type = @(x,y,z) z <= sin(x/10)+cos(y/10);
    for x=1:CHUNK_SIZE
        for y=1:CHUNK_SIZE
            for z=1:CHUNK_SIZE
                bid = type(x+cx*CHUNK_SIZE,y+cy*CHUNK_SIZE,z+cz*CHUNK_SIZE);
                if bid
                    blocks(x,y,z) = bid;
                end
            end
        end
    end
    chunk_map(cx) = containers.Map("KeyType","int64","ValueType","any");
    cmx = chunk_map(cx);
    cmx(cy) = containers.Map("KeyType","int64","ValueType","any");
    cmy = cmx(cy);
    cmy(cz) = blocks;
end

function blocks = get_chunk(cx,cy,cz)
    %fprintf("getting chunk (%d,%d,%d)\n",cx,cy,cz);
    if ~isKey(chunk_map,cx)
        blocks = generate_chunk(cx,cy,cz);
        return;
    end
    cmx = chunk_map(cx);
    if ~isKey(cmx,cy)
        blocks = generate_chunk(cx,cy,cz);
        return;
    end
    cmy = cmx(cy);
    if ~isKey(cmy,cz)
        blocks = generate_chunk(cx,cy,cz);
        return;
    end
    blocks = cmy(cz);
end

function [xl,yl,zl,cl] = get_block_faces(b,bd,bu,bl,br,bf,bb)
    xl=[];yl=[];zl=[];cl=[];
    % bottom face
    if ~bd
        xl = [xl [-1;-1;0;0]];
        yl = [yl [-1;0;0;-1]];
        zl = [zl [-1;-1;-1;-1]];
        cl = [cl [0.4;0.4;0.4;0.4]];
    end
    % top face
    if ~bu
        xl = [xl [-1;-1;0;0]];
        yl = [yl [-1;0;0;-1]];
        zl = [zl [0;0;0;0]];
        cl = [cl [1;1;1;1]];
    end
    % left face
    if ~bl
        xl = [xl [-1;-1;-1;-1]];
        yl = [yl [-1;0;0;-1]];
        zl = [zl [0;0;-1;-1]];
        cl = [cl [0.2;0.2;0.2;0.2]];
    end
    % right face
    if ~br
        xl = [xl [0;0;0;0]];
        yl = [yl [-1;0;0;-1]];
        zl = [zl [0;0;-1;-1]];
        cl = [cl [0.8;0.8;0.8;0.8]];
    end
    % front face
    if ~bf
        xl = [xl [-1;-1;0;0]];
        yl = [yl [-1;-1;-1;-1]];
        zl = [zl [0;-1;-1;0]];
        cl = [cl [0;0;0;0]];
    end
    % back face
    if ~bb
        xl = [xl [-1;-1;0;0]];
        yl = [yl [0;0;0;0]];
        zl = [zl [0;-1;-1;0]];
        cl = [cl [0.6;0.6;0.6;0.6]];
    end
end

function [xl,yl,zl,cl] = get_chunk_faces(cx,cy,cz)
    xl=[];yl=[];zl=[];cl=[];
    bs = get_chunk(cx,cy,cz);
    nd = get_chunk(cx,cy,cz-1);
    nu = get_chunk(cx,cy,cz+1);
    nl = get_chunk(cx-1,cy,cz);
    nr = get_chunk(cx+1,cy,cz);
    nf = get_chunk(cx,cy-1,cz);
    nb = get_chunk(cx,cy+1,cz);
    % inner chunk
    for x=1:CHUNK_SIZE
        for y=1:CHUNK_SIZE
            for z=1:CHUNK_SIZE
                b = bs(x,y,z);
                if b
                    % get l/r blocks
                    if x==1
                        bl = nl(CHUNK_SIZE,y,z);
                        br = bs(x+1,y,z);
                    elseif x==CHUNK_SIZE
                        bl = bs(x-1,y,z);
                        br = nr(1,y,z);
                    else
                        bl = bs(x-1,y,z);
                        br = bs(x+1,y,z);
                    end
                    % get f/b blocks
                    if y==1
                        bf = nf(x,CHUNK_SIZE,z);
                        bb = bs(x,y+1,z);
                    elseif y==CHUNK_SIZE
                        bf = bs(x,y-1,z);
                        bb = nb(x,1,z);
                    else
                        bf = bs(x,y-1,z);
                        bb = bs(x,y+1,z);
                    end
                    % get u/d blocks
                    if z==1
                        bu = bs(x,y,z+1);
                        bd = nd(x,y,CHUNK_SIZE);
                    elseif z==CHUNK_SIZE
                        bu = nu(x,y,1);
                        bd = bs(x,y,z-1);
                    else
                        bu = bs(x,y,z+1);
                        bd = bs(x,y,z-1);
                    end
                    % append faces
                    [fx,fy,fz,fc] = get_block_faces(b,bd,bu,bl,br,bf,bb);
                    xl = [xl x+fx];
                    yl = [yl y+fy];
                    zl = [zl z+fz];
                    cl = [cl fc];
                end
            end
        end
    end
    xl = xl+CHUNK_SIZE*(cx-0.5);
    yl = yl+CHUNK_SIZE*(cy-0.5);
    zl = zl+CHUNK_SIZE*(cz-0.5);
end

function generate
    x=[];y=[];z=[];c=[];
    for x=-2:2
        for y=-2:2
            for z=-2:2
                [fx,fy,fz,fc] = get_chunk_faces(x,y,z);
                patch(ax,fx,fy,fz,fc,"EdgeColor","none");
            end
        end
    end
end

%% State Updating
function update(dt)
    % update game state
    %fprintf("ut:\t%d\n",dt*1000);
    if mouse_offset(1) || mouse_offset(2)
        rot = rot + ROTATION_SPEED*[mouse_offset(2) mouse_offset(1)];
        T = rotate(rot);
        forward = (T*([0;0;-1]))';
        up = (T*([0;1;0]))';
        right = cross(forward,up);
    end
    vel = forward*(keymap(KEY_W)-keymap(KEY_S))+...
        right*(keymap(KEY_D)-keymap(KEY_A))+...
        up*(keymap(KEY_SPACE)-keymap(KEY_SHIFT));
    if vel(1)~=0 || vel(2)~=0 || vel(3)~=0
        pos = pos + MOVE_SPEED*dt*vel/norm(vel);
    end
end

%% Event Binders
function on_close(~,~)
    % figure termination event
    running = false;
    delete(window);
end

function on_key_down(~,event)
    % key down event
    key = string(event.Key);
    if ismember(key,keys) && ~keymap(key)
        keymap(key) = true;
    elseif (strcmp(key,KEY_ESCAPE) || strcmp(key,KEY_ALT)) && lock_mouse
        lock_mouse = false;
        mouse_offset = [0 0];
        mouse.mouseMove(window_center(1),window_center(2));
        set(window,"Pointer","arrow");
    end
end

function on_key_up(~,event)
    % key up event
    key = string(event.Key);
    if ismember(key,keys) && keymap(key)
        keymap(key) = false;
    end
end

function on_mouse_down(~,~)
    % figure mouse down event
    if ~lock_mouse
        lock_mouse = true;
        set(window,"Pointer","custom","PointerShapeCData",NaN(16,16));
    end
end

function on_mouse_up(src,~)
    % figure mouse up event
    on_resize(src);
end

function on_resize(src,~)
    % figure resize event
    w_pos = get(src,"Position");
    window_center = floor([w_pos(1)+w_pos(3)/2 screen_h-w_pos(2)-w_pos(4)/2]);
end
end