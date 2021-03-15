function minecraft
%% Initialization
clear; clc; close all; opengl hardware; import java.awt.Robot;
% global variables
UPDATE_LIMIT = 0;       % minimum time before game state is updated (performance)
DRAW_LIMIT = 0;         % minimum time before frame is redrawn (performance)
CHUNK_SIZE = 15;        % chunk width in block lengths (requires int)
SEED = tic;             % unique game seed
RANGE = 3;              % camera range in chunk lengths
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
LIGHT_COLOR = [255 255 251]; % sunlight color
LIGHT_ANGLE = [0 35];   % sunlight angle (az, el)
AMBIENT = 1;            % ambient light strength
SPECULAR = 0.1;         % specular reflection strength
DIFFUSE = 1;            % diffuse reflection strength

% block ids
B_AIR =     0;
B_GRASS =   1;
B_DIRT =    2;
B_STONE =   3;

% block colors (bottom, top, left, right, front, back)
BLOCK_COLORS(B_GRASS,:,:) = [155 118 83;154 205 50;155 118 83;155 118 83;155 118 83;155 118 83];
BLOCK_COLORS(B_DIRT,:,:) =  [155 118 83;155 118 83;155 118 83;155 118 83;155 118 83;155 118 83];
BLOCK_COLORS(B_STONE,:,:) = [136 140 141;136 140 141;136 140 141;136 140 141;136 140 141;136 140 141];
BLOCK_COLORS = BLOCK_COLORS/255;

% helper functions
rotate_full = @(cp,sp,cy,sy) [cy cp*sy -sp*sy; -sy cp*cy -cy*sp; 0 sp cp];
rotate = @(r) rotate_full(cos(r(1)),sin(r(1)),cos(r(2)),sin(r(2)));
chunk_pos = @(p) ceil(p/CHUNK_SIZE-0.5); % get block position from chunk position

% global variables
running = true;     % whether the program is running
lock_mouse = false; % whether to lock the mouse to the center of the figure
pos = [0 0 0];      % camera position (x,y,z) in block lengths
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
chunk_surface_map = containers.Map("KeyType","int64","ValueType","any");

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
xlim([-RANGE,RANGE]*CHUNK_SIZE+pos(1));
ylim([-RANGE,RANGE]*CHUNK_SIZE+pos(2));
zlim([-RANGE,RANGE]*CHUNK_SIZE+pos(3));

% initialize camera directions
T = rotate(rot);
forward = (T*(forward'))';
up = (T*(up'))';
right = cross(forward,up);

% initialize camera
camproj("perspective");
camva(FOV);
set(ax,"cameraposition",pos);
set(ax,"cameratarget",pos+forward);
set(ax,"cameraupvector",up);
light = camlight(LIGHT_ANGLE(1),LIGHT_ANGLE(2));
set(light,"style","infinite");
set(light,"Color",LIGHT_COLOR/255);

% initialize mouse
mouse = Robot;
mouse_offset = [0 0];

%% Game Loop
generate_chunks;
tic;
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
        pos_last = get(ax,"cameraposition");
        % update position
        if ~isequal(pos_last,pos)
            set(ax,"cameraposition",pos);
            cpos = chunk_pos(pos);
            cpos_last = chunk_pos(pos_last);
            % crossed chunk border
            if ~isequal(cpos,cpos_last)
                draw_chunks(cpos,cpos_last);
                ra = ([-RANGE,RANGE]'+cpos+0.5)*CHUNK_SIZE;
                xlim(ra(:,1));
                ylim(ra(:,2));
                zlim(ra(:,3));
            end
        end
        % update look vector
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

%% Chunking
function blocks = get_chunk(cx,cy,cz)
    % returns the block data for the chunk at chunk position (cx,cy,cz)
    % assure cx is in map
    if ~isKey(chunk_map,cx)
        chunk_map(cx) = containers.Map("KeyType","int64","ValueType","any");
    end
    % assure cy is in map
    cmx = chunk_map(cx);
    if ~isKey(cmx,cy)
        cmx(cy) = containers.Map("KeyType","int64","ValueType","any");
    end
    % assure cz is in map
    cmy = cmx(cy);
    if ~isKey(cmy,cz)
        % generate block data
        ox = CHUNK_SIZE*(cx-0.5)+0.5;
        oy = CHUNK_SIZE*(cy-0.5)+0.5;
        oz = CHUNK_SIZE*(cz-0.5)+0.5;
        r = 1:CHUNK_SIZE;
        x = r+ox;
        y = (r+oy)';
        z = permute(r+oz,[3,1,2]);
        % assign block data
        blocks = z+2<=sin(x)+cos(y);
        cmy(cz) = blocks;
        return;
    end
    blocks = cmy(cz);
end

function surface = generate_chunk_faces(cx,cy,cz)
    % returns a patch object for the chunk at chunk position (cx,cy,cz)
    potential = 6*ceil((CHUNK_SIZE^3)/2);
    xl=zeros(4,potential);
    yl=zeros(4,potential);
    zl=zeros(4,potential);
    cl=zeros(1,potential,3);
    faces = 0;
    % get chunk and chunk neighbors
    bs = get_chunk(cx,cy,cz);
    nd = get_chunk(cx,cy,cz-1);
    nu = get_chunk(cx,cy,cz+1);
    nl = get_chunk(cx-1,cy,cz);
    nr = get_chunk(cx+1,cy,cz);
    nf = get_chunk(cx,cy-1,cz);
    nb = get_chunk(cx,cy+1,cz);
    % get chunk faces
    r = 1:CHUNK_SIZE;
    for x=r
        for y=r
            for z=r
                b = bs(x,y,z);
                if b
                    % get left/right neighbors
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
                    % get front/back neighbors
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
                    % get top/bottom neighbors
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
                    % bottom face
                    if ~bd
                        faces = faces+1;
                        xl(:,faces) = x+[-1;-1;0;0];
                        yl(:,faces) = y+[-1;0;0;-1];
                        zl(:,faces) = z+[-1;-1;-1;-1];
                        cl(:,faces,:) = BLOCK_COLORS(b,1,:);
                    end
                    % top face
                    if ~bu
                        faces = faces+1;
                        xl(:,faces) = x+[-1;-1;0;0];
                        yl(:,faces) = y+[-1;0;0;-1];
                        zl(:,faces) = z+[0;0;0;0];
                        cl(:,faces,:) = BLOCK_COLORS(b,2,:);
                    end
                    % left face
                    if ~bl
                        faces = faces+1;
                        xl(:,faces) = x+[-1;-1;-1;-1];
                        yl(:,faces) = y+[-1;0;0;-1];
                        zl(:,faces) = z+[0;0;-1;-1];
                        cl(:,faces,:) = BLOCK_COLORS(b,3,:);
                    end
                    % right face
                    if ~br
                        faces = faces+1;
                        xl(:,faces) = x+[0;0;0;0];
                        yl(:,faces) = y+[-1;0;0;-1];
                        zl(:,faces) = z+[0;0;-1;-1];
                        cl(:,faces,:) = BLOCK_COLORS(b,4,:);
                    end
                    % front face
                    if ~bf
                        faces = faces+1;
                        xl(:,faces) = x+[-1;-1;0;0];
                        yl(:,faces) = y+[-1;-1;-1;-1];
                        zl(:,faces) = z+[0;-1;-1;0];
                        cl(:,faces,:) = BLOCK_COLORS(b,5,:);
                    end
                    % back face
                    if ~bb
                        faces = faces+1;
                        xl(:,faces) = x+[-1;-1;0;0];
                        yl(:,faces) = y+[0;0;0;0];
                        zl(:,faces) = z+[0;-1;-1;0];
                        cl(:,faces,:) = BLOCK_COLORS(b,6,:);
                    end
                end
            end
        end
    end
    % no faces to render
    if ~faces
        surface = patch(ax,[],[],[],[],"EdgeColor","none");
        return
    end
    % reduce set
    xl = xl(:,1:faces);
    yl = yl(:,1:faces);
    zl = zl(:,1:faces);
    cl = cl(:,1:faces,:);
    % shift faces by global chunk position
    xl = xl+CHUNK_SIZE*(cx-0.5);
    yl = yl+CHUNK_SIZE*(cy-0.5);
    zl = zl+CHUNK_SIZE*(cz-0.5);
    % initialize chunk surface
    surface = patch(ax,xl,yl,zl,cl,"EdgeColor","none");
    set(surface,"DiffuseStrength",DIFFUSE);
    set(surface,"SpecularStrength",SPECULAR);
    set(surface,"AmbientStrength",AMBIENT);
end

function surface = get_chunk_faces(cx,cy,cz)
    % returns the block data for the chunk at chunk position (cx,cy,cz)
    if ~isKey(chunk_surface_map,cx)
        chunk_surface_map(cx) = containers.Map("KeyType","int64","ValueType","any");
    end
    cmx = chunk_surface_map(cx);
    if ~isKey(cmx,cy)
        cmx(cy) = containers.Map("KeyType","int64","ValueType","any");
    end
    cmy = cmx(cy);
    if ~isKey(cmy,cz)
        surface = generate_chunk_faces(cx,cy,cz);
        cmy(cz) = surface;
        return;
    end
    surface = cmy(cz);
end

function draw_chunks(cpos_new,cpos_last)
    % creates and draws chunks around the player
    rval = RANGE-1;
    r = -rval:rval;
    % unload previous
    for x=cpos_last(1)+r
        for y=cpos_last(2)+r
            for z=cpos_last(3)+r
                if x<cpos_new(1)-rval || x>cpos_new(1)+rval || ...
                        y<cpos_new(2)-rval || y>cpos_new(2)+rval ||...
                        z<cpos_new(3)-rval || z>cpos_new(3)+rval
                    s = get_chunk_faces(x,y,z);
                    set(s,"Visible","off");
                end
            end
        end
    end
    % load new
    for x=cpos_new(1)+r
        for y=cpos_new(2)+r
            for z=cpos_new(3)+r
                if x<cpos_last(1)-rval || x>cpos_last(1)+rval || ...
                        y<cpos_last(2)-rval || y>cpos_last(2)+rval ||...
                        z<cpos_last(3)-rval || z>cpos_last(3)+rval
                    s = get_chunk_faces(x,y,z);
                    set(s,"Visible","on");
                end
            end
        end
    end
end

function generate_chunks
    % initializes the block data for ...
    r = RANGE-1;
    r = -r:r;
    for x=r
        for y=r
            for z=r
                get_chunk_faces(x,y,z);
            end
        end
    end
end

%% State Updating
function update(dt)
    % update game state
    %fprintf("ut:\t%d\n",dt*1000);
    % update player rotation
    if mouse_offset(1) || mouse_offset(2)
        rot = rot + ROTATION_SPEED*[mouse_offset(2) mouse_offset(1)];
        T = rotate(rot);
        forward = (T*([0;0;-1]))';
        up = (T*([0;1;0]))';
        right = cross(forward,up);
    end
    % update player location
    vel = forward*(keymap(KEY_W)-keymap(KEY_S))+...
        right*(keymap(KEY_D)-keymap(KEY_A))+...
        [0 0 1]*(keymap(KEY_SPACE)-keymap(KEY_SHIFT));
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