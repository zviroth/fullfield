% fullFieldGrating.m
%      usage: fullFieldGrating
%         by: zvi
%       date: 04/6/2017
%    purpose: test for attentional modulation of orientation and spatial
%    frequency tuning curves and population responses
%
%
%    TODO: 1) add target stimulus and task
%           2 
%
%
function [] = contrastPatch2(varargin);

% check arguments`
if ~any(nargin == [0:10])
    help fullFieldGrating
    return
end

getArgs(varargin,[],'verbose=1');


if ieNotDefined('contrast'), contrast = 1; end
% if ieNotDefined('easyTilt'), easyTilt = 20;

% init screen and open up the window
myscreen.background = [127 127 127];
myscreen.autoCloseScreen = 0;
myscreen.saveData = 1;
myscreen = initScreen(myscreen);



% init the task
task{1}.waitForBacktick = 1;
% 8 segments, 200ms stimulus, 150ms blank, etc.
% task{1}.seglen = [0.2 0.15 0.2 0.15 0.2 0.15 0.2 0.15]; 
task{1}.seglen = 0.25*ones(1,6); 
stimSegs = [0.25 0.05 0.25 0.05 0.25 0.05 0.25];
respSeg = 2 - sum(stimSegs) - 0.2;
task{1}.seglen = [stimSegs respSeg];
task{1}.getResponse = ones*length(task{1}.seglen);
% we test threshold relative to reference

% Stim1 is the LEFT stimulus
task{1}.randVars.uniform.tiltStim1 = [-1 1];
% Stim2 is the RIGHT stimulus
task{1}.randVars.uniform.tiltStim2 = [-1 1];
task{1}.random = 1;
task{1}.tiltDeg = 15;

for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
global stimulus;



myscreen = initStimulus('stimulus',myscreen);
stimulus.contrast = contrast;
% stimulus.easyTilt = easyTilt;
stimulus = myInitStimulus(stimulus,myscreen,task);
task{1}.numTrials = stimulus.ntrials;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % run the eye calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);
mglClearScreen;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end


% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus;


stimulus.fixColor = [255 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus;
% clear the screen
mglClearScreen;
if any(task.thistrial.thisseg == stimulus.stimulusSegments)
    mglBltTexture(stimulus.tex(task.trialnum, ceil(task.thistrial.thisseg/2)), [0 0], 0, 0, 0);
    targetX = stimulus.targetLocation(task.trialnum,1);
    targetY = stimulus.targetLocation(task.trialnum,2);
    
    %     task.thistrial.thisOrientationStim1 = task.thistrial.tiltStim1;
    %     task.thistrial.thisOrientationStim2 = task.thistrial.tiltStim2;
    %     lineCenterRight = stimulus.targetLocation(task.trialnum,:);
    %     lineCenterLeft = [-lineCenterRight(1) lineCenterRight(2)];
    %     rightRotation =task.thistrial.thisOrientationStim2 * task.tiltDeg * pi/180;
    %     leftRotation = task.thistrial.thisOrientationStim1 * task.tiltDeg * pi/180;
    %     drawOriLine(lineCenterRight, stimulus.targetLength, stimulus.targetWidth, rightRotation, [0 1 0])
    %     drawOriLine(lineCenterLeft, stimulus.targetLength, stimulus.targetWidth, leftRotation, [0 1 0])
    
%     if task.thistrial.thisseg == stimulus.targetSeg(task.trialnum) %present target
        for hemi=1:2
            mglBltTexture(stimulus.targetTex(hemi,task.trialnum,ceil(task.thistrial.thisseg/2)), [0 0], 0, 0, 0);
        end
%     end
end


mglFixationCross(0.5, 2, stimulus.fixColor, [0 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets subject  response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task, myscreen)
global stimulus;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the grating stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stimulus = myInitStimulus(stimulus, myscreen, task)


% keep an array that lists which of the segments we
% are presenting the stimulus in.
% stimulus.stimulusSegments = [1:2:8];`
stimulus.stimulusSegments = [1:2:length(task{1}.seglen)];

%create stimulus order
ntrials=32; %should be a whole multiple of the number of stimulus (background) types
stimulus.ntrials = ntrials;
stimulus.randOrder = randperm(ntrials);
stimulus.randTargetOrder = randperm(ntrials);
phases = 1:4;
oris = [ones(1,8) 2*ones(1,8) 3*ones(1,8) 4*ones(1,8)];
% oris = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4];
contrasts = [ones(1,4) 2*ones(1,4)];
contrasts = repmat(contrasts,1,4);
% contrasts = [1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2];
freqs = 1:4;
freqs = repmat(freqs,1,8);
% freqs = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];



% %create stimulus order
% ntrials=8; %should be a whole multiple of the number of stimulus (background) types
% stimulus.ntrials = ntrials;
% stimulus.randOrder = randperm(ntrials);
% stimulus.randTargetOrder = randperm(ntrials);
% phases = 1:4;
% oris = [ones(1,8)];
% % oris = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4];
% contrasts = [ones(1,4) 2*ones(1,4)];
% % contrasts = repmat(contrasts,1,4);
% % contrasts = [1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2];
% freqs = 1:4;
% freqs = repmat(freqs,1,2);
% % freqs = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];



%change from ordinal to meaningful stimulus parameters
phases = (phases-1) * 360/max(phases);
oris = (oris-1) * 180/max(oris);
targetContrastDiff = 1/(1+max(contrasts));
contrasts = contrasts * 1/(1+max(contrasts));
% targetContrastDiff=0.2;
% contrasts = targetContrastDiff+(contrasts-1) * 1/(max(contrasts));
% contrasts = contrasts*0.8 + 0.2;
% contrasts = contrasts.^2;
freqs = freqs * 1/max(freqs) + 0.6;

% target locations, in the right hemifield
X = 1:1:4;
i=0;
for x=X
    for y= (-x+1):2:(x-1)
        i=i+1;
        locations(i,:) = [x y];
    end
end
xscale=1.3;
yscale=1;
locations(:,1) = locations(:,1) * xscale;
locations(:,2) = locations(:,2) * yscale;
%replicate the locations
locations = [locations; locations; locations; locations]; %should have at least ntrials locations
nlocations = size(locations,1);%this should be made to equal the number of trials
%randomize location order
locations = locations(randperm(nlocations),:);
%%
% figure
% scatter(locations(:,1), locations(:,2),'filled');
% keyboard
%%
% background size
stimulus.textureHeight = 18;
stimulus.textureWidth = 25;

% target size properties (in deg)
stimulus.targetSize = 1.5;
stimulus.transition = 0.25;`
% stimulus.textureSize = 1.5;

% make the window for the raised cosine
xDeg2pix = mglGetParam('xDeviceToPixels');
yDeg2pix = mglGetParam('yDeviceToPixels');
widthPixels = round(stimulus.textureWidth*xDeg2pix);%width of the fullfield
heightPixels = round(stimulus.textureHeight*yDeg2pix);%height of the fullfield
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);
%% Make a disc for every location.
% remember, x and y are flipped
for i=1:nlocations
    yLocationPix = locations(i,1)*xDeg2pix;% !!!!!
    xLocationPix = locations(i,2)*yDeg2pix;% !!!!!
    % [locationPix] = locations(i,:) .* [xDeg2pix yDeg2pix];
    disc{i,1}(:,:) = mkDisc([heightPixels widthPixels], ... % matrix size
        stimulus.targetSize*xDeg2pix/2, ... % radius
        [heightPixels widthPixels]/2 + [xLocationPix yLocationPix], ... % origin
        stimulus.transition * xDeg2pix); % transition
    
    disc{i,2}(:,:) = mkDisc([heightPixels widthPixels], ... % matrix size
        stimulus.targetSize*xDeg2pix/2, ... % radius
        [heightPixels widthPixels]/2 + [xLocationPix -yLocationPix], ... % origin
        stimulus.transition * xDeg2pix); % transition
end


%%
stimulus.targetLocation = locations;
stimulus.targetContrastSign1 = sign(rand(nlocations,1)-0.5);
stimulus.targetContrastSign2 = sign(rand(nlocations,1)-0.5);
% stimulus.targetContrastDiff = 0.2;
stimulus.targetContrastDiff = targetContrastDiff/2;
% make a full-field grating for each set of parameters
tic
for i=1:stimulus.ntrials
    randPhaseOrder = randperm(length(phases));
    for j=1:length(phases)
        % make a gabor patch
        grating = mglMakeGrating(stimulus.textureWidth, stimulus.textureHeight, freqs(stimulus.randOrder(i)), oris(stimulus.randOrder(i)), phases(randPhaseOrder(j)));

        %create the target gratings with higher/lower contrast
        targetGrating1 = grating*(contrasts(stimulus.randOrder(i)) + stimulus.targetContrastDiff*stimulus.targetContrastSign1(i))^2;
        targetGrating2 = grating*(contrasts(stimulus.randOrder(i)) + stimulus.targetContrastDiff*stimulus.targetContrastSign2(i))^2;
        grating = grating*contrasts(stimulus.randOrder(i))^2;

        % scale it to be between 0 and 255
        grating = 255*(grating+1)/2;
        targetGrating1 = 255*(targetGrating1+1)/2;
        targetGrating2 = 255*(targetGrating2+1)/2;

        % make it into a texture
        stimulus.tex(i,j) = mglCreateTexture(grating);
        %change to RGB + alpha format
        alphaGrating1 = repmat(targetGrating1,1,1,3);
        alphaGrating1(:,:,4) = 255*disc{i,1};
        alphaGrating2 = repmat(targetGrating2,1,1,3);
        alphaGrating2(:,:,4) = 255*disc{i,2};

        stimulus.targetTex(1,i,j) = mglCreateTexture(alphaGrating1);
        stimulus.targetTex(2,i,j) = mglCreateTexture(alphaGrating2);

    end
    %choose segments to present the target patches
    temp = randperm(length(stimulus.stimulusSegments));
    stimulus.targetSeg(i) = stimulus.stimulusSegments(temp(1));
end
toc
