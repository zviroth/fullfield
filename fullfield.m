% fullfield.m
%      usage: fullfield
%         by: zvi
%       date: 04/6/2017
%    purpose: test for attentional modulation of orientation and spatial
%    frequency tuning curves and population responses
%
%       example: fullfield('targetContrastFactor = 1.8','attendRight=1');
%       targetContrastFactor should be <1.7`
%
%    TODO: 1) add subject responses and feedback
%           2) check whether targetContrastFactor given by user is ok
%   
%
%
function [] = fullfield(varargin);

% check arguments`
if ~any(nargin == [0:10])
    help fullFieldGrating`
    return
end

getArgs(varargin,[],'verbose=1');


if ieNotDefined('attendRight'), attendRight = 1; end %should the subject attend to the right target (or the left)
if ieNotDefined('target1segment'), target1segment = 1; end %should the target appear in only 1 segment per trial
if ieNotDefined('targetContrastFactor'), targetContrastFactor = 1.5; end 
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
% task{1}.seglen = 0.25*ones(1,6); 
% stimSegs = [0.25 0.05 0.25 0.05 0.25 0.05 0.25];
% stimSegs = [0.2 0.02 0.2 0.02 0.2 0.02 0.2];
stimSegs = [0.2 0.1 0.2 0.1 0.2 0.1 0.2];
respSeg = 2 - sum(stimSegs) - 0.2;
task{1}.seglen = [stimSegs respSeg];

task{1}.getResponse = ones(1,length(task{1}.seglen));%subject can respond anytime during the trial
% task{1}.getResponse = zeros(1,length(task{1}.seglen));
% task{1}.getResponse(length(task{1}.seglen)) = 1;%only accept a response during the blank period


% Stim1 is the LEFT stimulus
task{1}.randVars.uniform.tiltStim1 = [-1 1];
% Stim2 is the RIGHT stimulus
task{1}.randVars.uniform.tiltStim2 = [-1 1];
task{1}.random = 1;
task{1}.respCorrect = 0;
task{1}.respWrong = 0;

for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

stimulus.attendRight = attendRight;
stimulus.target1segment = target1segment;
stimulus.targetContrastFactor = targetContrastFactor;
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

stimulus.fixColor = [0 0 0];
% stimulus.fixColor = [0 0 255];


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

    if ~stimulus.target1segment || task.thistrial.thisseg == stimulus.targetSeg(task.trialnum) %present target
%     if task.thistrial.thisseg == stimulus.targetSeg(task.trialnum) %present target
        for hemi=1:2
            mglBltTexture(stimulus.targetTex(hemi,task.trialnum,ceil(task.thistrial.thisseg/2)), [0 0], 0, 0, 0);
        end
    end
end


mglFixationCross(0.5, 2, stimulus.fixColor, [0 0]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets subject  response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = responseCallback(task, myscreen)
global stimulus;
% get contrast change of attended target
if stimulus.attendRight
    attContrast = stimulus.targetContrastSign1(task.trialnum);% +1 or -1
else
    attContrast = stimulus.targetContrastSign2(task.trialnum);% +1 or -1
end
% convert contrast change (-1,1) into button responses (1,2)
correctButton = attContrast/2 + 1.5;

% make sure we have not already received a response
if task.thistrial.gotResponse==0
    if (task.thistrial.whichButton==correctButton)
        % play the correct sound
        mglPlaySound('Tink');
        task.respCorrect = task.respCorrect+1;
        stimulus.fixColor = [0 255 0];
    else
        % play the incorrect sound
        mglPlaySound('Pop');
        task.respWrong = task.respWrong+1;
        stimulus.fixColor = [255 0 0];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the grating stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stimulus = myInitStimulus(stimulus, myscreen, task)


% keep an array that lists which of the segments we
% are presenting the stimulus in.
% stimulus.stimulusSegments = [1:2:8];`
stimulus.stimulusSegments = [1:2:length(task{1}.seglen)];



%create stimulus order
ntrials=160; %should be a whole multiple of the number of stimulus (background) types
phases = 1:4;
numOris = 8;
numFreqs = 5; 
numContrasts =4;

% %create stimulus order
% ntrials=8; %should be a whole multiple of the number of stimulus (background) types
% phases = 1:4;
% numOris = 2;
% numFreqs = 2; 
% numContrasts =2;

stimulus.ntrials = ntrials;
stimulus.randOrder = randperm(ntrials);
stimulus.randTargetOrder = randperm(ntrials);

oris = [];
for i=1:numOris
   oris = [oris i*ones(1,ntrials/numOris)]; 
end
contrasts = [];
for i=1:numContrasts
    contrasts = [contrasts i*ones(1,numOris/numContrasts)];
end
contrasts = repmat(contrasts,1,ntrials/numOris);




freqs = 1:numFreqs;
freqs = repmat(freqs,1,ntrials/numFreqs);

%change from ordinal to meaningful stimulus parameters
phases = (phases-1) * 360/max(phases);
oris = (oris-1) * 180/max(oris);


% contrasts = 0.1 * 2 .^ (contrasts-1);
contrasts = 0.1 * 1.86 .^ (contrasts-1);
freqs = 0.1 * 2.5 .^ (freqs-1);



% check if the target contrasts is above 1 or below 0.06
if stimulus.targetContrastFactor * max(contrasts) > 1
   fprintf('Target contrast will be too high!\n');
end
if min(contrasts)/stimulus.targetContrastFactor <0.06
   fprintf('Target contrast will be too low!\n');
end

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
%replicate the locations to be equal to the number of trials (10*16=160)
temp=locations;
for i=1:(ntrials/length(temp)-1)
   locations =  [temp; locations];
end
nlocations = size(locations,1);

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
stimulus.targetSize = 2;%diameter?`
stimulus.transition = 0.25;
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
    xLocationPix = locations(i,1)*xDeg2pix;% 
    yLocationPix = locations(i,2)*yDeg2pix;% 
    % [locationPix] = locations(i,:) .* [xDeg2pix yDeg2pix];
    disc{i,1}(:,:) = mkDisc([heightPixels widthPixels], ... % matrix size
        stimulus.targetSize*xDeg2pix/2, ... % radius
        [heightPixels widthPixels]/2 + [yLocationPix xLocationPix], ... % !!!! Right origin
        stimulus.transition * xDeg2pix); % transition
    
    disc{i,2}(:,:) = mkDisc([heightPixels widthPixels], ... % matrix size
        stimulus.targetSize*xDeg2pix/2, ... % radius
        [heightPixels widthPixels]/2 + [yLocationPix -xLocationPix], ... % !!!! Left origin
        stimulus.transition * xDeg2pix); % transition
end


%%
stimulus.targetLocation = locations;
stimulus.targetContrastSign1 = sign(rand(nlocations,1)-0.5);
stimulus.targetContrastSign2 = sign(rand(nlocations,1)-0.5);
% stimulus.targetContrastDiff = 0.2;
% stimulus.targetContrastDiff = targetContrastDiff;

% make a full-field grating for each set of parameters
tic
for i=1:stimulus.ntrials
    bgcontrast = contrasts(stimulus.randOrder(i)); %the background contrast
    randPhaseOrder = randperm(length(phases));
    for j=1:length(phases)
        % make a gabor patch
        grating = mglMakeGrating(stimulus.textureWidth, stimulus.textureHeight, freqs(stimulus.randOrder(i)), oris(stimulus.randOrder(i)), phases(randPhaseOrder(j)));

        %create the target gratings with higher/lower contrast
        %for positive targetContrastSign increase bgcontrast by factor of
        %targetContrastFactor, for negative decrease by same factor, so the
        %target contrasts multiplied by targetContrastFactor equals bgcontrast.
        targetGrating1 = grating*(bgcontrast * stimulus.targetContrastFactor^stimulus.targetContrastSign1(i)); %Right
        targetGrating2 = grating*(bgcontrast * stimulus.targetContrastFactor^stimulus.targetContrastSign2(i)); %Left
        grating = grating*bgcontrast;

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
