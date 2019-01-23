%
%
%  Batch Process Function for Imaris 7.3.0
%
%  Copyright Bitplane AG 2011
%
%
%  Installation:
%
%  - Copy this file into the folder containing the XTensions you want to use
%  - You will find this function in the Image Processing menu
%
%    <CustomTools>
%      <Menu>
%        <Item name="Batch Process" icon="Matlab" tooltip="Automatically process images with a custom script">
%          <Command>MatlabXT::XTBatchProcess(%i)</Command>
%        </Item>
%      </Menu>
%    </CustomTools>
%  
%
%  Description:
%
%   This XTension batch processes images.
%
%

function PM_Spotsinsideaxon1(aImarisApplicationID)

%Get script directory
scriptname = mfilename;
scriptfullname = mfilename('fullpath');
idx = strfind(scriptfullname, scriptname);
scriptpath = scriptfullname(1:idx(size(idx,2))-1);

%Get .m files (provided functions names equal filenames)
matlabfiles = what(scriptpath);
[junk,mfiles] = cellfun(@fileparts,matlabfiles.m,'UniformOutput',0); %#ok
clear junk;

%Get the current function name
[junk,curfun] = fileparts(scriptname); %#ok
clear junk;

%Remove it from the functions list
curfunidx = ismember(mfiles,curfun);
mfiles(curfunidx) = [];

%Get the images folder
folder = uigetdir;
files = [folder '\*.ims'];
listing = dir(files);
nfiles = size(listing,1);

 %Connect to Imaris based on Imaris 7.6 XT Interface documentation
    if isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
        vImarisApp = aImarisApplicationID;
    else
        if exist('ImarisLib','class') == 0
            javaaddpath ImarisLib.jar
        end
        vImarisLib = ImarisLib;
        aServer = vImarisLib.GetServer;
        numObjects = aServer.GetNumberOfObjects;
        for vIndex = 0:numObjects-1
            aObjectId = aServer.GetObjectID(vIndex);
            break
        end
        if ischar(aObjectId)
            aObjectId = round(aObjectId);
        end
        vImarisApp = vImarisLib.GetApplication(aObjectId);
    end

    %Process the files
    for i=1:nfiles
        filename = [folder '\' listing(i).name];
        vImarisApp.FileOpen(filename, '');
        
        %Get Surpass Scene
        XT = EasyXT(0);
        XT.CreateNewScene()
        
        vImarisApp.GetSurpassCamera().SetOrientationAxisAngle([0,0,1],0);
        vImarisApp.GetSurpassCamera().Fit();
        %Get image
        aDataSet = vImarisApp.GetDataSet();
        
%Gaussian filter
        vSigma = 0.0517;
        vDataSetOut = vImarisApp.GetImageProcessing.GaussFilterDataSet(aDataSet,vSigma);
        vImarisApp.SetDataSet(vDataSetOut)
    
%Background substraction
        vSigma=13.2;
        vDataSetOut=vImarisApp.GetImageProcessing.SubtractBackgroundDataSet(aDataSet,vSigma);
        vImarisApp.SetDataSet(vDataSetOut);
        
%Normalize layers
         vDataSetIn=vImarisApp.GetDataSet; %% Current dataset
         vChannelIndex=1;
         vImarisApp.GetImageProcessing.NormalizeLayersChannel(vDataSetIn,vChannelIndex)
        vChannelIndex=2;
         vImarisApp.GetImageProcessing.NormalizeLayersChannel(vDataSetIn,vChannelIndex);
        
%Localize spots
        XT = EasyXT(0);
        channel = 2;
        newSpots = XT.DetectSpots(channel, 'Name', 'GAD65', ...
                                  'Diameter XY', 0.6, ...
                                  'Spots Filter', '"Quality" above 9.16' ...
                         );
                     XT.SetColor(newSpots, [0 255 255]);
                     XT.AddToScene(newSpots);
%Localize spots
        XT = EasyXT(0);
        channel = 1;
        newSpots = XT.DetectSpots(channel, 'Name', 'geph', ...
                                  'Diameter XY', 0.3, ...
                                  'Spots Filter', '"Quality" above 16.1' ...
                         );
                     XT.SetColor(newSpots, [0 255 0]);
                     XT.AddToScene(newSpots);
         
%Detect Surface
        XT = EasyXT(0);
        channel =3;
        surface = XT.DetectSurfaces(channel, 'Name', 'mcherry', ...
                                   'Smoothing', 0.05, ...
                                   'Threshold', 46.2, ...
                                  'Filter', '"Volume" above 0.02 um^3'...
                                   ); 
        XT.SetColor(surface, [255 0 0]);
        XT.AddToScene(surface);
%% Colocalize Spots
% get the spots
vSpots = vImarisApp.GetSurpassSelection;
vSpotsSelected = vImarisApp.GetFactory.IsSpots(vSpots);

if vSpotsSelected
    vScene = vSpots.GetParent;
else
    vScene = vImarisApp.GetSurpassScene;
end
vNumberOfSpots = 0;
vSpotsList{vScene.GetNumberOfChildren} = [];
vNamesList{vScene.GetNumberOfChildren} = [];
for vChildIndex = 1:vScene.GetNumberOfChildren
    vDataItem = vScene.GetChild(vChildIndex - 1);
    if vImarisApp.GetFactory.IsSpots(vDataItem)
        vNumberOfSpots = vNumberOfSpots+1;
        vSpotsList{vNumberOfSpots} = vImarisApp.GetFactory.ToSpots(vDataItem);
        vNamesList{vNumberOfSpots} = char(vDataItem.GetName);
    end
end

if vNumberOfSpots<2
    msgbox('Please create at least 2 spots objects!');
    return;
end

vSpots1 = XT.GetObject('Name', 'GAD65');
vSpots2 = XT.GetObject('Name', 'geph');
% ask for threshold
%if nargin<2
    %vThreshold = 0.38;
%else
    %vThreshold = aThreshold;
%end
vThreshold = 0.38;
vThresholdSquare = vThreshold.^2;

vProgressDisplay = waitbar(0,'Colocalizing spots');

vSpotsXYZ1 = vSpots1.GetPositionsXYZ;
vTime1 = vSpots1.GetIndicesT;
vRadius1 = vSpots1.GetRadiiXYZ;

vSpotsXYZ2 = vSpots2.GetPositionsXYZ;
vTime2 = vSpots2.GetIndicesT;
vRadius2 = vSpots2.GetRadiiXYZ;

% initialize coloc to zero
vColoc1 = false(numel(vTime1), 1);
vColoc2 = false(numel(vTime2), 1);

vTime1 = double(vTime1);
vTime2 = double(vTime2);

vStart = max([min(vTime1), min(vTime2)]);
vEnd = min([max(vTime1), max(vTime2)]);
for vTime = vStart:vEnd
    vValid1 = find(vTime1 == vTime);
    vValid2 = find(vTime2 == vTime);

    vXYZ = vSpotsXYZ2(vValid2, :);
    for vSpot1 = 1:numel(vValid1)
        vColocated1 = vValid1(vSpot1);
        
        vX = vXYZ(:, 1) - vSpotsXYZ1(vColocated1, 1);
        vY = vXYZ(:, 2) - vSpotsXYZ1(vColocated1, 2);
        vZ = vXYZ(:, 3) - vSpotsXYZ1(vColocated1, 3);
        vDistanceList = vX.^2 + vY.^2 + vZ.^2 <= vThresholdSquare;
        
        vColocated2 = vValid2(vDistanceList);
        if ~isempty(vColocated2)
            vColoc1(vColocated1) = true;
            vColoc2(vColocated2) = true;
        end
    end
    
    waitbar((vTime-vStart+1)/(vEnd-vStart+1), vProgressDisplay);
end

close(vProgressDisplay);

if isempty(find(vColoc1, 1))
    msgbox('There is no colocated spots.');
    return
end
vNonColoc1 = ~vColoc1;
vNonColoc2 = ~vColoc2;
% create new group
%vSpotsGroup = vImarisApp.GetFactory.CreateDataContainer;
%vSpotsGroup.SetName(sprintf('Coloc[%.2f] %s | %s', ...
    %vThreshold, char(vSpots1.GetName), char(vSpots2.GetName)));

vNewSpots1 = vImarisApp.GetFactory.CreateSpots;
vNewSpots1.Set(vSpotsXYZ1(vColoc1, :), vTime1(vColoc1), zeros(sum(vColoc1),1));
vNewSpots1.SetRadiiXYZ(vRadius1(vColoc1,:));
vNewSpots1.SetName([char(vSpots1.GetName), ' colocated']);
vRGBA = vSpots1.GetColorRGBA;
vNewSpots1.SetColorRGBA(vRGBA);
XT.AddToScene(vNewSpots1);

vNewSpots2 = vImarisApp.GetFactory.CreateSpots;
vNewSpots2.Set(vSpotsXYZ2(vColoc2, :), vTime2(vColoc2), zeros(sum(vColoc2),1));
vNewSpots2.SetRadiiXYZ(vRadius2(vColoc2,:));
vNewSpots2.SetName([char(vSpots2.GetName),' colocated']);
vRGBA = vSpots2.GetColorRGBA;
vNewSpots2.SetColorRGBA(vRGBA);
XT.AddToScene(vNewSpots2);

vSpots1.SetVisible(0);
vSpots2.SetVisible(0);
   %%
        channelOutputImage = strcat(filename, '_V1.ims');
        vImarisApp.FileSave(channelOutputImage, 'writer="Imaris5"');
      
    end

end