%
%  Note that you need to unify the surface before running this second script!
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

function PM_Spotsinsideaxon2  (aImarisApplicationID)

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
%%
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
    %%
    for i=1:nfiles
        filename = [folder '\' listing(i).name];
        vImarisApp.FileOpen(filename, ...
                            'reader="Imaris5"' ...
                                );
        vImarisApp.GetSurpassCamera().SetOrientationAxisAngle([0,0,1],0);
        vImarisApp.GetSurpassCamera().Fit();
       
%% SPLIT SPOTS 1
XT = EasyXT (0);
bout1 = XT.GetObject('Name', 'GAD65');
    
for aObjectId = 0;
    vSpots = vImarisApp.GetFactory.ToSpots(bout1);
vSurfaces = vImarisApp.GetFactory.ToSurfaces(vImarisApp.GetSurpassSelection);

vSurfacesSelected = vImarisApp.GetFactory.IsSurfaces(vSurfaces);

if vSurfacesSelected
    vScene = vSurfaces.GetParent;
else
    vScene = vImarisApp.GetSurpassScene;
end
vNumberOfSurfaces = 0;
vSurfacesList{vScene.GetNumberOfChildren} = [];
vNamesList{vScene.GetNumberOfChildren} = [];
for vChildIndex = 1:vScene.GetNumberOfChildren
    vDataItem = vScene.GetChild(vChildIndex - 1);
    if vImarisApp.GetFactory.IsSurfaces(vDataItem)
        vNumberOfSurfaces = vNumberOfSurfaces+1;
        vSurfacesList{vNumberOfSurfaces} = vImarisApp.GetFactory.ToSurfaces(vDataItem);
        vNamesList{vNumberOfSurfaces} = char(vDataItem.GetName);
    end
end

if vNumberOfSurfaces==0
    msgbox('Please create at a surfaces object!');
    return;
end

vNamesList = vNamesList(1:vNumberOfSurfaces);
%%
%Choose the surfaces
if vNumberOfSurfaces > 1
    vPair = 1;
    vSurfaces1 = vSurfacesList{vPair(1)}; 
else
    vSurfaces1 = vSurfacesList{1}; 
end
vNumberOfSurfaces1=vSurfaces1.GetNumberOfSurfaces;

%%
% Skip if boutons doesn't exist
if ~vImarisApp.GetFactory.IsSpots(vSpots)  
  continue;
end

% get the spots coordinates
vSpotsXYZ = vSpots.GetPositionsXYZ;
vSpotsTime = vSpots.GetIndicesT;
vSpotsRadius = vSpots.GetRadiiXYZ;
vSpotsName = char(vSpots.GetName);

vSpots.SetVisible(false);
vTimeInterval = min(vSpotsTime):max(vSpotsTime);
vIndicesSpotsTime = cell(numel(vTimeInterval), 1);
for vTime = vTimeInterval
  vIndicesSpotsTime{vTime - vTimeInterval(1) + 1} = find(vSpotsTime == vTime);
end
[vSpotsXYZN,vSpotsD] = size(vSpotsXYZ);
% mask volume (AH Changed to fit one bouton)
vMin = min(vSpotsXYZ);
vMax = max(vSpotsXYZ);
    if vSpotsXYZN == 1;
        vMin = vSpotsXYZ;
        vMax = vSpotsXYZ;
    end
% add 1% border to be sure to include all spots (avoid edge effects)
vDelta = vMax - vMin;
if ~any(vDelta > 0)
  vDelta = [1, 1, 1];
end
vDelta(vDelta == 0) = mean(vDelta(vDelta > 0));
vMin = vMin - vDelta*0.005;
vMax = vMax + vDelta*0.005;

vMaskSize = 350; 

vMaxMaskSize = vMaskSize / max(vMax-vMin);

% spots coordinates on the mask
vSpotsOnMaskXYZ = zeros(size(vSpotsXYZ));
for vDim = 1:3
  vSpotsOnMaskXYZ(:, vDim) = ceil((vSpotsXYZ(:, vDim)-vMin(vDim))*vMaxMaskSize);
end

% the zeros belongs to the first interval
vSpotsOnMaskXYZ = max(vSpotsOnMaskXYZ, 1); 

vMaskSize = ceil((vMax-vMin)*vMaxMaskSize);

% loop through each surface
for vSurfaceIndex = 0:vNumberOfSurfaces1-1
  vAllIndices = vSpotsTime;
  vAllIndicesSize = 0;

  vMask = vSurfaces1.GetSingleMask(vSurfaceIndex, ...
    vMin(1), vMin(2), vMin(3), vMax(1), vMax(2), vMax(3),...
    int32(vMaskSize(1)), int32(vMaskSize(2)), int32(vMaskSize(3)));

  vTimeIndex = vSurfaces1.GetTimeIndex(vSurfaceIndex);
  vMaskImage = vMask.GetDataVolumeAs1DArrayBytes(0, 0);
  vMaskImage = reshape(vMaskImage, vMaskSize);

  % search the element of the spot that lies inside the surface
  vIndexSpotsTime = vIndicesSpotsTime{vTimeIndex - vTimeInterval(1) + 1};
  vSpotsCoords = vSpotsOnMaskXYZ(vIndexSpotsTime, :);
  vIndexSpotsInside = vMaskImage(vSpotsCoords(:, 1) + ...
    (vSpotsCoords(:, 2)-1)*vMaskSize(1) + ...
    (vSpotsCoords(:, 3)-1)*vMaskSize(1)*vMaskSize(2)) == 1;
  vIndexSpotsInside = vIndexSpotsTime(vIndexSpotsInside);

  % copy to complete list
  vSize = numel(vIndexSpotsInside);
  vAllIndices(vAllIndicesSize + (1:vSize)) = vIndexSpotsInside;
  vAllIndicesSize = vAllIndicesSize + vSize;

  vAllIndices = vAllIndices(1:vAllIndicesSize);

  vSpotsInside = vImarisApp.GetFactory.CreateSpots;
  vSpotsInside.Set(vSpotsXYZ(vAllIndices, :), vSpotsTime(vAllIndices), ...
    zeros(sum(vAllIndices~=0),1));
  vSpotsInside.SetRadiiXYZ(vSpotsRadius(vAllIndices,:));
  vSpotsInside.SetName(sprintf('%s inside %s [%i] t%i', ...
    vSpotsName, char(vSurfaces1.GetName), vSurfaceIndex + 1, vTimeIndex));
  vNumberOfSurfaces = max(vNumberOfSurfaces, 2);
  %vRed = (1-vSurfaceIndex/(vNumberOfSurfaces-1)) * 255;
  %vGreen = (1-vChildIndex/(vNumberOfChildren-1)) * 255;
  %vBlue = (vSurfaceIndex/(vNumberOfSurfaces-1)) * 255;
  vSpotsInside.SetColorRGBA((rand(1, 1)) * 256 * 256 * 256 );

  XT = EasyXT(0);
  XT.AddToScene(vSpotsInside);
end
end
%% SPLIT SPOTS 2
XT = EasyXT (0);
bout1 = XT.GetObject('Name', 'GAD65 colocated');
    
for aObjectId = 0;
    vSpots = vImarisApp.GetFactory.ToSpots(bout1);
vSurfaces = vImarisApp.GetFactory.ToSurfaces(vImarisApp.GetSurpassSelection);

vSurfacesSelected = vImarisApp.GetFactory.IsSurfaces(vSurfaces);

if vSurfacesSelected
    vScene = vSurfaces.GetParent;
else
    vScene = vImarisApp.GetSurpassScene;
end
vNumberOfSurfaces = 0;
vSurfacesList{vScene.GetNumberOfChildren} = [];
vNamesList{vScene.GetNumberOfChildren} = [];
for vChildIndex = 1:vScene.GetNumberOfChildren
    vDataItem = vScene.GetChild(vChildIndex - 1);
    if vImarisApp.GetFactory.IsSurfaces(vDataItem)
        vNumberOfSurfaces = vNumberOfSurfaces+1;
        vSurfacesList{vNumberOfSurfaces} = vImarisApp.GetFactory.ToSurfaces(vDataItem);
        vNamesList{vNumberOfSurfaces} = char(vDataItem.GetName);
    end
end

if vNumberOfSurfaces==0
    msgbox('Please create at a surfaces object!');
    return;
end

vNamesList = vNamesList(1:vNumberOfSurfaces);
%%
%Choose the surfaces
if vNumberOfSurfaces > 1
    vPair = 1;
    vSurfaces1 = vSurfacesList{vPair(1)}; 
else
    vSurfaces1 = vSurfacesList{1}; 
end
vNumberOfSurfaces1=vSurfaces1.GetNumberOfSurfaces;

%%
% Skip if boutons doesn't exist
if ~vImarisApp.GetFactory.IsSpots(vSpots)  
  continue;
end

% get the spots coordinates
vSpotsXYZ = vSpots.GetPositionsXYZ;
vSpotsTime = vSpots.GetIndicesT;
vSpotsRadius = vSpots.GetRadiiXYZ;
vSpotsName = char(vSpots.GetName);

vSpots.SetVisible(false);
vTimeInterval = min(vSpotsTime):max(vSpotsTime);
vIndicesSpotsTime = cell(numel(vTimeInterval), 1);
for vTime = vTimeInterval
  vIndicesSpotsTime{vTime - vTimeInterval(1) + 1} = find(vSpotsTime == vTime);
end
[vSpotsXYZN,vSpotsD] = size(vSpotsXYZ);
% mask volume (AH Changed to fit one bouton)
vMin = min(vSpotsXYZ);
vMax = max(vSpotsXYZ);
    if vSpotsXYZN == 1;
        vMin = vSpotsXYZ;
        vMax = vSpotsXYZ;
    end
% add 1% border to be sure to include all spots (avoid edge effects)
vDelta = vMax - vMin;
if ~any(vDelta > 0)
  vDelta = [1, 1, 1];
end
vDelta(vDelta == 0) = mean(vDelta(vDelta > 0));
vMin = vMin - vDelta*0.005;
vMax = vMax + vDelta*0.005;

vMaskSize = 350; 

vMaxMaskSize = vMaskSize / max(vMax-vMin);

% spots coordinates on the mask
vSpotsOnMaskXYZ = zeros(size(vSpotsXYZ));
for vDim = 1:3
  vSpotsOnMaskXYZ(:, vDim) = ceil((vSpotsXYZ(:, vDim)-vMin(vDim))*vMaxMaskSize);
end

% the zeros belongs to the first interval
vSpotsOnMaskXYZ = max(vSpotsOnMaskXYZ, 1); 

vMaskSize = ceil((vMax-vMin)*vMaxMaskSize);

% loop through each surface
for vSurfaceIndex = 0:vNumberOfSurfaces1-1
  vAllIndices = vSpotsTime;
  vAllIndicesSize = 0;

  vMask = vSurfaces1.GetSingleMask(vSurfaceIndex, ...
    vMin(1), vMin(2), vMin(3), vMax(1), vMax(2), vMax(3),...
    int32(vMaskSize(1)), int32(vMaskSize(2)), int32(vMaskSize(3)));

  vTimeIndex = vSurfaces1.GetTimeIndex(vSurfaceIndex);
  vMaskImage = vMask.GetDataVolumeAs1DArrayBytes(0, 0);
  vMaskImage = reshape(vMaskImage, vMaskSize);

  % search the element of the spot that lies inside the surface
  vIndexSpotsTime = vIndicesSpotsTime{vTimeIndex - vTimeInterval(1) + 1};
  vSpotsCoords = vSpotsOnMaskXYZ(vIndexSpotsTime, :);
  vIndexSpotsInside = vMaskImage(vSpotsCoords(:, 1) + ...
    (vSpotsCoords(:, 2)-1)*vMaskSize(1) + ...
    (vSpotsCoords(:, 3)-1)*vMaskSize(1)*vMaskSize(2)) == 1;
  vIndexSpotsInside = vIndexSpotsTime(vIndexSpotsInside);

  % copy to complete list
  vSize = numel(vIndexSpotsInside);
  vAllIndices(vAllIndicesSize + (1:vSize)) = vIndexSpotsInside;
  vAllIndicesSize = vAllIndicesSize + vSize;

  vAllIndices = vAllIndices(1:vAllIndicesSize);

  vSpotsInside = vImarisApp.GetFactory.CreateSpots;
  vSpotsInside.Set(vSpotsXYZ(vAllIndices, :), vSpotsTime(vAllIndices), ...
    zeros(sum(vAllIndices~=0),1));
  vSpotsInside.SetRadiiXYZ(vSpotsRadius(vAllIndices,:));
  vSpotsInside.SetName(sprintf('%s inside %s [%i] t%i', ...
    vSpotsName, char(vSurfaces1.GetName), vSurfaceIndex + 1, vTimeIndex));
  vNumberOfSurfaces = max(vNumberOfSurfaces, 2);
  %vRed = (1-vSurfaceIndex/(vNumberOfSurfaces-1)) * 255;
  %vGreen = (1-vChildIndex/(vNumberOfChildren-1)) * 255;
  %vBlue = (vSurfaceIndex/(vNumberOfSurfaces-1)) * 255;
  vSpotsInside.SetColorRGBA((rand(1, 1)) * 256 * 256 * 256 );

  XT = EasyXT(0);
  XT.AddToScene(vSpotsInside);
end
end
%% Statistics
%% Volume
   XT = EasyXT(0);

   s1 = XT.GetObject('Name', 'mcherry');
   stat1 = xtgetstats(vImarisApp, s1,'All');
   volume = stat1(71).Values;
   
%% GAD65 mcherry
   for ObjectId = 0;
   gad1 = XT.GetObject('Name', 'GAD65 inside mcherry [1] t0');
   if isempty(gad1);
       numbergad1 = 0;
       continue
   end
   statgad1 = xtgetstats(vImarisApp, gad1,'All');
   numbergad1 = statgad1(32).Values;
   end
   
%%  Spots in mCherry
   for ObjectId = 0;
   col1 = XT.GetObject('Name', 'GAD65 colocated inside mcherry [1] t0');
   if isempty(col1);
       insidecol1 = 0;
       continue
   end
   statcol1 = xtgetstats(vImarisApp, col1,'All');
   insidecol1 = statcol1(32).Values;
   end
   %%
   volume = num2cell(volume);
   numbergad1 = num2cell(numbergad1);
   insidecol1 = num2cell(insidecol1);
   data3 = {filename}; 
   alldata = horzcat (data3, volume, numbergad1, insidecol1); 
   A = cell2dataset(alldata, 'VarNames', {'filename', 'mCherry_volume', 'GAD65_in_mCherry', 'GAD65COL_IN_mCherry'}); 
   %%
   Outputname = strcat(filename, '_output.xls'); 
   export (A, 'File', Outputname, 'Delimiter', ',');
   %Eliminates variables for loop
   varlist = {'s1', 'stat1','volume','gad1','statgad1', 'numbergad1','col1', 'statcol1', 'insidecol1', 'data3'};
   clear(varlist{:})

%% Save file              
        channelOutputImage = strcat(filename, '_V1.ims');
        vImarisApp.FileSave(channelOutputImage, 'writer="Imaris5"');
        %vImarisApp.FileClose(filename, '');
    end

end