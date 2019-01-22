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

function MB2v2_BatchSomaticSyt2boutons(aImarisApplicationID)
aImarisApplicationID=0;
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
%% Spots Close to Surface        
        % the user has to create a scene with some spots and surface
        vSurpassScene = vImarisApp.GetSurpassScene;
        if isequal(vSurpassScene, [])
            msgbox('Please create some Spots and Surface in the Surpass scene!')
            return
        end
% get the spots and the surface object

vSpots = vImarisApp.GetFactory.ToSpots(vImarisApp.GetSurpassSelection);
vSurfaces = vImarisApp.GetFactory.ToSurfaces(vImarisApp.GetSurpassSelection);

vSpotsSelected = ~isequal(vSpots, []);
vSurfaceSelected = ~isequal(vSurfaces, []);
if vSpotsSelected
    vParent = vSpots.GetParent;
elseif vSurfaceSelected
    vParent = vSurfaces.GetParent;
else
    vParent = vSurpassScene;
end

% get the spots and surfaces
vSpotsSelection = 1;
vSurfaceSelection = 1;
vNumberOfSpots = 0;
vNumberOfSurfaces = 0;
vSpotsList = [];
vSurfacesList = [];
vSpotsName = {};
vSurfacesName = {};
%Index=Lists of names
for vIndex = 1:vParent.GetNumberOfChildren
    XT = EasyXT (0);
    vItem = vParent.GetChild(vIndex-1);
    if vImarisApp.GetFactory.IsSpots(vItem)
        vNumberOfSpots = vNumberOfSpots + 1;
        vSpotsList(vNumberOfSpots) = vIndex;
        vSpotsName{vNumberOfSpots} = char(vItem.GetName);
        
        if vSpotsSelected && isequal(vItem.GetName, vSpots.GetName)
            vSpotsSelection = vNumberOfSpots; 
        end
    elseif vImarisApp.GetFactory.IsSurfaces(vItem)
        vNumberOfSurfaces = vNumberOfSurfaces + 1;
        vSurfacesList(vNumberOfSurfaces) = vIndex;
        vSurfacesName{vNumberOfSurfaces} = char(vItem.GetName);
        
        if vSurfaceSelected && isequal(vItem.GetName, vSurfaces.GetName)
            vSurfaceSelection = vNumberOfSurfaces;
        end
    end
end

if min(vNumberOfSpots,vNumberOfSurfaces) == 0
    msgbox('Please create some spots AND a surface object!')
    return
end

if vNumberOfSpots>1
    [vSpotsSelection,vOk] = listdlg('ListString',vSpotsName, ...
        'InitialValue', vSpotsSelection, 'SelectionMode','multiple', ...
        'ListSize',[300 300], 'Name','Find Spots Close To Surface', ...
        'PromptString',{'Please select the spots:'});
    if vOk<1, return, end
end
if vNumberOfSurfaces>1
   vSurfaceSelection = [2,3,4,5];
end
%Threshold
vThreshold = 0.4;

vProgressDisplay = waitbar(0,'Finding Spots Close To Surface');

% compute the distances and create new spots objects
vNumberOfSurfacesSelected = numel(vSurfaceSelection);
vNumberOfSpotsSelected = numel(vSpotsSelection);

for vSurfaceIndex = 1:vNumberOfSurfacesSelected
     vItem = vParent.GetChild(vSurfacesList( ...
        vSurfaceSelection(vSurfaceIndex)) - 1);
    vSurface = vImarisApp.GetFactory.ToSurfaces(vItem);

    vSurfaceVertices = [];
    for vIndex = 0:vSurface.GetNumberOfSurfaces - 1
      vSurfaceVertices = [vSurfaceVertices; vSurface.GetVertices(vIndex)];
    end
    vNumberOfVertices = size(vSurfaceVertices, 1);
    
    % limit the memory usage to 3*10000*vNumberOfSpots for each block
    vBlockLimit = 10000;
    vBlockIndices = 1:vBlockLimit:vNumberOfVertices;
    vNumberOfBlock = size(vBlockIndices,2);
    
    for vSpotsIndex = 1:vNumberOfSpotsSelected
        vItem = vParent.GetChild(vSpotsList( ...
            vSpotsSelection(vSpotsIndex)) - 1);
        vSpots = vImarisApp.GetFactory.ToSpots(vItem);
        
		vSpotsPosition = vSpots.GetPositionsXYZ;
		vSpotsTime = vSpots.GetIndicesT;
		vSpotsRadius = vSpots.GetRadiiXYZ;
        vNumberOfSpots = size(vSpotsPosition, 1);
        
        vDistancesMin = [];
        
        for vBlockIndex = 1:vNumberOfBlock
            vBlockMin = vBlockIndices(vBlockIndex);
            vBlockMax = min(vBlockMin+vBlockLimit-1,vNumberOfVertices);
            vBlockSize = vBlockMax - vBlockMin + 1;
            vSparseVertices = 1:vBlockSize;
            % sparse diagonals
            vSurfX = sparse(vSparseVertices,vSparseVertices, double(vSurfaceVertices(vBlockMin:vBlockMax,1)));
            vSurfY = sparse(vSparseVertices,vSparseVertices, double(vSurfaceVertices(vBlockMin:vBlockMax,2)));
            vSurfZ = sparse(vSparseVertices,vSparseVertices, double(vSurfaceVertices(vBlockMin:vBlockMax,3)));
            % build the distance matrix
            vOnes = ones(vBlockSize, vNumberOfSpots);
            vDistX = vOnes*diag(vSpotsPosition(:, 1)) - vSurfX*vOnes;
            vDistY = vOnes*diag(vSpotsPosition(:, 2)) - vSurfY*vOnes;
            vDistZ = vOnes*diag(vSpotsPosition(:, 3)) - vSurfZ*vOnes;
            vDistances = sqrt( min( vDistX.^2+vDistY.^2+vDistZ.^2, [], 1 ) );
            if (vBlockIndex == 1)
                vDistancesMin = vDistances;
            else
                vDistancesMin = min(vDistances,vDistancesMin);
            end
        end
        
        vSpots.SetVisible(false);
        
        vSpotsClose = vDistancesMin <= vThreshold;
        if any(vSpotsClose)
          vNewSpotsClose = vImarisApp.GetFactory.CreateSpots;
		  % The radius can be setted to zero because it is overwritten on the next line
		  % Only an array with the right size is needed
          vNewSpotsClose.Set(vSpotsPosition(vSpotsClose, :), ...
              vSpotsTime(vSpotsClose), zeros(sum(vSpotsClose),1));
          vNewSpotsClose.SetRadiiXYZ(vSpotsRadius(vSpotsClose,:));
          vNewSpotsClose.SetColorRGBA(hex2dec('ff00ff'));
          vNewSpotsClose.SetName(sprintf('%s close to %s [%.2f]', ...
              char(vSpots.GetName), char(vSurface.GetName), vThreshold));
          vParent.AddChild(vNewSpotsClose, -1);
        end
       
        waitbar((vNumberOfSpotsSelected*(vSurfaceIndex-1)+vSpotsIndex) / ...
            (vNumberOfSpotsSelected*vNumberOfSurfacesSelected));
    end
end
close(vProgressDisplay)
clear vSurfacesList
%% SPLIT SPOTS 1
XT = EasyXT (0);
bout1 = XT.GetObject('Name', 'Syt2 close to 1 [0.40]');
    
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
bout1 = XT.GetObject('Name', 'Syt2 close to 2 [0.40]');
    
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
%% SPLIT SPOTS 3
XT = EasyXT (0);
bout1 = XT.GetObject('Name', 'Syt2 close to 3 [0.40]');
    
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
%% SPLIT SPOTS 4
XT = EasyXT (0);
bout1 = XT.GetObject('Name', 'Syt2 close to 4 [0.40]');
    
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
%%Statistics
%% Area
   XT = EasyXT(0);

   s1 = XT.GetObject('Name', '1');
   stat1 = xtgetstats(vImarisApp, s1,'All');
   area(1) = stat1(1).Values;
   
   
  %%Spots close to Surface
   for ObjectId = 0;
   sp1 = XT.GetObject('Name', 'Syt2 close to 1 [0.40]');
   if isempty(sp1);
       numbersp(1) = 0;
       continue
   end
   statp1 = xtgetstats(vImarisApp, sp1,'All');
   numbersp(1) = statp1(32).Values;
   end
   
  %%Spots in mCherry
   for ObjectId = 0;
   sps1 = XT.GetObject('Name', 'Syt2 close to 1 [0.40] inside mcherry [1] t0');
   if isempty(sps1);
       insidesp(1) = 0;
       continue
   end
   statps1 = xtgetstats(vImarisApp, sps1,'All');
   insidesp(1) = statps1(32).Values;
   end
   %%
   area = area';
   numbersp = numbersp';
   insidesp = insidesp';
   
   area = num2cell(area);
   numbersp = num2cell(numbersp);
   insidesp = num2cell(insidesp);
   data = cell(1,1); %change depending on how many cells you have; 
   data = {'cell1'}';
   data3 = cell(1,1); 
   data3(:,:,:,:) = {filename}; 
   alldata = horzcat (data3, data, area, numbersp, insidesp); 
   A = cell2dataset(alldata, 'VarNames', {'filename', 'cell_number', 'area', 'Boutons', 'Boutons_in_mCherry'}); 
   %%
   Outputname = strcat(filename, '_output.xls'); 
   export (A, 'File', Outputname, 'Delimiter', ',');
   %Eliminates variables for loop
   varlist = {'s1','sp1','sps1','stat1','statp1','statps1'};
   varlist2 = {'area', 'numbersp', 'insidesp'};
   clear(varlist{:})
   clear(varlist2{:})
%% Save file              
        channelOutputImage = strcat(filename, '_V1.ims');
        vImarisApp.FileSave(channelOutputImage, 'writer="Imaris5"');
        %vImarisApp.FileClose(filename);
    end

end