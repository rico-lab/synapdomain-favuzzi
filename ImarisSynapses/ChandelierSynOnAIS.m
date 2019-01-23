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

function Emilia_ChSynapses7(aImarisApplicationID)

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
        vImarisApp.FileOpen(filename, ...
                            'reader="Imaris5"' ...
                                );
                            
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
       
 %Localize surface
        %Localize spots
        XT = EasyXT(0);
        channel = 1;
        newSpots = XT.DetectSpots(channel, 'Name', 'mCherry', ...
                                  'Diameter XY', 0.8, ...
                                  'Spots Filter', '"Quality" above 20.21' ...
                         );
                     XT.SetColor(newSpots, [255 0 0]);
                     XT.AddToScene(newSpots);
 %Detect Surface
        XT = EasyXT(0);
        channel = 2;
        surface = XT.DetectSurfaces(channel, 'Name', 'AnkG', ...
                                   'Smoothing', 0.103, ...
                                   'Threshold', 29.44, ...
                                   'Filter', '"Volume" above 3.37 um^3'...
                                   ); 
        XT.SetColor(surface, [0 255 0]);
        XT.AddToScene(surface);
 
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
          vNewSpotsClose.SetName(sprintf('ChSynapses', ...
              char(vSpots.GetName), char(vSurface.GetName), vThreshold));
          vParent.AddChild(vNewSpotsClose, -1);
        end
       
        waitbar((vNumberOfSpotsSelected*(vSurfaceIndex-1)+vSpotsIndex) / ...
            (vNumberOfSpotsSelected*vNumberOfSurfacesSelected));
    end
end
close(vProgressDisplay)
 %%Statistics Spots close to Surface
   for ObjectId = 0;
   sp1 = XT.GetObject('Name', 'kvin');
   if isempty(sp1);
       numbersp(1) = 0;
       continue
   end
   statp1 = xtgetstats(vImarisApp, sp1,'All');
   numbersp(1) = statp1(32).Values;
   end
   numbersp = numbersp';
   numbersp = num2cell(numbersp);
   data3 = cell(1,1); 
   data3(:) = {1};
alldata = horzcat (data3, numbersp); 
   A = cell2dataset(alldata, 'VarNames', {'filename', 'Boutons'}); 
   %%
filename2 = sprintf('%s_%d','filename',i);
xlswrite(filename2,numbersp);
   %Eliminates variables for loop
clear sp1
clear statp1
clear numbersp
   %%
        channelOutputImage = strcat(filename, '_V1.ims');
        vImarisApp.FileSave(channelOutputImage, 'writer="Imaris5"');
    end

end