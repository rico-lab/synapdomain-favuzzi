%
%
%  Batch Process Function for Imaris 7.3.0
%
%  Copyright Bitplane AG 2011
%
%
%  Installation:
%(
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

function MB1v2_BatchSomaticSyt2boutons(aImarisApplicationID)

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
%Localize spots
        XT = EasyXT(0);
        channel = 3;
        newSpots = XT.DetectSpots(channel, 'Name', 'Syt2', ...
                                  'Diameter XY', 0.8, ...
                                  'Spots Filter', '"Quality" above 5.2' ...
                         );
                     XT.SetColor(newSpots, [0 255 0]);
                     XT.AddToScene(newSpots);
         
%Detect Surface
        XT = EasyXT(0);
        channel = 1;
        surface = XT.DetectSurfaces(channel, 'Name', 'mcherry', ...
                                   'Smoothing', 0.103, ...
                                   'Threshold', 15.0, ...
                                  'Filter', '"Number of Voxels" above 100'...
                                   ); 
        XT.SetColor(surface, [255 0 0]);
        XT.AddToScene(surface);
        

   %Detect Surface
        XT = EasyXT(0);
        channel = 2;
        surface = XT.DetectSurfaces(channel, 'Name', 'volume', ...
                                   'Smoothing', 0.103, ...
                                   'Threshold', 100,0 ...
                                  'Filter', '"Volume" above 5.0 um^3'...
                                   ); 
        XT.SetColor(surface, [0 0 255]);
        XT.AddToScene(surface);
        
        %%
        channelOutputImage = strcat(filename, '_V1.ims');
        vImarisApp.FileSave(channelOutputImage, 'writer="Imaris5"');
      
    end

end