function transectViewer(varargin)
%TRANSECTVIEWER - display and edit bathy data
%
%   TRANSECTVIEWER loads and displays all supported
%   raw (*.RAW) files in a specified catalog (.LOG) file.
%   The user can navigate through transects using forward/ back buttons.
%   Filtering options are provided to remove outliers and smooth
%   data.
%
%   FILTERING:
%   Filtering is performed in the following way.  First, the data are
%   broken into chunks that are roughly equal to the value in the "Win.
%   Len." (window length) edit box.  The mean and standard deviation of
%   each chunk is calculated.  Outliers are removed using a threshold
%   number of standard deviations from the mean defined in the "Strength"
%   edit box (ie. the smaller "Strength" leads to more data considered
%   outliers). Once outliers are removed, a running filter is performed.
%   The type and length of the filter are controlled using the "Type"
%   and "Smoothness" control objects.  In addition to automatic outlier
%   rejection, bad data can manually be removed with the "edit manually"
%   button.  In general, the best results are achieved by first removing
%   obviously erroneous data with the manual tool.  Afterwhich, perform the
%   automatic filtering options.  Where bad data have been removed, gaps
%   are filled using linear interpolation.
%
%   WAVE HEIGHT AND PERIOD:
%   You can estimate the wave height and period based on the measured GPS
%   antenna height (Signal Processing Toolbox required). The wave period
%   is relative to the boat (ie- I don't try to adjust for the speed and
%   direction of the boat relative to direction of wave propagation).
%
%   SYNTAX:
%       transectViewer - promts user to supply .LOG file
%       containing names of files to be viewed.
%
%       transectViewer('logname','C:\bathyData\Raw\RAW0420.LOG') - will
%       open
%       the specified log file instead of invoking uigetdir.
%
%
% Andrew Stevens, 5/25/2007
% astevens@usgs.gov

gdata.tv_ver=2.84;
gdata.modified='1/17/2017';


%defaults
logname=[];
ncfile=[];
gdata.invert=1;
%pre-filtering options
gdata.fDist=5; % m
gdata.fStrength=1.5;
%Filter Type
gdata.fType='meanf';
%running filter length
gdata.wLen=10;
%good depths must be between...
gdata.minz=0.2;
gdata.maxz=50;
%maximum offset
gdata.maxoffset=inf;
%maximum interpolation distance
gdata.maxGap=25;
%fix hydrobox echosounder flag (derelict should remove this stuff)
gdata.echofix=0;
gdata.fixlen=200;
%raw data options (odom echosounder)
gdata.rawflag=0;


%parse inputs (if any)
%this is a bit complex, but is easily expandible
if nargin>0
    [m1,n1]=size(varargin); %#ok
    opts={'logname','invert','echofix'};
    
    for i=1:n1;
        indi=strcmpi(varargin{i},opts);
        ind=find(indi==1);
        if isempty(ind)~=1
            switch ind
                case 1
                    logname=varargin{i+1};
                    [pathname,filename] = fileparts(logname);
                    if isempty(pathname)
                        pathname=[pwd,filesep];
                    end
                    if exist(logname,'file')==0
                        errordlg('Log file not found. Try again');
                        return
                    end
                case 2
                    ncfile=varargin{i+1};
                    [pathname,filename] = fileparts(nfile);
                    if isempty(pathname)
                        pathname=[pwd,filesep];
                    end
                    if exist(ncfile,'file')==0
                        errordlg('File not found. Try again');
                        return
                    end
                case 3
                    gdata.invert=0;
                case 4
                    gdata.echofix=1;
            end
        else
        end
    end
end

if (isempty(logname) && isempty(ncfile))
    [filename, pathname, fidx] = uigetfile( ...
        {'*.log', 'LOG Files (*.log)';...
        '*.nc', 'netCDF Files (*.nc)'},...
        'Select a file');
    
end

if filename==0
    return
end

if strcmp(pathname(end),filesep)==0
    pathname=[pathname,filesep];
end

%files in the specified directory
switch fidx
    case 1
        logname=[pathname,filename];
        fnames=textread(logname,'%s'); %#ok
        [~,gdata.logfile] = fileparts(logname);
        d=cellfun(@(x)(dir([pathname,x])),fnames);
        gdata.filesize=arrayfun(@(x)(sprintf('%.1f',x.bytes./1000)),...
            d,'uni',0);
        fInd=1:length(fnames);
        
        %load and show first profile
        % try
        if gdata.invert==1
            [gdata.bdata,gdata.hdr,gdata.rdata]=readRAW(...
                'filename',[pathname,fnames{1}],'invert',...
                'minz',gdata.minz,'maxz',gdata.maxz);
        else
            [gdata.bdata,gdata.hdr,gdata.rdata]=readRAW(...
                'filename',[pathname,fnames{1}],...
                'minz',gdata.minz,'maxz',gdata.maxz);
        end
        
        if isfield(gdata.bdata,'zc')==0
            gdata.bdata.zc=gdata.bdata.zraw;
            gdata.bdata.tide=zeros(size(gdata.bdata.zraw));
            gdata.rdata.tide=zeros(size(gdata.rdata.x));
            
        end
        
    case 2 %netcdf
        ncfile=[pathname,filename];
        ncdata=nc2tv_struct(ncfile);
        fields=fieldnames(ncdata);
        for i=1:length(fields)
            if ~strcmpi(fields{i},'gsetts')
                gdata.(fields{i})=ncdata.(fields{i});
            end
        end
        
end

%set up figure window
set(0,'units','pixels');
ssize=get(0,'screensize');
lp=round(ssize(3)*0.1);
bp=round(ssize(4)*0.05);
wp=round(ssize(3)*0.8);
hp=round(ssize(4)*0.9);

hfig=figure('position',[lp bp wp hp]);
set(hfig,'name','CPS Transect Viewer',...
    'menubar','none','numbertitle','off',...
    'KeyReleaseFcn',@cpfcn,...
    'WindowScrollWheelFcn',@scrollfcn,...
    'renderer','zbuf')
gdata.fh=gcf;




if gdata.echofix==1
    gdata.bdata=echofix(gdata.fixlen,gdata.bdata);
end



figure(gdata.fh)
clf
gdata.ax=axes;
gdata.l1=plot(gdata.bdata.distance,gdata.bdata.zc);
hold on
if isfield(gdata,'cdata')
    gdata.gg=plot(gdata.bdata.distance,...
        gdata.cdata.zc,'r-','linewidth',2);
    gdata.flag=zeros(length(gdata.cdata.zc),1);
end


pos=get(gca,'position');
pos(2)=pos(2)+0.1;
pos(4)=pos(4)-0.2;
set(gca,'position',pos)
gdata.ax1=gca;

gdata.xlimo=[min(gdata.bdata.distance),...
    max(gdata.bdata.distance)];
gdata.ylimo=[min(gdata.bdata.zc),...
    max(gdata.bdata.zc)];
gdata.xlims=[min(gdata.bdata.distance),...
    max(gdata.bdata.distance)];
gdata.ylims=[min(gdata.bdata.zc),...
    max(gdata.bdata.zc)];

goodtrans=1;


%gui elements
%navigation elements
gdata.push1=uicontrol('style','pushbutton',...
    'units','normalized','position',[0.01 0.05 0.085 0.05],...
    'string','<<<','fontsize',14,'callback',@movePhoto);

gdata.push2=uicontrol('style','pushbutton',...
    'units','normalized','position',[0.105 0.05 0.085 0.05],...
    'string','>>>','fontsize',14,'callback',@movePhoto);

gdata.edit1=uicontrol('style','edit',...
    'units','normalized','position',[0.20 0.05 0.085 0.05],...
    'string',num2str(1),'fontsize',14,'callback',@movePhoto);

gdata.text1=uicontrol('style','text',...
    'fontsize',8,'units','normalized','position',...
    [0.01 0.01 0.275 0.03],'horizontalalignment','left');

switch fidx
    case 1
        if goodtrans==1;
            textstr= ['Displaying transect 1 of ',num2str(length(fnames)) ,...
                ' - ',fnames{1}];
            set(gdata.text1,'string',textstr)
            
        else
            textstr= ['Problems loading transect 1 of ',num2str(length(fnames)) ,...
                ' - ',fnames{1}];
            set(gdata.text1,'string',textstr,'foregroundcolor','r')
        end
        gdata.fnames=fnames;
        gdata.fInd=fInd;
        gdata.pointer=1;
    case 2
        set(gdata.push1,'enable','off')
        set(gdata.push2,'enable','off')
        set(gdata.edit1,'enable','off')
        set(gdata.text1,'string',['Displaying ',filename])
end



%filter controls
gdata.hp = uipanel('Title','Filter Controls','FontSize',12,...
    'BackgroundColor','white',...
    'Position',[.295 .01 .45 .1]);

gdata.list1=uicontrol('parent',gdata.hp,'style','popup',...
    'units','normalized','position',[0.01 0.1 0.188 0.45],...
    'string',{'mean';'median';'min';'max'},...
    'callback',@xyzFilt);

gdata.text2=uicontrol('parent',gdata.hp,'Style','text',...
    'units','normalized','position',[0.01 0.65 0.188 0.3],...
    'string','Type');

gdata.edit2=uicontrol('parent',gdata.hp,'Style','edit',...
    'units','normalized','position',[0.208 0.1 0.188 0.45],...
    'string',num2str(gdata.fDist),...
    'callback',@xyzFilt);

gdata.text3=uicontrol('parent',gdata.hp,'Style','text',...
    'units','normalized','position',[0.208 0.65 0.188 0.3],...
    'string','Win. Len.');

gdata.edit3=uicontrol('parent',gdata.hp,'Style','edit',...
    'units','normalized','position',[0.406 0.1 0.188 0.45],...
    'string',num2str(gdata.fStrength),...
    'callback',@xyzFilt);

gdata.text5=uicontrol('parent',gdata.hp,'Style','text',...
    'units','normalized','position',[0.406 0.65 0.188 0.3],...
    'string','Strength');

gdata.edit4=uicontrol('parent',gdata.hp,'Style','edit',...
    'units','normalized','position',[0.604 0.1 0.188 0.45],...
    'string',num2str(gdata.wLen),...
    'callback',@xyzFilt);

gdata.text4=uicontrol('parent',gdata.hp,'Style','text',...
    'units','normalized','position',[0.604 0.65 0.188 0.3],...
    'string','Smoothing');

gdata.push5=uicontrol('parent',gdata.hp,'style','pushbutton',...
    'units','normalized','position',[0.802 0.05 0.188 0.85],...
    'string','Edit Manually','callback',@cleanProfile);

%raw data controls
gdata.blanking=0.5;
gdata.sf1=35000;
gdata.bwin=0.2;

gdata.rawe1=uicontrol('parent',gdata.hp,'Style','edit',...
    'units','normalized','position',[0.01 0.1 0.188 0.45],...
    'string',sprintf('%.1f',gdata.blanking),...
    'visible','off','callback',@digitizebtm);

gdata.rawt1=uicontrol('parent',gdata.hp,'Style','text',...
    'units','normalized','position',[0.01 0.65 0.188 0.3],...
    'string','Blanking Distance (m)','visible','off');

gdata.rawt2=uicontrol('parent',gdata.hp,'Style','text',...
    'units','normalized','position',[0.208 0.65 0.188 0.3],...
    'string','Bottom Threshold','visible','off');

gdata.rawe2=uicontrol('parent',gdata.hp,'Style','edit',...
    'units','normalized','position',[0.208 0.1 0.188 0.45],...
    'string',sprintf('%d',gdata.sf1),'visible','off',...
    'callback',@digitizebtm);

gdata.rawt3=uicontrol('parent',gdata.hp,'Style','text',...
    'units','normalized','position',[0.406 0.65 0.188 0.3],...
    'string','Min. Bottom Width (m)','visible','off');

gdata.rawe3=uicontrol('parent',gdata.hp,'style','edit',...
    'units','normalized','position',[0.406 0.1 0.188 0.45],...
    'string',sprintf('%.2f',gdata.bwin),'visible','off',...
    'callback',@digitizebtm);

gdata.rawp1=uicontrol('parent',gdata.hp,'style','pushbutton',...
    'units','normalized','position',[0.6040 0.05 0.188 0.85],...
    'string','Hand Digitize','callback',@handdigitize,...
    'visible','off');



%view controls
gdata.hp2 = uipanel('Title','View Controls','FontSize',12,...
    'BackgroundColor','white',...
    'Position',[.75 .01 .24 .1]);

gdata.push3=uicontrol('parent',gdata.hp2,'Style','pushbutton',...
    'units','normalized','position',[0.03 0.05 0.3 0.85],...
    'string','zoom','callback',@zoomto);

gdata.toggle1=uicontrol('parent',gdata.hp2,'style','togglebutton',...
    'units','normalized','position',[0.35 0.05 0.3 0.85],...
    'string','Pan','callback',@panIm);

gdata.push4=uicontrol('parent',gdata.hp2,'Style','pushbutton',...
    'units','normalized','position',[0.67 0.05 0.3 0.85],...
    'string','Full Extents','callback',@fullExtents);

%menu controls
gdata.menu1=uimenu('label','File');
gdata.open=uimenu(gdata.menu1,'Label','Open');
uimenu(gdata.open,'label','XYZ File','callback',@addTopo);
uimenu(gdata.open,'label','netCDF','callback',@open_nc)
uimenu(gdata.open,'label','META File','callback',@addMeta);
gdata.logf=uimenu(gdata.menu1,'label','Log File Contents',...
    'callback',@showLog);
if ~isdeployed
    uimenu(gdata.menu1,'Label','Send data to workspace',...
        'callback',@export1);
end

uimenu(gdata.menu1,'label','Load Metadata',...
    'callback',@run_read_info);
uimenu(gdata.menu1,'label','Load Transducer Settings',...
    'callback',@load_xducer);
gdata.ppkmenu=uimenu(gdata.menu1,'label','Post-processed GPS');
uimenu(gdata.ppkmenu,'label','Import',...
    'callback',@run_ppk_gui)
gdata.ppkview=uimenu(gdata.ppkmenu,'label','View PPK GPS Data',...
    'callback',@view_ppk_data,'enable','off');

gdata.menu2=uimenu(gdata.menu1,'Label','Export File');
uimenu(gdata.menu2,'label','To .xyz','callback',@export2)
uimenu(gdata.menu2,'label','To .nc','callback',@export6);
uimenu(gdata.menu2,'Label','To .kml',...
     'callback',@export3);
uimenu(gdata.menu2,'label','To .shp','callback',@export_shp);



uimenu(gdata.menu2,'Label','Batch','callback',@exportAll);

gdata.configexport=uimenu(gdata.menu1,'label','Configure');
uimenu(gdata.configexport,'label','Google Earth Export Options',...
    'callback',@gesetup);
uimenu(gdata.configexport,'label','Batch Export Options',...
    'callback',@tv_batch_gui);



gdata.menu10=uimenu('label','Edit');
gdata.menu11=uimenu(gdata.menu10,'label','Edit Manually',...
    'callback',@cleanProfile);
gdata.lfmenu=uimenu(gdata.menu10,'label','Local Filter');
gdata.lfmenu1=uimenu(gdata.lfmenu,'label','Local Filter Options',...
    'callback',@filt_local);
gdata.lfmenu2=uimenu(gdata.lfmenu,'label','Apply local Filter',...
    'callback',@applylf);
gdata.menu12=uimenu(gdata.menu10,'label','Undo',...
    'visible','off','callback',@undo);
gdata.menu13=uimenu(gdata.menu10,'label','Clear Edits',...
    'visible','off','callback',@clearEdits);


gdata.menu7=uimenu('label','View');
gdata.menu8=uimenu(gdata.menu7,'label','Digitized Data',...
    'callback',@showCleanData,'visible','off');
gdata.menu9=uimenu(gdata.menu7,'label','Raw Data',...
    'callback',@showRawData);
uimenu(gdata.menu7,'label','Line Info','callback',@showLineInfo);
gdata.xdmenu=uimenu(gdata.menu7,'label','Display Transducer Settings',...
    'callback',@showXducerInfo);
set(gdata.xdmenu,'visible','off');
gdata.showmeta=uimenu(gdata.menu7,'label','Display Metadata',...
    'callback',@show_meta_info,'visible','off');
if isfield(gdata,'info')
    set(gdata.showmeta,'visible','on')
end


uimenu(gdata.menu7,'label','Linefile','callback',@showLineFile);
uimenu(gdata.menu7,'label','Offline','callback',@showOffline);
uimenu(gdata.menu7,'label','GPS Info','callback',@showgps);


gdata.menu14=uimenu(gdata.menu7,'label','Along Line Distance',...
    'callback',@viewAlong);
if isempty(gdata.bdata.adist)
    set(gdata.menu14,'visible','off')
end
gdata.menu15=uimenu(gdata.menu7,'label','Absolute Distance',....
    'visible','off','callback',@viewAbsolute);



gdata.menu4=uimenu('label','Utilities');
gdata.wv=uimenu(gdata.menu4,'label','Estimate wave height','callback',@calcHs);
uimenu(gdata.menu4,'label','Speed of Sound Correction',...
    'callback',@run_sos_gui)
uimenu(gdata.menu4,'label','Compare Profiles','callback',...
    @bathy_comp)
gdata.plotrtkppk=uimenu(gdata.menu4,'label','PPK - RTK GPS Comparison',...
    'callback',@plot_ppk_rtk,'visible','off');

bmenu=uimenu(gdata.menu4,'label','NetCDF Batch Processing');

uimenu(bmenu,'label','Add metadata',...
    'callback',@batch_add_meta_gui);
uimenu(bmenu,'label','Add vertical offset',...
    'callback',@batch_add_offset);



gdata.menu5=uimenu('label','Options');
gdata.depthlims=uimenu(gdata.menu5,'label','Define depth limits');
gdata.zlims=uimenu(gdata.depthlims,'label','Digitized Data','callback',@zlims);
gdata.chrawz=uimenu(gdata.depthlims,'label','Raw Data','callback',@setmaxraw);
gdata.set_file_date=uimenu(gdata.menu5,'label','Define File Date',...
    'callback',@set_file_date);
gdata.off1=uimenu(gdata.menu5,'label','Apply Offset');
gdata.interpDist=uimenu(gdata.menu5,'label','Max Interpolation Distance',...
    'callback',@interpDist);
gdata.manOffset=uimenu(gdata.off1,'label','Manual Offset',...
    'callback',@manualOffset);
uimenu(gdata.menu5,'label','Tide Options','callback',@run_tidedlg);
if gdata.echofix==1;
    gdata.efmenu=uimenu(gdata.menu5,'label','Hyrdobox Fix Window Length',...
        'callback',@runEchoFix);
end

switch fidx
    case 1
        if exist([pathname,strtok(fnames{1},'.'),'.bin'],'file')==0;
            set(gdata.menu9,'visible','off')
            set(gdata.chrawz,'visible','off')
        end
        gdata.bindata=[];
    case 2
        if ~isfield(gdata,'bindata')
            set(gdata.menu9,'visible','off')
            set(gdata.chrawz,'visible','off')
        end
        set(gdata.logf,'visible','off')
end

a=ver;
tboxes={a.Name}';
if any(strcmpi('signal processing toolbox',tboxes))~=1;
    set(gdata.wv,'enable','off')
end
gdata.menu16=uimenu(gdata.menu5,'label','Raw Data Color Limits',...
    'callback',@setrclims,'visible','off');

gdata.apply_ppk_menu=uimenu(gdata.menu5,'label',...
    'Apply Post-processed GPS (PPK)',...
    'Callback',@apply_ppk_menu,...
    'enable','off');



gdata.menu6=uimenu('label','About');
uimenu(gdata.menu6,'label','About','callback',@dispHelp);
uimenu(gdata.menu6,'label','Keyboard Shortcuts',...
    'callback',@dispKeys)



if isfield(gdata,'cdata')
    
    switch ncdata.gsetts.smoothing_type
        case 'mean'
            gdata.fType='mean';
            set(gdata.list1,'value',1)
        case 'median'
            gdata.fType='median';
            set(gdata.list1,'value',2)
        case 'min'
            gdata.fType='min';
            set(gdata.list1,'value',3)
        case 'max'
            gdata.fType='max';
            set(gdata.list1,'value',4)
    end
    set(gdata.edit2,'string',...
        sprintf('%0.0f',ncdata.gsetts.smoothing_length));
    set(gdata.edit3,'string',...
        sprintf('%0.0f',ncdata.gsetts.fstrength))
    set(gdata.edit4,'string',...
        sprintf('%0.0f',ncdata.gsetts.window_length))
    gdata.maxGap=ncdata.gsetts.max_interp_dist;
else
    gdata.cdata=[];
end

gdata.numedits=0;
gdata.offFig=[];
gdata.tcorr=[];
if ~isfield(gdata,'manoff')
    gdata.manoff=0;
end
gdata.pan=0;
gdata.newInd=1;
gdata.logname=filename;
gdata.filepath=pathname;

gdata.outpath=gdata.filepath;
gdata.outpath1=gdata.filepath;
gdata.topo=[];
gdata.topopath=gdata.filepath;
gdata.metapath=gdata.filepath;

gdata.maxrawz=-20;

gdata.ja = java.awt.Robot;

gdata.alongflag=0;
gdata.tlineh=[];

if ~isfield(gdata,'flag')
    gdata.flag=[];
end

gdata.ge.thin=20;
gdata.ge.cmin=-10;
gdata.ge.cmax=3;
gdata.ge.scale=0.3;
gdata.ge.cmap='jet';

if ~isfield(gdata,'info')
    gdata.info=[];
end
gdata.rclims=[];

gdata.lfd.lftype=1;
gdata.lfd.lflen=3;

gdata.applyppk=0;

gdata.batch.out_xyz=1;
gdata.batch.out_nc=1;
gdata.batch.out_kml=1;
gdata.batch.out_shp=1;

gdata.tideopt.method='none';
gdata.tideopt.maxgap=gdata.maxGap;
gdata.tideopt.badtide=[];

guidata(hfig,gdata);

end

%%%%Callbacks--------------------------------------------------------------
%netcdf functionality

function open_nc(hfig,evnt) %#ok

gdata=guidata(hfig);
gdata.cdata=[];
[filename, pathname] = uigetfile( ...
    {'*.nc', 'netCDF Files (*.nc)'},...
    'Select a file',gdata.filepath);
if filename==0
    return
end

ncfile=[pathname,filename];
ncdata=nc2tv_struct(ncfile);
fields=fieldnames(ncdata);
for i=1:length(fields)
    if ~strcmpi(fields{i},'gsetts')
        gdata.(fields{i})=ncdata.(fields{i});
    end
end

%deal with tide options
gdata.tideopt.badtide=[];
guidata(hfig,gdata);
applytidecorr(hfig);
gdata=guidata(hfig);

hold off
if isfield(gdata.bdata,'zc')==0
    gdata.bdata.zc=gdata.bdata.zraw;
    gdata.bdata.tide=zeros(size(gdata.bdata.zraw));
    gdata.rdata.tide=zeros(size(gdata.rdata.antennaH));
end


% if isfield(gdata,'sv');
%     if gdata.invert==1
%         zraw=-gdata.bdata.zraw;
%     end
%     zsos=apply_sos_prof(zraw,gdata.sv.sos_orig,...
%         [gdata.sv.depth gdata.sv.sos],gdata.sv.use_mean_sos,...
%         gdata.sv.mean_vel)-gdata.manoff;
%     if gdata.invert==1
%         zsos=-zsos;
%     end
%     
%     gdata.bdata.zc=zsos-gdata.bdata.tide;
% else
%     gdata.bdata.zc=(gdata.bdata.zraw+gdata.manoff)-...
%         gdata.bdata.tide;
% end
% if isfield(gdata,'sos_corr');
%     zsos=(gdata.bdata.zraw.*gdata.sos_corr.ratio)-gdata.manoff;
%     gdata.bdata.zc=zsos-gdata.bdata.tide;
% else
%     gdata.bdata.zc=(gdata.bdata.zraw-gdata.manoff)-...
%         gdata.bdata.tide;
% end

%what about using ppk GPS
% if gdata.applyppk==2
%     guidata(hfig,gdata)
%     applyppk(hfig)
%     gdata=guidata(hfig);
% else
%     if isfield(gdata.info,'ppk_filename')
%         gdata.info=rmfield(gdata.info,'ppk_filename');
%     end

gdata.applyppk=0;


if isfield(gdata.bdata,'tide')
    set(gdata.wv,'enable','on')
    
else
    set(gdata.wv,'enable','off')
    
end

if isfield(gdata,'gg')
    gdata=rmfield(gdata,'gg');
end

figure(gdata.fh)
hold off
gdata.l1=plot(gdata.bdata.distance,gdata.bdata.zc);
hold on
if ~isempty(gdata.cdata)
    gdata.gg=plot(gdata.bdata.distance,...
        gdata.cdata.zc,'r-','linewidth',2);
    gdata.flag=zeros(length(gdata.cdata.zc),1);
else
    gdata.flag=[];
end

fun=@(x)(sort([min(x)-(0.05*(max(x)-min(x))),...
    max(x)+(0.05*(max(x)-min(x)))]));
zidx=isfinite(gdata.bdata.zraw);
gdata.xlimo=fun(gdata.bdata.distance(zidx));
gdata.ylimo=fun(gdata.bdata.zc);
gdata.xlims=fun(gdata.bdata.distance(zidx));
gdata.ylims=fun(gdata.bdata.zc);
set(gca,'xlim',gdata.xlimo,...
    'ylim',gdata.ylimo)


if ~isempty(gdata.cdata)
    
    switch ncdata.gsetts.smoothing_type
        case 'mean'
            gdata.fType='mean';
            set(gdata.list1,'value',1)
        case 'median'
            gdata.fType='median';
            set(gdata.list1,'value',2)
        case 'min'
            gdata.fType='min';
            set(gdata.list1,'value',3)
        case 'max'
            gdata.fType='max';
            set(gdata.list1,'value',4)
    end
    set(gdata.edit2,'string',...
        sprintf('%0.0f',ncdata.gsetts.smoothing_length));
    set(gdata.edit3,'string',...
        sprintf('%0.0f',ncdata.gsetts.fstrength))
    set(gdata.edit4,'string',...
        sprintf('%0.0f',ncdata.gsetts.window_length))
    gdata.maxGap=ncdata.gsetts.max_interp_dist;
end



set(gdata.push1,'enable','off')
set(gdata.push2,'enable','off')
set(gdata.edit1,'enable','off')
set(gdata.text1,'string',['Displaying ',filename],...
    'foregroundcolor','k')


if isfield(gdata,'bindata')
    if ~isempty(gdata.bindata)
        if ~isfield(gdata,'menu7')
            gdata.menu7=uimenu('label','View');
            gdata.menu8=uimenu(gdata.menu7,'label','Digitized Data',...
                'callback',@showCleanData,'visible','off');
            gdata.menu9=uimenu(gdata.menu7,'label','Raw Data',...
                'callback',@showRawData);
        else
            set(gdata.menu7,'visible','on');
            set(gdata.menu8,'visible','off');
            set(gdata.menu9,'visible','on');
        end
        set(gdata.chrawz,'visible','on')
    end
else
    gdata.bindata=[];
    set(gdata.chrawz,'visible','off')
    if isfield(gdata,'menu7')
        set(gdata.menu9,'visible','off');
    end
    
end

if isfield(gdata,'info')
    set(gdata.showmeta,'visible','on')
end

if isfield(gdata,'eraw')
    gdata=rmfield(gdata,{'eraw';'ezc'});
end

gdata.alongflag=0;
gdata.rawflag=0;
set(gdata.hp,'Title','Filter Controls');
set(gdata.list1,'visible','on');
set(gdata.text2,'visible','on');
set(gdata.edit2,'visible','on');
set(gdata.text3,'visible','on');
set(gdata.edit3,'visible','on');
set(gdata.text5,'visible','on');
set(gdata.edit4,'visible','on');
set(gdata.text4,'visible','on');
set(gdata.push5,'callback',@cleanProfile);
set(gdata.menu11,'callback',@cleanProfile);
set(gdata.logf,'visible','off')
set(gdata.rawe1,'visible','off')
set(gdata.rawt1,'visible','off')
set(gdata.rawe2,'visible','off')
set(gdata.rawt2,'visible','off')
set(gdata.rawe3,'visible','off')
set(gdata.rawt3,'visible','off')
set(gdata.rawp1,'visible','off')

set(gdata.menu15,'visible','off');
if ~isempty(gdata.bdata.adist)
    set(gdata.menu14,'visible','on')
end

gdata.numedits=0;
gdata.edits=[];
set(gdata.menu12,'visible','off')
set(gdata.menu13,'visible','off')
gdata.topo=[];
if isfield(gdata,'bb')
    gdata=rmfield(gdata,'bb');
end

gdata.filepath=pathname;
guidata(hfig,gdata)
setFocus(hfig);

end

%%%%------------------------------------------------------------------

function varargout=raw2nc(varargin)
% RAW2NC - convert Hypack .Raw file to netCDF

p=inputParser;
opts={'filename',   [],     {'char'},    {};...
    'fpath',        [],     {'char'},    {};...
    'outfile',      [],     {'char'},    {};...
    'bdata',        [],     {'struct'},  {};...
    'fdata',        [],     {'struct'},  {};...
    'rdata',        [],     {'struct'},  {};...
    'hdr',          [],     {'struct'},  {};...
    'bindata',      [],     {'struct'},  {};...
    'readbin',      1,      {'numeric'}, {};...
    'info',         [],     {'struct'},  {};...
    'sv',           []      {'struct'},  {}};

cellfun(@(x)(p.addParamValue(x{1},x{2},...
    @(y)(validateattributes(y, x{3},x{4})))),num2cell(opts,2));

p.KeepUnmatched = true;
p.parse(varargin{:})
opt=p.Results;

%process the inputs
%if bdata,rdata,hdr,bindata absent, read from file

if all([isempty(opt.filename);...
        isempty(opt.bdata);...
        isempty(opt.hdr);...
        isempty(opt.rdata)])
    [fname, dname] = uigetfile( ...
        '*.Raw', 'Raw Files (*.RAW)',...
        'Select a  file');
    if fname==0
        varargout{1}=[];
        return
    else
        opt.filename=[dname,fname];
    end
else
    if isempty(opt.filename)
        opt.filename=opt.bdata.filename;
    end
    if ~isempty(opt.fpath)
        opt.filename=[opt.fpath,opt.filename];
    end
end

[dname,fname]=fileparts(opt.filename);

%look for bin file
if opt.readbin
    if isempty(opt.bindata)
        binfile=[dname,filesep,fname,'.bin'];
        if ~exist(binfile,'file')
            opt.readbin=0;
            warning('Raw data (*.bin file) not found.')
        else
            opt.bindata=readBinSolo(binfile);
            opt.bindata.range=-opt.bindata.range;
        end
    end
end

%read the raw file,if not supplied
if any([isempty(opt.bdata);...
        isempty(opt.hdr);...
        isempty(opt.rdata)])
    [opt.bdata,opt.hdr,opt.rdata]=readRAW('filename',opt.filename,...
        'invert');
end

%output file
if isempty(opt.outfile)
    [fout,dout]=uiputfile('*.nc','Save As');
    
    if fout==0
        return
    else
        opt.outfile=[dout,fout];
    end
end

%create nc
ncid = netcdf.create(opt.outfile,'NETCDF4');


%write global atrributes
gp=netcdf.getConstant('NC_GLOBAL');

%first the atts supplies by the user
hdr_coord=1; %show coord sys info in hypack hdr
if ~isempty(opt.info);
    fields=fieldnames(opt.info);
    for i=1:length(fields)
        netcdf.putAtt(ncid,gp,fields{i},...
            opt.info.(fields{i}));
    end
    if isfield(opt.info,'ppk_filename');
        [~,~,ext]=fileparts(opt.info.ppk_filename);
        if strcmpi(ext,'.txt')
            hdr_coord=0;
        end
    end
end

%now the stuff in the hypack header
%hypack version
netcdf.putAtt(ncid,gp,...
    'hypack_version',opt.hdr.hypack_ver);
%raw filename
netcdf.putAtt(ncid,gp,...
    'raw_filename',opt.bdata.filename);
%Line Number
netcdf.putAtt(ncid,gp,...
    'line_number',opt.hdr.lineNum);
%Line Coordinates
netcdf.putAtt(ncid,gp,...
    'line_x',opt.hdr.linex);
netcdf.putAtt(ncid,gp,...
    'line_y',opt.hdr.liney);
%coordinate system
if hdr_coord %hypack coord system info not necessarily correct if ppk is applied
    netcdf.putAtt(ncid,gp,...
        'ellipsoid',opt.hdr.ellipsoid{1});
    netcdf.putAtt(ncid,gp,...
        'semimajor_axis',str2double(opt.hdr.ellipsoid{2}));
    netcdf.putAtt(ncid,gp,...
        'inverse_flattening',str2double(opt.hdr.ellipsoid{3}));
    netcdf.putAtt(ncid,gp,...
        'projection',opt.hdr.projection{1});
    netcdf.putAtt(ncid,gp,...
        'reference_lon',str2double(opt.hdr.projection{2}));
    netcdf.putAtt(ncid,gp,...
        'reference_lat',str2double(opt.hdr.projection{4}));
    netcdf.putAtt(ncid,gp,...
        'scale_factor',str2double(opt.hdr.projection{3}));
    netcdf.putAtt(ncid,gp,...
        'north_parallel',str2double(opt.hdr.projection{5}));
    netcdf.putAtt(ncid,gp,...
        'south_parallel',str2double(opt.hdr.projection{6}));
    netcdf.putAtt(ncid,gp,...
        'false_easting',str2double(opt.hdr.projection{7}));
    netcdf.putAtt(ncid,gp,...
        'false_northing',str2double(opt.hdr.projection{8}));
    if isfield(opt.hdr,'geoid')
        netcdf.putAtt(ncid,gp,...
            'geoid',opt.hdr.geoid);
    end
    netcdf.putAtt(ncid,gp,...
        'orthometric_height_corr',...
        str2double(opt.hdr.orthometric_height_corr));
end

if isfield(opt.rdata,'antenna_offset')
    offset=unique(opt.rdata.antenna_offset);
    offset=offset(isfinite(offset));
    netcdf.putAtt(ncid,gp,...
        'gps_antenna_offset',offset);
end


%need to add other Global Atts, sounder settings, benchinfo, etc------

%add groups for gps and sonar data
%start with raw gps (rdata)
gps=netcdf.defGrp(ncid,'gps');

gpsdim=netcdf.defDim(gps,'time',numel(opt.rdata.hygtime));
gfields={'hygtime','hypack_time','hypack time','seconds past midnight';...
    'ptime','mtime','matlab date number','days past midnight Jan 1, 0000';...
    'x','x','Easting','meters';...
    'y','y','Northing','meters';...
    'hdop','hdop','GPS dilution of precision','-';...
    'numsats','numsats','number of satellites','-';...
    'gpsmode','mode','Mode of GPS Solution','-';...
    'antennaH','height','Antenna height','meters';...
    'gpsTime','time','GPS Time','HHMMSS.FF';...
    'lat','latitude','Latitude','decimal degrees north';...
    'lon','longitude','Longitude','decimal degrees east';...
    'undulation','undulation','Geoid Height','meters';...
    'tide','tide','Tide Correction','meters'};
%def. variables
for i=1:length(gfields)
    if isfield(opt.rdata,gfields{i,1})
        varid=netcdf.defVar(gps,gfields{i,2},'double',gpsdim);
        netcdf.putAtt(gps,varid,'long_name',gfields{i,3});
        netcdf.putAtt(gps,varid,'units',gfields{i,4});
    end
end

%sonar data
if ~isempty(opt.bindata)
    [junk,ia,ib] = intersect(opt.bdata.vtime,...
        opt.bindata.vtime);%#ok
    
    
    opt.bdata.mtime=opt.bdata.mtime(ia);
    opt.bdata.vtime=opt.bdata.vtime(ia);
    opt.bdata.x=opt.bdata.x(ia);
    opt.bdata.y=opt.bdata.y(ia);
    opt.bdata.distance=opt.bdata.distance(ia);
    if ~isempty(opt.bdata.adist)
        opt.bdata.adist=opt.bdata.adist(ia);
        opt.bdata.offline=opt.bdata.offline(ia);
    end
    if isfield(opt.bdata,'lon')
        opt.bdata.lon=opt.bdata.lon(ia);
        opt.bdata.lat=opt.bdata.lat(ia);
    end
    opt.bdata.zraw=opt.bdata.zraw(ia);
    opt.bdata.tide=opt.bdata.tide(ia);
    opt.bdata.zc=opt.bdata.zc(ia);
    
    if ~isempty(opt.fdata)
        opt.fdata.zf=opt.fdata.zf(ia);
    end
    
    opt.bindata.vtime=opt.bindata.vtime(ib);
    opt.bindata.vals=opt.bindata.vals(:,ib);
end

sonar=netcdf.defGrp(ncid,'sonar');
sonardim=netcdf.defDim(sonar,'time',numel(opt.bdata.mtime));
sfields={'mtime','mtime','matlab date number','days past midnight Jan 1, 0000';...
    'vtime','hypack_time','hypack time','seconds past midnight';...
    'zraw','zraw','Depth Below Transducer','meters'};

%def. variables
for i=1:size(sfields,1)
    if isfield(opt.bdata,sfields{i,1})
        varid=netcdf.defVar(sonar,sfields{i,2},'double',sonardim);
        netcdf.putAtt(sonar,varid,'long_name',sfields{i,3});
        netcdf.putAtt(sonar,varid,'units',sfields{i,4});
    end
end
%raw data fields
if ~isempty(opt.bindata)
    %define the raw data dimension
    depthdim=netcdf.defDim(sonar,'depth',numel(opt.bindata.range));
    varid=netcdf.defVar(sonar,'ping_num','double',...
        sonardim);
    netcdf.putAtt(sonar,varid,'long_name','Ping Number');
    netcdf.putAtt(sonar,varid,'units','-')
    
    varid=netcdf.defVar(sonar,'amplitude','NC_FLOAT',...
        [depthdim sonardim]);
    netcdf.putAtt(sonar,varid,'long_name','Backscatter Amplitude');
    netcdf.putAtt(sonar,varid,'units','-')
    
    varid=netcdf.defVar(sonar,'range','double',depthdim);
    netcdf.putAtt(sonar,varid,'long_name','Range to Transducer');
    netcdf.putAtt(sonar,varid,'units','meters')
end

%merged dataset
bathy=netcdf.defGrp(ncid,'bathy');
bathydim=netcdf.defDim(bathy,'time',numel(opt.bdata.mtime));
bfields={'mtime','mtime','matlab date number','days past midnight Jan 1, 0000';...
    'vtime','hypack_time','hypack time','seconds past midnight';...
    'x','x','Easting','meters';...
    'y','y','Northing','meters';...
    'distance','distance','Distance','meters';...
    'adist','x_dist','Along-line Distance','meters';...
    'offline','offline','Distance off line','meters';...
    'lat','latitude','Latitude','decimal degrees north';...
    'lon','longitude','Longitude','decimal degrees east';...
    'zraw','zraw','Depth Below Transducer','meters';...
    'tide','tide','Tide Correction','meters';...
    'zc','zc','Corrected Depth','meters'};

if ~isempty(opt.fdata)
    bfields2={'zf','zf','Smoothed, Corrected Depth','meters'};
    bfields=cat(1,bfields,bfields2);
    zfatts={'smoothing_type','ftype';...
        'smoothing_length','fdist';...
        'window_length','wlen';...
        'fstrength','fstrength';...
        'max_interp_dist','maxgap'};
end


%def. variables
for i=1:length(bfields)
    if isfield(opt.bdata,bfields{i,1})
        varid=netcdf.defVar(bathy,bfields{i,2},'double',bathydim);
        netcdf.putAtt(bathy,varid,'long_name',bfields{i,3});
        netcdf.putAtt(bathy,varid,'units',bfields{i,4});
    end
    if strcmpi('zf',bfields{i,1})
        varid=netcdf.defVar(bathy,bfields{i,2},'double',bathydim);
        netcdf.putAtt(bathy,varid,'long_name',bfields{i,3});
        netcdf.putAtt(bathy,varid,'units',bfields{i,4});
        for j=1:length(zfatts)
            netcdf.putAtt(bathy,varid,zfatts{j,1},...
                opt.fdata.(zfatts{j,2}));
        end
    end
end

if ~isempty(opt.sv) %if svel cast exists
    svel=netcdf.defGrp(ncid,'svel');
    if isempty(opt.sv.depth)
        sveldim=netcdf.defDim(svel,'cast_depth',1);
        varid=netcdf.defVar(svel,'svel','double',sveldim);
        netcdf.putAtt(svel,varid,'long_name','Sound Velocity');
        netcdf.putAtt(svel,varid,'units','m/s')
    else
        
        sveldim=netcdf.defDim(svel,'cast_depth',numel(opt.sv.depth));
        
        svfields={'depth','depth','cast water depth','meters';...
            'sos','svel','sound velocity','m/s'};
        %to do - add cast coordinates, not sure if using avg of mutilple
        %casts
        
        for i=1:size(svfields,1)
            if isfield(opt.sv,svfields{i,1})
                varid=netcdf.defVar(svel,svfields{i,2},'double',sveldim);
                netcdf.putAtt(svel,varid,'long_name',svfields{i,3});
                netcdf.putAtt(svel,varid,'units',svfields{i,4});
            end
        end
    end
end

netcdf.endDef(ncid);

%write gps data
[~,nvars]=netcdf.inq(gps);
for i=1:nvars
    varname=netcdf.inqVar(gps,i-1);
    midx=find(strcmpi(varname,gfields(:,2)));
    netcdf.putVar(gps,i-1,opt.rdata.(gfields{midx,1})); %#ok
end

%write sonar data
[~,nvars]=netcdf.inq(sonar);
for i=1:nvars
    varname=netcdf.inqVar(sonar,i-1);
    midx=find(strcmpi(varname,sfields(:,2)));
    if ~isempty(midx)
        netcdf.putVar(sonar,i-1,opt.bdata.(sfields{midx,1}))
    end
end
if ~isempty(opt.bindata)
    if numel(opt.bindata.vtime)>numel(opt.bdata.mtime);
        
        ind=find(opt.bindata.vtime>=min(opt.bdata.vtime) & ...
            opt.bindata.vtime<=max(opt.bdata.vtime));
    else
        ind=1:numel(opt.bindata.vtime);
    end
    
    varid=netcdf.inqVarID(sonar,'ping_num');
    netcdf.putVar(sonar,varid,(1:numel(ind))')
    
    varid=netcdf.inqVarID(sonar,'amplitude');
    netcdf.putVar(sonar,varid,opt.bindata.vals(:,ind))
    
    varid=netcdf.inqVarID(sonar,'range');
    netcdf.putVar(sonar,varid,opt.bindata.range)
end





%write bathy data
[~,nvars]=netcdf.inq(bathy);
for i=1:nvars
    varname=netcdf.inqVar(bathy,i-1);
    midx=find(strcmpi(varname,bfields(:,2)));
    if strcmpi('zf',varname)
        netcdf.putVar(bathy,i-1,opt.fdata.zf);
    else
        if isempty(opt.bdata.(bfields{midx,1}));
            opt.bdata.(bfields{midx,1})=opt.bdata.mtime.*NaN;
        end
        netcdf.putVar(bathy,i-1,opt.bdata.(bfields{midx,1}));
    end
end



%write svel data
if ~isempty(opt.sv)
    if ~isempty(opt.sv.depth)
        [~,nvars]=netcdf.inq(svel);
        for i=1:nvars
            varname=netcdf.inqVar(svel,i-1);
            midx=find(strcmpi(varname,svfields(:,2)));
            netcdf.putVar(svel,i-1,opt.sv.(svfields{midx,1})); %#ok
        end
    else
        netcdf.putVar(svel,0,opt.sv.mean_vel)
        
    end
end


netcdf.close(ncid);

if nargout>0
    varargout{1}=rawnc2mat(opt.outfile);
    if nargout==2
        varargout{2}=opt;
    end
end
end

%%%%---------------------------------------------------------------------
function data=rawnc2mat(filename)
% RAWNC2MAT- Read netCDF Hypack file

%read the global attributes
info=ncinfo(filename);
ngatts = length(info.Attributes);
gnames=genvarname(arrayfun(@(x)(x.Name),info.Attributes,'un',0)');
gvals=arrayfun(@(x)(x.Value),info.Attributes,'un',0)';
gatts=cell2struct(gvals,gnames);

ncid=netcdf.open(filename);
grps=netcdf.inqGrps(ncid);
grp_names=cellfun(@(x)(netcdf.inqGrpName(x)),...
    num2cell(grps'),'un',0);

%read the group data
for i=1:length(grps);
    [ndims,nvars]=netcdf.inq(grps(i));
    [vnames,~,~,natts]=cellfun(@(x)(netcdf.inqVar(...
        grps(i),x)),num2cell((0:nvars-1)'),'un',0);
    vdata=cellfun(@(x)(netcdf.getVar(grps(i),x)),...
        num2cell((0:nvars-1)'),'un',0);
    data.(grp_names{i})=cell2struct(vdata,vnames);
    
    for j=1:nvars
        attnames=cellfun(@(x)(netcdf.inqAttName(grps(i),...
            j-1,x)),num2cell((0:natts{j}-1)'),'un',0);
        attvals=cellfun(@(x)(netcdf.getAtt(grps(i),...
            j-1,x)),attnames,'un',0);
        data.(grp_names{i}).vatts.(vnames{j})=...
            cell2struct(attvals,attnames);
    end
    
    
end

data.gatts=gatts;
netcdf.close(ncid);

end

%%%%-----------------------------------------------------------------------
function bindata = readBinSolo(fname)
% READBIN - Read raw acoustic data from .BIN file
%
%   READBIN reads binary files with the extention .BIN
%   containing raw waveform acoustic data from a single-
%   beam echosounder recorded with HYPACK. A .RAW file with
%   the same name must be present in the same directory
%   as the .BIN file. The code is currently only configured
%   for a single-frequency sonar system.
%
%   INPUTS: Two optional input arguements are available
%       filename - If a string is among the input arguments,
%                  the program will use this as the file to
%                   process.  If no filename is provided,
%                   the user will be prompted to select a file.
%       plotflag - 0 (default) or 1. If plotflag is set to 1,
%                  a simple plot will be made.
%
%   OUTPUT: The data are returned in a structure with
%   the following fileds:
%       'filename' -  name of data file
%       'vtime'    -  Hypack time tag (millisecs since midnight)
%       'range'    -  range from transducer (m)
%       'vals'     -  backscatter values (unknown units)
%       'mtime'    -  Matlab datenum
%       'lat'      -  latitude (dec. degrees, WGS84)
%       'lon'      -  longitude (dec. degrees, WGS84)
%       'x'        -  projected x-coordinate
%       'y'        -  projected y-coordinate
%       'tide'     -  tide correction (m)
%       'zraw'     -  digitized depth, uncorrected
%       'zc'       -  digitized depth, corrected
%
%   EXAMPLES:
%
%   data=readBIN; %user selects file, no plot made
%
%   data=readBIN('d:\data\foo.bin'); %user supplies file, no plot
%
%   data=readBIN(1); %user selects file, plot is made
%
% SEE ALSO readRAW transectViewer

% Andrew Stevens
% astevens@usgs.gov
% 8/12/2009

%check inputs


% h=waitbar(0);

%
% set(h,'name','Reading .BIN file');

fid=fopen(fname,'r');
fseek(fid,0,'eof');
numbytes=ftell(fid);
frewind(fid)

pos=ftell(fid);
numread=1;
while pos<numbytes
    
    dummy=fread(fid,1,'uint32'); %#ok
    vers=fread(fid,1,'uint8'); %#ok
    dep_no=fread(fid,1,'uint8'); %#ok
    timer=fread(fid,1,'uint32')./1000;
    draft=fread(fid,1,'ushort'); %#ok
    unit=fread(fid,1,'char'); %#ok
    sw(numread)=fread(fid,1,'ushort'); %#ok
    eos(numread)=fread(fid,1,'ushort'); %#ok
    num_samp=fread(fid,1,'ushort');
    if isempty(num_samp)
        break
    end
    fseek(fid,16,'cof');
    samp=fread(fid,num_samp,'ushort');
    
    if numread==1
        bytes_per_ping=35+(num_samp*2);
        num_pings=floor(numbytes/bytes_per_ping);
        
        bindata=struct('filename',fname,'vtime',zeros(1,num_pings),...
            'range',zeros(num_samp,num_pings),...
            'vals',zeros(num_samp,num_pings));
        
    end
    
    
    bindata.range(:,numread)=linspace(eos(numread)-sw(numread),...
        eos(numread),num_samp);
    
    if isempty(samp), break, end
    
    bindata.vtime(numread)=timer;
    bindata.vals(:,numread)=samp;
    
    
    numread=numread+1;
    pos=ftell(fid);
    
    %     waitbar(pos/numbytes,h,sprintf('%d%% complete...',...
    %         round((pos/numbytes)*100)));
    
end
fclose(fid);


if (numel(unique(sw))>1 || numel(unique(eos))>1);
    %     set(h,'name','Processing Bin File.');
    rangei=linspace(min(eos)-max(sw),...
        max(eos),num_samp);
    
    for i = 1:numread-1;
        bindata.vals(:,i)=interp1(bindata.range(:,i),...
            bindata.vals(:,i),rangei);
        
        %         waitbar(i/(numread-1),h,sprintf('%d%% complete...',...
        %                 round((i/(numread-1))*100)));
    end
    bindata.range=rangei;
else
    bindata.range=bindata.range(:,1);
end
% close(h);
end
%%%%----------------------------------------------------------------------
function gdata=nc2tv_struct(ncfile)
%this function imports data from netcdf into data structures
%consistent with transectviewer.

ncdata=rawnc2mat(ncfile);

%need to convert fieldnames
gfields={'hygtime','hypack_time';...
    'ptime','mtime';...
    'x','x';...
    'y','y';...
    'hdop','hdop';...
    'numsats','numsats';...
    'gpsmode','mode';...
    'antennaH','height';...
    'gpsTime','time';...
    'lat','latitude';...
    'lon','longitude';...
    'undulation','undulation';...
    'tide','tide'};

%rdata
for i=1:size(gfields,1)
    if isfield(ncdata.gps,gfields{i,2})
        rdata.(gfields{i,1})=ncdata.gps.(gfields{i,2});
    end
end

%bdata
bfields={'mtime','mtime';...
    'vtime','hypack_time';...
    'x','x';...
    'y','y';...
    'distance','distance';...
    'adist','x_dist';...
    'offline','offline';...
    'lat','latitude';...
    'lon','longitude';...
    'zraw','zraw';...
    'tide','tide';...
    'zc','zc'};

bdata.filename=ncdata.gatts.raw_filename;
for i=1:size(bfields,1)
    if isfield(ncdata.bathy,bfields{i,2})
        bdata.(bfields{i,1})=ncdata.bathy.(bfields{i,2});
    end
end

if isfield(ncdata.bathy,'zf')
    cdata.filename=ncdata.gatts.raw_filename;
    cdata.xc=ncdata.bathy.x;
    cdata.yc=ncdata.bathy.y;
    cdata.lat=ncdata.bathy.latitude;
    cdata.lon=ncdata.bathy.longitude;
    cdata.zc=ncdata.bathy.zf;
    cdata.mtime=ncdata.bathy.mtime;
    gsetts=ncdata.bathy.vatts.zf;
end

%hdr- what a pain
hdr.hypack_ver=ncdata.gatts.hypack_version;

if isfield(ncdata.gatts,'ellipsoid')
    hdr.ellipsoid={ncdata.gatts.ellipsoid,...
        sprintf('%0.9f',ncdata.gatts.semimajor_axis),...
        sprintf('%0.9f',ncdata.gatts.inverse_flattening)};
end
if isfield(ncdata.gatts,'projection')
    hdr.projection={ncdata.gatts.projection,...
        sprintf('%0.6f',ncdata.gatts.reference_lon),...
        sprintf('%0.6f',ncdata.gatts.scale_factor),...
        sprintf('%0.6f',ncdata.gatts.reference_lat),...
        sprintf('%0.6f',ncdata.gatts.north_parallel),...
        sprintf('%0.6f',ncdata.gatts.south_parallel),...
        sprintf('%0.4f',ncdata.gatts.false_easting),...
        sprintf('%0.4f',ncdata.gatts.false_northing)};
    hdr.orthometric_height_corr=ncdata.gatts.orthometric_height_corr;
end
if isfield(ncdata.gatts,'geoid')
    hdr.geoid=ncdata.gatts.geoid;
end

hdr.datestr={datestr(bdata.mtime(1),'HH:MM:SS'),...
    datestr(bdata.mtime(1),'mm/dd/yyyy'), '0'};
hdr.mtime=bdata.mtime(1);
hdr.lineNum=ncdata.gatts.line_number;
hdr.linex=ncdata.gatts.line_x;
hdr.liney=ncdata.gatts.line_y;

%deal with the metadata
params={'project';...
    'ops_area';...
    'sci_pi';...
    'cruise_id';...
    'vessel_id';...
    'vessel_operator';...
    'receiver_type';...
    'antenna_type';...
    'gps_antenna_offset';...
    'ppk_filename';...
    'map_projection';...
    'zone';...
    'datum';...
    'geoid';...
    'echosounder_type';...
    'orig_speed_of_sound';...
    'applied_speed_of_sound';...
    'manual_offset';...
    'xducer_frequency';...
    'xducer_blanking';...
    'xducer_gain';...
    'xducer_tx';...
    'base_stn_id';...
    'data_analyst';...
    'note';...
    'xducer_pitch_applied_deg'};
for i=1:length(params)
    if isfield(ncdata.gatts,params{i})
        info.(params{i})=ncdata.gatts.(params{i});
    end
end

%bin data
if isfield(ncdata.sonar,'amplitude')
    bindata.vtime=ncdata.sonar.hypack_time;
    bindata.range=ncdata.sonar.range;
    bindata.vals=ncdata.sonar.amplitude;
end


%speed of sound corrections
if isfield(ncdata.gatts,'orig_speed_of_sound');
    sv.sos_orig=str2double(ncdata.gatts.orig_speed_of_sound);
    sv.depth=[];
end
if isfield(ncdata.gatts,'applied_speed_of_sound');
    if strcmpi(ncdata.gatts.applied_speed_of_sound,'profile')
        sv.depth=ncdata.svel.depth;
        sv.sos=ncdata.svel.svel;
        sv.mean_vel=mean(sv.sos);
        sv.use_prof=1;
        sv.use_mean_sos=0;
    else
        if isfield(ncdata,'svel')
            if isfield(ncdata.svel,'depth')
                sv.depth=ncdata.svel.depth;
                sv.sos=ncdata.svel.svel;
            else
                sv.sos=ncdata.svel.svel;
            end
        else
            sv.depth=[];
            sv.sos=[];
        end
        sv.mean_vel=str2double(ncdata.gatts.applied_speed_of_sound);
        sv.use_prof=0;
        sv.use_mean_sos=1;
    end
    
    
end


%collect the output
gdata.rdata=rdata;
gdata.bdata=bdata;
gdata.hdr=hdr;
if exist('bindata','var')
    gdata.bindata=bindata;
end
if exist('info','var')
    gdata.info=info;
end
if exist('sos_corr','var')
    gdata.sos_corr=sos_corr;
end
if isfield(ncdata.gatts,'manual_offset')
    if ischar(ncdata.gatts.manual_offset)
        gdata.manoff=str2double(ncdata.gatts.manual_offset);
    else
    gdata.manoff=ncdata.gatts.manual_offset;
    end
end
if exist('cdata','var')
    gdata.cdata=cdata;
    gdata.gsetts=gsetts;
end
if exist('sv','var');
    gdata.sv=sv;
end

end
%%%%----------------------------------------------------------------------


function showgps(hfig,evnt) %#ok
gdata=guidata(hfig);

figure
subplot(131)
plot(gdata.rdata.x./1000,...
    gdata.rdata.y./1000,'b.')
hold on
plot(gdata.rdata.x(1)./1000,...
    gdata.rdata.y(1)./1000,'rs',...
    'markerfacecolor','r')
axis equal

xl=get(gca,'xlim');
yl=get(gca,'ylim');
xint=diff(xl)/2;
yint=diff(yl)/3;
xt=xl(1):xint:xl(2);
yt=yl(1):yint:yl(2);
set(gca,'xtick',xt,'ytick',yt)
xlabel('Easting (km)','fontsize',10,...
    'fontweight','bold')
ylabel('Northing (km)','fontsize',10,...
    'fontweight','bold')

minx=min(gdata.rdata.ptime);
maxx=max(gdata.rdata.ptime);
tint=(maxx-minx)/4;

subplot(4,3,2:3)
plot(gdata.rdata.ptime,gdata.rdata.antennaH)
set(gca,'xtick',(minx:tint:maxx),'xticklabel',[],...
    'yaxislocation','right','xlim',[minx maxx])
ylabel('Antenna Height (m)','fontsize',10,...
    'fontweight','bold');

subplot(4,3,5:6)
plot(gdata.rdata.ptime,gdata.rdata.tide)
set(gca,'xtick',(minx:tint:maxx),'xticklabel',[],...
    'yaxislocation','right','xlim',[minx maxx])
ylabel('Tide Correction (m)','fontsize',10,...
    'fontweight','bold');

subplot(4,3,8:9)
plot(gdata.rdata.ptime,gdata.rdata.numsats)
set(gca,'xtick',(minx:tint:maxx),'xticklabel',[],...
    'yaxislocation','right','ylim',...
    [min(gdata.rdata.numsats)-1 max(gdata.rdata.numsats)+1],...
    'xlim',[minx maxx])
ylabel('Satellites','fontsize',10,...
    'fontweight','bold');


subplot(4,3,11:12)
plot(gdata.rdata.ptime,gdata.rdata.gpsmode)
set(gca,'xtick',(minx:tint:maxx),'xticklabel',[],...
    'yaxislocation','right','ylim',[0 5],...
    'ytick',(1:4),'xlim',[minx maxx])
ylabel('GPS mode','fontsize',10,...
    'fontweight','bold');

datetick('x',15,'keepticks','keeplimits')


end
%%%%-----------------------------------------------------------------------
function viewAlong(hfig,evnt) %#ok
gdata=guidata(hfig);
set(gdata.l1,'xdata',gdata.bdata.adist)
set(get(gdata.ax1,'xlabel'),'string','Along-line Distance (m)')
set(get(gdata.ax1,'ylabel'),'string','Elevation (m)')


if ~isempty(gdata.flag);
    if isfield(gdata,'bb')
        set(gdata.bb,'xdata',gdata.bdata.adist(gdata.flag==1));
    end
    if isfield(gdata,'gg')
        
        set(gdata.gg,'xdata',gdata.bdata.adist);
    end
end

if ~isempty(gdata.topo)
    if ishandle(gdata.tlineh)
        set(gdata.tlineh,'visible','on')
    else
        gdata.tlineh=plot(gdata.topo.dist,gdata.topo.z,...
            'g-','linewidth',2);
        hold on
        gdata.l1=plot(gdata.bdata.adist,gdata.bdata.zc);
    end
    gdata.xlimo=[min([gdata.bdata.adist;gdata.topo.dist]),...
        max([gdata.bdata.adist;gdata.topo.dist])];
    gdata.ylimo=[min([gdata.bdata.zc;gdata.topo.z]),...
        max([gdata.bdata.zc;gdata.topo.z])];
    gdata.xlims=[min([gdata.bdata.adist;gdata.topo.dist]),...
        max([gdata.bdata.adist;gdata.topo.dist])];
    gdata.ylims=[min([gdata.bdata.zc;gdata.topo.z]),...
        max([gdata.bdata.zc;gdata.topo.z])];
    
else
    gdata.xlimo=[min(gdata.bdata.adist),...
        max(gdata.bdata.adist)];
    gdata.xlims=[min(gdata.bdata.adist),...
        max(gdata.bdata.adist)];
    gdata.ylims=[min(gdata.bdata.zc),...
        max(gdata.bdata.zc)];
end
set(gca,'xlim',gdata.xlims,'ylim',gdata.ylims)

gdata.alongflag=1;

set(gdata.menu15,'visible','on')
set(gdata.menu14,'visible','off');



guidata(hfig,gdata)
end
%%%%%---------------------------------------------------------------------
function viewAbsolute(hfig,evnt) %#ok
gdata=guidata(hfig);
set(gdata.l1,'xdata',gdata.bdata.distance)
set(get(gdata.ax1,'xlabel'),'string','Distance (m)')
set(get(gdata.ax1,'ylabel'),'string','Elevation (m)')

if ~isempty(gdata.flag);
    if isfield(gdata,'bb')
        set(gdata.bb,'xdata',gdata.bdata.distance(gdata.flag==1));
    end
    set(gdata.gg,'xdata',gdata.bdata.distance);
end

gdata.xlimo=[min(gdata.bdata.distance),...
    max(gdata.bdata.distance)];
gdata.xlims=[min(gdata.bdata.distance),...
    max(gdata.bdata.distance)];
gdata.ylims=[min(gdata.bdata.zc),...
    max(gdata.bdata.zc)];
set(gca,'xlim',gdata.xlims,'ylim',gdata.ylims)

gdata.alongflag=0;

set(gdata.menu15,'visible','off')
set(gdata.menu14,'visible','on');

if ~isempty(gdata.topo)
    set(gdata.tlineh,'visible','off')
end

guidata(hfig,gdata)
end
%%%%%---------------------------------------------------------------------
function ppk_data=run_ppk_gui(hfig,evnt) %#ok

gdata=guidata(hfig);


%load the ppk data output from grafnav
gdata.ppk_info=ppkgui;
drawnow

h = waitbar(0,'Reading PPK GPS file, Please wait...');
set(h,'name','Import PPK GPS data');
if gdata.ppk_info.ftype==2;
    gdata.ppk=readppk(gdata.ppk_info,gdata.hdr);
elseif gdata.ppk_info.ftype==1
    gdata.ppk=readgnavtxt(gdata.ppk_info);
    gdata.ppk_info.base_date=floor(min(gdata.ppk.mtime));
else
end

waitbar(1,h,'Done!');
close(h)

set(gdata.ppkview,'enable','on')
set(gdata.apply_ppk_menu,'enable','on')
set(gdata.plotrtkppk,'visible','on')


guidata(hfig,gdata)

end

%%%%----------------------------------------------------------------------
function ppk = ppkgui


hf = figure('units','normalized',...
    'position',[0.253 0.462 0.251 0.264],...
    'menubar','none','name','ppkgui','numbertitle',...
    'off','color',[0.941 0.941 0.941]);


pushbutton2 = uicontrol(hf,'style','pushbutton',...
    'units','normalized','position',[0.0841 0.785 0.236 0.148],...
    'string','Open','backgroundcolor',[0.941 0.941 0.941],...
    'callback',@getppkfile);
ppkd.file = uicontrol(hf,'style','text','units','normalized',...
    'position',[0.362 0.813 0.602 0.1],'string','Select File to Open',...
    'backgroundcolor',[0.941 0.941 0.941],'horizontalalign','left');


text3=uicontrol(hf,'style','text','units','normalized',...
    'position',[0.094 0.63 0.223 0.124],'string','Base Date',...
    'backgroundcolor',[0.941 0.941 0.941]);
ppkd.date = uicontrol(hf,'style','edit','units','normalized',...
    'position',[0.392 0.65 0.259 0.119],...
    'string',datestr(floor(now)),'backgroundcolor',[1 1 1]);
ppkd.changedate=uicontrol(hf,'style','pushbutton','units','normalized',...
    'position',[0.66 0.65 0.259 0.119],...
    'string','Change Date','backgroundcolor',[1 1 1],...
    'callback',@changedate);

text2 = uicontrol(hf,'style','text','units','normalized',...
    'position',[0.1094 0.47 0.223 0.124],'string','Antenna Height',...
    'backgroundcolor',[0.941 0.941 0.941]);
ppkd.edit1 = uicontrol(hf,'style','edit','units','normalized',...
    'position',[0.392 0.47 0.259 0.119],...
    'string',sprintf('%0.2f',0),'backgroundcolor',[1 1 1]);

ppkd.check1=uicontrol(hf,'style','checkbox','units','normalized',...
    'position',[0.66 0.47 0.259 0.119],...
    'string','Use Ellipsoid Height','backgroundcolor',[1 1 1]);


text1 = uicontrol(hf,'style','text','units','normalized',...
    'position',[0.104 0.266 0.223 0.124],'string','Units',...
    'backgroundcolor',[0.941 0.941 0.941]);
ppkd.pop1 = uicontrol(hf,'style','popupmenu','units','normalized',...
    'position',[0.392 0.261 0.262 0.124],...
    'string',{'Meters';'Survey Feet'},'backgroundcolor',[1 1 1]);



ppkd.done = uicontrol(hf,'style','pushbutton',...
    'units','normalized','position',[0.392 0.081 0.259 0.133],...
    'string','Done','backgroundcolor',[0.941 0.941 0.941],...
    'enable','off','callback',@ppkd_close);


guidata(hf,ppkd)


uiwait
ppkd=guidata(hf);
ppk.filename=ppkd.filename;
ppk.pathname=ppkd.pathname;
ppk.ftype=ppkd.ftype;
ppk.antenna_height=ppkd.antenna_height;
ppk.base_date=ppkd.base_date;
ppk.use_ellipsoid=ppkd.use_ellipsoid;
ppk.units=ppkd.units;

close(hf)
end

%%%%%-------------------------------------------------------------------
function getppkfile(hf,evnt) %#ok

ppkd=guidata(hf);

[filename, pathname, fmt] = uigetfile( ...
    {'*.txt', 'TXT Files (*.txt)';...
    '*.gga', 'GGA Files (*.gga)'},...
    'Select a PPK file');

if filename==0
    return
end

set(ppkd.file,'string',filename)

ppkd.pathname=pathname;
ppkd.filename=filename;
ppkd.ftype=fmt;
if ppkd.ftype==1
    set(ppkd.date,'enable','off')
    set(ppkd.changedate,'enable','off');
else
    set(ppkd.date,'enable','on')
    set(ppkd.changedate,'enable','on');
end


set(ppkd.done,'enable','on')


guidata(hf,ppkd)
end
%%%%----------------------------------------------------------------------
function varargout=readgps(varargin)
%READGPS - reads gps strings from a text file
%
% READGPS(FILENAME) - reads nmea formatted strings contained
%   in the file, FILENAME. Currently only RMC and DPT strings
%   are supported. If no filename is supplied, the user will
%   be prompted interactively to supply a file.
%
%   GPS = READGPS(FILENAME) - returns the gps data in a
%   structure.
%
%   READGPS(...,FILEOUT) - also outputs a text file with
%   parsed output
%
% EXAMPLE 1 - specify both input and output files, request output
%   fin = 'foo.txt';
%   fout = 'gpsout.txt'
%   gps = readgps(fin,fout);
%
% EXAMPLE 2 - specify only output file, no matlab output
%   fout = 'gpsout.txt';
%   readgps([],fout)
%
% Andrew Stevens 05/26/2011
% astevens@usgs.gov

%check inputs
error(nargchk(0,2,nargin,'struct'));
if nargin>0;
    fpath=varargin{1};
    
    if isempty(fpath)
        [filename, pathname] = ...
            uigetfile('*.txt', 'Pick a text file');
        fpath=[pathname,filename];
    else
        if exist(fpath,'file')==0
            error('File not found. Try again.');
        end
    end
    if nargin==2
        base_date=varargin{2};
        
    end
else
    [filename, pathname] = ...
        uigetfile('*.txt', 'Pick a text file');
    fpath=[pathname,filename];
end

%open file and check its length
fid=fopen(fpath);
data = textscan(fid,'%s%[^\n]','delimiter',',');

%define codes and associated formats
codes={'$GPRMC';'$SDDPT';'$GPGGA'};
formats={['%s %*s %f %*s %f %*s',...
    '%f %f %s %*[^\n]'];...
    '%f %*[^\n]';...
    '%f %f %*s %f %*s %f %f %f %f %*s %f %*[^\n]'};


for i=1:length(codes);
    
    data2={data{2}(strcmpi(codes{i},data{1}))}';
    if ~isempty(data2{:})
        r=cellfun(@(x)(textscan(x,formats{i},...
            'delimiter',',')),...
            data2{:},'un',0);
        
        switch codes{i}
            case '$GPRMC'
                %time
                times=cellfun(@(x)(str2double(x{1})),r);
                
                hr=fix(times./10000);
                minute=fix((times-(hr*10000))/100);
                sec=times-(hr*10000+minute*100);
                
                
                dempty=cellfun(@(x)(isempty(x{6})),r);
                r{dempty}{6}='000000';
                
                dates=cellfun(@(x)(textscan(char(x{6}),...
                    '%2.0f%2.0f%2.0f')),r,'un',0);
                dmat=cell2mat(fliplr(cat(1,dates{:})));
                dn=datenum(dmat(:,1)+2000,dmat(:,2),dmat(:,3),...
                    hr,minute,sec);
                dn(dempty)=NaN;
                
                %lat
                lat=cellfun(@(x)(x{2}),r);
                lat1=fix(lat/100);
                lat2= (lat-(lat1*100))/60;
                latd=lat1+lat2;
                
                %lon
                lon=cellfun(@(x)(x{3}),r);
                lon1=fix(lon/100);
                lon2= (lon-(lon1*100))/60;
                lond=-(lon1+lon2);
            case '$SDDPT'
                depth=cell2mat(cat(1,r{:}));
            case '$GPGGA'
                
                
                if ~exist('times','var')
                    %                     times=cellfun(@(x)(textscan(char(x{1}),...
                    %                         '%2.0f%2.0f%2.0f.%2.0f')),r,'un',0);
                    
                    times=cellfun(@(x)(x{1}),r);
                    hr=fix(times./10000);
                    minute=fix((times-(hr*10000))/100);
                    sec=times-(hr*10000+minute*100);
                    
                    
                    
                    if ~exist('base_date','var')
                        warning(['No base date supplied with GGA string.\n',...
                            'Using today''s date.']);
                        base_date=datevec(floor(now));
                    else
                        base_date=datevec(base_date);
                        
                    end
                    
                    dvec=[repmat(base_date(1:3),length(times),1),...
                        hr minute sec];
                    dn=datenum(dvec);
                    
                    %in case dates change
                    idh = find(diff(dn)<-(1/24));
                    if ~isempty(idh)
                        nh = length(idh);
                        for jh = 1:nh
                            dn(idh(jh)+1:end) = dn(idh(jh)+1:end)+1;
                        end
                    end
                    
                    %lat
                    lat=cellfun(@(x)(x{2}),r);
                    lat1=fix(lat/100);
                    lat2= (lat-(lat1*100))/60;
                    latd=lat1+lat2;
                    
                    %lon
                    lon=cellfun(@(x)(x{3}),r);
                    lon1=fix(lon/100);
                    lon2= (lon-(lon1*100))/60;
                    lond=-(lon1+lon2);
                end
                
                %gps quality
                gmode=cellfun(@(x)(x{4}),r); %1-stand alone, 2-diff, 3-PPS fix
                nsats=cellfun(@(x)(x{5}),r); % 4-float 5-RTK fixed
                hdop=cellfun(@(x)(x{6}),r);
                
                %antenna height
                elev=cellfun(@(x)(x{7}),r);
                geoid=cellfun(@(x)(x{8}),r);
                
                
                
        end
    end
end
fclose(fid);

%collect output
gps.mtime=dn;
gps.lon=lond;
gps.lat=latd;

%some vars are optional
if exist('depth','var')
    gps.depth=depth;
end
if exist('gmode','var');
    gps.mode=gmode;
    gps.nsats=nsats;
    gps.hdop=hdop;
    gps.elev=elev;
    gps.geoid=geoid;
end


if nargout==1
    varargout={gps};
end
end

%%%%----------------------------------------------------------------------
function date=changedate(hf,evnt)

ppkd=guidata(hf);
new_date=uical(datenum(get(ppkd.date,'string')));

set(ppkd.date,'string',datestr(new_date));

end

%%%%----------------------------------------------------------------------
function ppk=ppkd_close(hf,evnt) %#ok

ppkd=guidata(hf);
ppkd.ftype=ppkd.ftype;
ppkd.antenna_height=str2double(get(ppkd.edit1,'string'));
ppkd.base_date=datenum(get(ppkd.date,'string'));

ppkd.use_ellipsoid=get(ppkd.check1,'value');
switch get(ppkd.pop1,'value');
    case 1
        ppkd.units='meters';
    case 2
        ppkd.units='Survey feet';
end

guidata(hf,ppkd)
uiresume

end

%%%%%---------------------------------------------------------------------
function view_ppk_data(hfig,evnt) %#ok
gdata=guidata(hfig);

figure
subplot(131)
plot(gdata.ppk.x./1000,...
    gdata.ppk.y./1000,'b.')
hold on
plot(gdata.ppk.x(1)./1000,...
    gdata.ppk.y(1)./1000,'rs',...
    'markerfacecolor','r')
axis equal

xl=get(gca,'xlim');
yl=get(gca,'ylim');
xint=diff(xl)/2;
yint=diff(yl)/3;
xt=xl(1):xint:xl(2);
yt=yl(1):yint:yl(2);
set(gca,'xtick',xt,'ytick',yt)
xlabel('Easting (km)','fontsize',10,...
    'fontweight','bold')
ylabel('Northing (km)','fontsize',10,...
    'fontweight','bold')

minx=min(gdata.ppk.mtime);
maxx=max(gdata.ppk.mtime);
tint=(maxx-minx)/4;

subplot(4,3,2:3)
plot(gdata.ppk.mtime,gdata.ppk.elev)
set(gca,'xtick',(minx:tint:maxx),'xticklabel',[],...
    'yaxislocation','right','xlim',[minx maxx])
ylabel('Orthometric Height (m)','fontsize',10,...
    'fontweight','bold');

subplot(4,3,5:6)
plot(gdata.ppk.mtime,gdata.ppk.geoid)
set(gca,'xtick',(minx:tint:maxx),'xticklabel',[],...
    'yaxislocation','right','xlim',[minx maxx])
ylabel('Geoid Height (m)','fontsize',10,...
    'fontweight','bold');

subplot(4,3,8:9)
plot(gdata.ppk.mtime,gdata.ppk.nsats)
set(gca,'xtick',(minx:tint:maxx),'xticklabel',[],...
    'yaxislocation','right','ylim',...
    [min(gdata.ppk.nsats)-1 max(gdata.ppk.nsats)+1],...
    'xlim',[minx maxx])
ylabel('Satellites','fontsize',10,...
    'fontweight','bold');


subplot(4,3,11:12)
plot(gdata.ppk.mtime,gdata.ppk.mode)
set(gca,'xtick',(minx:tint:maxx),'xticklabel',[],...
    'yaxislocation','right',...
    'ylim',[min(gdata.ppk.mode)-1 max(gdata.ppk.mode)+1],...
    'ytick',(min(gdata.ppk.mode)-1:max(gdata.ppk.mode)+1),...
    'xlim',[minx maxx])
ylabel('GPS mode','fontsize',10,...
    'fontweight','bold');

datetick('x','keepticks','keeplimits')

end

%%%%-----------------------------------------------------------------------
function plot_ppk_rtk(hfig,evnt) %#ok
gdata=guidata(hfig);


if gdata.applyppk~=0;
    if gdata.invert==1
        [gdata.bdata,gdata.hdr,gdata.rdata]=readRAW(...
            'filename',[gdata.filepath,gdata.fnames{gdata.newInd}],...
            'invert','minz',gdata.minz,'maxz',gdata.maxz,...
            'maxoffset',gdata.maxoffset);
    else
        [gdata.bdata,gdata.hdr,gdata.rdata]=readRAW(...
            'filename',[gdata.filepath,gdata.fnames{gdata.newInd}],...
            'minz',gdata.minz,'maxz',gdata.maxz,...
            'maxoffset',gdata.maxoffset);
    end
end




%convert the time in the rdata field (raw file)
if isfield(gdata,'info')
    fixgpstime=1;
    if isfield(gdata.info,'ppk_filename')
        fixgpstime=0;
    end
else
    fixgpstime=1;
end

if fixgpstime
    %convert the time in the rdata field (raw file)
    hr=fix(gdata.rdata.gpsTime./10000);
    minute=fix((gdata.rdata.gpsTime-(hr*10000))/100);
    sec=gdata.rdata.gpsTime-(hr*10000+minute*100);
    
    base_date=datevec(gdata.ppk_info.base_date);
    
    dvec=[repmat(base_date(1:3),length(hr),1),...
        hr minute sec];
    dn=datenum(dvec);
    
    %in case dates change
    idh = find(diff(dn)<-(1/24));
    if ~isempty(idh)
        nh = length(idh);
        for jh = 1:nh
            dn(idh(jh)+1:end) = dn(idh(jh)+1:end)+1;
        end
    end
    
    
else
    dn=gdata.rdata.gpsTime;
end


idx=(gdata.ppk.mtime>=dn(1) & gdata.ppk.mtime<=dn(end));

figure
plot(gdata.ppk.x(idx),...
    gdata.ppk.y(idx),'b.')
hold on
plot(gdata.rdata.x,...
    gdata.rdata.y,'r.')
axis equal

fun=@(x)([min(x) max(x)]);
xl=fun([gdata.rdata.x(:);gdata.ppk.x(idx)]);
yl=fun([gdata.rdata.y(:);gdata.ppk.y(idx)]);

xint=diff(xl)/2;
yint=diff(yl)/3;
xt=xl(1):xint:xl(2);
yt=yl(1):yint:yl(2);
set(gca,'xtick',xt,'ytick',yt,...
    'xlim',xl,'ylim',yl)
xlabel('Easting (m)','fontsize',10,...
    'fontweight','bold')
ylabel('Northing (m)','fontsize',10,...
    'fontweight','bold')

minx=min(dn);
maxx=max(dn);
tint=(maxx-minx)/4;

offset=-mode(gdata.rdata.antennaH+gdata.rdata.tide);

figure
plot(gdata.ppk.mtime(idx),gdata.ppk.elev(idx)-gdata.ppk_info.antenna_height)
hold on
plot(dn,gdata.rdata.antennaH+offset,'r');
% plot(dn,gdata.rdata.antennaH,'r');
leg=legend('PPK','RTK');
set(leg,'location','best')

set(gca,'xtick',(minx:tint:maxx),'xticklabel',[],...
    'yaxislocation','right','xlim',[minx maxx])
ylabel('Orthometric Height (m)','fontsize',10,...
    'fontweight','bold');

datetick('x','keepticks','keeplimits')


end
%%%%%----------------------------------------------------------------------
function set_file_date(hfig,evnt) %#ok

gdata=guidata(hfig);

if isfield(gdata,'ppk_info');
    gdata.raw_base_date=uical(gdata.ppk_info.base_date);
else
    gdata.raw_base_data=uical;
end

guidata(hfig,gdata);

end

%%%%-----------------------------------------------------------------------
function apply_ppk_menu(hfig,evnt) %#ok
gdata=guidata(hfig);

gdata.applyppk=ppkapp_menu(hfig);

guidata(hfig,gdata)

if gdata.applyppk>0
    applyppk(hfig)
end


end
%%%%----------------------------------------------------------------------
function applyppk(hfig,evnt) %#ok


gdata=guidata(hfig);


if isfield(gdata,'raw_base_date')
    base_date=gdata.raw_base_date;
else
    base_date=gdata.ppk_info.base_date;
end

%convert the time in the rdata field (raw file)
if isfield(gdata,'info')
    fixgpstime=1;
    if isfield(gdata.info,'ppk_filename')
        fixgpstime=0;
    end
else
    fixgpstime=1;
end

if fixgpstime
    hr=fix(gdata.rdata.gpsTime./10000);
    minute=fix((gdata.rdata.gpsTime-(hr*10000))/100);
    sec=gdata.rdata.gpsTime-(hr*10000+minute*100);
    
    base_date=datevec(base_date);
    
    dvec=[repmat(base_date(1:3),length(hr),1),...
        hr minute sec];
    dn=datenum(dvec);
    
    %in case dates change
    idh = find(diff(dn)<-(1/24));
    if ~isempty(idh)
        nh = length(idh);
        for jh = 1:nh
            dn(idh(jh)+1:end) = dn(idh(jh)+1:end)+1;
        end
    end
    
    gdata.rdata.gpsTime=dn;
end
%replace position information in raw file with post-process

%interpolate onto rdata.ptime
%first check to make sure ppk covers the time frame
if (gdata.rdata.gpsTime(1)<gdata.ppk.mtime(1) || ...
        gdata.rdata.gpsTime(end)>gdata.ppk.mtime(end))
    warning(['Post-process GPS does not cover entire ',...
        'data record. Tide strings will be truncated.'])
end

fields={'ptime','mtime';...
    'x','x';...
    'y','y';...
    'hdop','hdop';...
    'numsats','nsats';...
    'gpsmode','mode';...
    'antennaH','elev';...
    'gpsTime','mtime';...
    'lat','lat';...
    'lon','lon';...
    'undulation','geoid'};


idx=(gdata.ppk.mtime>=gdata.rdata.gpsTime(1) & ...
    gdata.ppk.mtime<=gdata.rdata.gpsTime(end));
[gpstime_u,ib]=unique(gdata.rdata.gpsTime);

gdata.rdata.hygtime=interp1(gpstime_u,...
    gdata.rdata.hygtime(ib),...
    gdata.ppk.mtime(idx));
for i=1:size(fields,1)
    if isfield(gdata.ppk,fields{i,2})
        gdata.rdata.(fields{i,1})=gdata.ppk.(fields{i,2})(idx);
    end
end

if gdata.ppk_info.use_ellipsoid==1
    if isfield(gdata.ppk,'geoid')
        gdata.rdata.undulation=zeros(size(gdata.rdata.x));
    end
end

%
%
% ptime=gdata.rdata.gpsTime;
% for i=1:size(fields,1)
%     gdata.rdata.(fields{i,1})=interp1(gdata.ppk.mtime,...
%         gdata.ppk.(fields{i,2}),ptime);
% end
gdata.rdata.tide=-(gdata.rdata.antennaH-gdata.ppk_info.antenna_height);



%now mesh ppk data with sounder data
[btime,bind]=unique(gdata.rdata.hygtime);

gdata.bdata.x=interp1(btime,gdata.rdata.x(bind),gdata.bdata.vtime);
gdata.bdata.y=interp1(btime,gdata.rdata.y(bind),gdata.bdata.vtime);

%calculate cumulative distance along line
x1=gdata.bdata.x(1:end-1);
x2=gdata.bdata.x(2:end);
y1=gdata.bdata.y(1:end-1);
y2=gdata.bdata.y(2:end);

xd = x2 - x1;
yd = y2 - y1;
dist=sqrt(xd.*xd +yd.*yd);
dist(isnan(dist))=0;
gdata.bdata.distance=[0,cumsum(dist)'];


% if it is a straight line followed, calculate offline distance
if isfield(gdata.hdr,'linex');
    if length(gdata.hdr.linex)==2
        m=(gdata.hdr.liney(end)-gdata.hdr.liney(1))/...
            (gdata.hdr.linex(end)-gdata.hdr.linex(1));
        b=gdata.hdr.liney(1)-(m*gdata.hdr.linex(1));
        yi=m.*gdata.bdata.x+b;
        
        %determine distance from line to points
        h=gdata.bdata.y-yi;
        ang=(pi/2)-abs(atan(m));
        d=h*sin(ang);
        gdata.bdata.offline=d;
        
        
        %determine along-line distance
        
        theta=atan(m);
        xt = gdata.bdata.x.*cos(theta) + gdata.bdata.y.*sin(theta);
        x0=gdata.hdr.linex.*cos(theta) + ...
            gdata.hdr.liney.*sin(theta);
        
        %calculate along line distance
        gdata.bdata.adist=xt-x0(1);
    else
        gdata.bdata.offline=[];
        gdata.bdata.adist=[];
    end
    
end

if isfield(gdata.rdata,'lon');
    gdata.bdata.lon=interp1(btime,gdata.rdata.lon(bind),gdata.bdata.vtime);
    gdata.bdata.lat=interp1(btime,gdata.rdata.lat(bind),gdata.bdata.vtime);
end


gdata.bdata.tide=interp1(btime,gdata.rdata.tide(bind),gdata.bdata.vtime);

if isfield(gdata,'sv')
    if gdata.invert==1
        zraw=-gdata.bdata.zraw;
    end
    zsos=apply_sos_prof(zraw,gdata.sv.sos_orig,...
        [gdata.sv.depth gdata.sv.sos],gdata.sv.use_mean_sos,...
        gdata.sv.mean_vel)-gdata.manoff;
    if gdata.invert==1
        zsos=-zsos;
    end
    
    gdata.bdata.zc=zsos-gdata.bdata.tide;
    gdata.info.orig_speed_of_sound=sprintf('%0.1f',...
        gdata.sv.sos_orig);
    if gdata.sv.use_mean_sos
        gdata.info.applied_speed_of_sound=sprintf('%0.1f',...
            gdata.sv.mean_vel);
    else
        gdata.info.applied_speed_of_sound='Profile';
    end
    
    
    
else
    gdata.bdata.zc=(gdata.bdata.zraw+gdata.manoff)-...
        gdata.bdata.tide;
end

gdata.xlimo=[nanmin(gdata.bdata.distance),...
    nanmax(gdata.bdata.distance)];
gdata.ylimo=[nanmin(gdata.bdata.zc),...
    nanmax(gdata.bdata.zc)];

%update the view
hold off
if gdata.alongflag==1
    gdata.alongflag=0;
end



gdata.l1=plot(gdata.bdata.distance,gdata.bdata.zc);
hold on

set(gca,'xlim',gdata.xlimo,'ylim',gdata.ylimo);


if isfield(gdata,'gg')
    gdata=rmfield(gdata,'gg');
end


gdata.info.ppk_filename=gdata.ppk_info.filename;
guidata(hfig,gdata);

end



%%%%-----------------------------------------------------------------------
function applyppk = ppkapp_menu(hfig,evnt) %#ok

gdata=guidata(hfig);
hf = figure('units','normalized','position',[0.253 0.456 0.141 0.168],...
    'menubar','none','name','ppkapp_menu',...
    'numbertitle','off','color',[0.941 0.941 0.941]);

appk.bg1 = uibuttongroup('parent',hf,'units','normalized',...
    'position',[0.107 0.247 0.754 0.6],'title','Apply PPK Data');
radiobutton1 = uicontrol(appk.bg1,'style','radiobutton',...
    'units','normalized','position',[0.0748 0.107 0.707 0.205],...
    'string','Off','backgroundcolor',[0.941 0.941 0.941]);
radiobutton2 = uicontrol(appk.bg1,'style','radiobutton',...
    'units','normalized','position',[0.0701 0.366 0.707 0.205],...
    'string','This File','backgroundcolor',[0.941 0.941 0.941]);
radiobutton3 = uicontrol(appk.bg1,'style','radiobutton',...
    'units','normalized','position',[0.0748 0.634 0.707 0.205],...
    'string','All Files in LOG','backgroundcolor',[0.941 0.941 0.941]);

if isfield(gdata,'applyppk')
    switch gdata.applyppk
        case 0
            set(appk.bg1,'selectedobject',radiobutton1)
        case 1
            set(appk.bg1,'selectedobject',radiobutton2)
        case 2
            set(appk.bg1,'selectedobject',radiobutton3)
    end
end


uicontrol(hf,'style','pushbutton','units',...
    'normalized','position',[0.343 0.0651 0.322 0.144],...
    'string','Done','backgroundcolor',[0.941 0.941 0.941],...
    'callback',@close_applyppk);



guidata(hf,appk)

uiwait
appk=guidata(hf);
applyppk=appk.type;

close(hf)
end

%%%%%--------------------------------------------------------------------
function close_applyppk(hf,evnt) %#ok

appk=guidata(hf);
str=get(get(appk.bg1,'selectedObject'),'string');
switch lower(str)
    case 'off'
        appk.type=0;
    case 'this file'
        appk.type=1;
    case 'all files in log'
        appk.type=2;
end
guidata(hf,appk);
uiresume
end
%%%%%---------------------------------------------------------------------
function ppk=readppk(ppk_info,hdr)
%READPPK - reads post-processed gps data in GGA format
%          Uses info in the .RAW file header to define the
%          projected coordinate system

fpath=[ppk_info.pathname,ppk_info.filename];
if exist(fpath,'file')==0
    errordlg('File not found.')
    return
else
    ppk=readgps(fpath,ppk_info.base_date);
    
end


%construct the transformation
figure
ptype=hdr.projection{1};
pdata=cellfun(@(x)(str2double(x)),hdr.projection(2:end));

%projection (only deals with transverse mercator or
% lamebert conformical
switch ptype
    case 'LCC'
        mstruct=gcm(axesm('lambertstd'));
        
    case 'TME'
        mstruct=gcm(axesm('tranmerc'));
end
close
mstruct.origin=[pdata(3) pdata(1)];
mstruct.mapparallels=[pdata(4) pdata(5)];
mstruct.falseeasting=pdata(6);
mstruct.falsenorthing=pdata(7);
mstruct.scalefactor=pdata(2);


%datum, only deal with wgs84 and grs80
datum=hdr.ellipsoid{1};
switch datum
    case 'WGS-84'
        ref=referenceEllipsoid('WGS84',ppk_info.units);
    case 'GRS-80'
        ref=referenceEllipsoid('GRS80',ppk_info.units);
end
mstruct.geoid=[ref.SemimajorAxis ref.Eccentricity];

%project the coordinates in the gga string
[ppk.x,ppk.y]=projfwd(mstruct,ppk.lat,ppk.lon);

if ppk_info.use_ellipsoid==1
    ppk.elev=ppk.elev+ppk.geoid;
end


end

%%%%----------------------------------------------------------------------
function ppk=readgnavtxt(ppk_info)
%reads custom gravnav text file output

fpath=[ppk_info.pathname,ppk_info.filename];
if exist(fpath,'file')==0
    errordlg('File not found.')
    return
else
    fid=fopen(fpath,'r');
    
end

hdrs=textscan(fgetl(fid),'%s',...
    'delimiter',',');

fmt=['%s%s',repmat('%f',1,10)];
data = textscan(fid,fmt,'delimiter',',');
obs=cell2struct(data',hdrs{1});

obs.mtime=datenum(strcat(cellfun(@(x)([x,' ']),obs.Date,'un',0),obs.Time),...
    'mm/dd/yyyy HH:MM:SS.FFF');

%transform the relavant parameters
ofields={'mtime';...
    'Latitude';...
    'Longitude';...
    'Quality';...
    'Num_Sat';...
    'PDOP';...
    'Ortho_height';...
    'Easting';...
    'Northing'};
nfields={'mtime',;...
    'lat';...
    'lon';...
    'mode';...
    'nsats';...
    'hdop';...
    'elev';...
    'x';...
    'y'};
for i=1:length(ofields)
    ppk.(nfields{i})=obs.(ofields{i});
end

ppk.geoid=obs.Ellip_Height-obs.Ortho_height;


if ppk_info.use_ellipsoid==1
    ppk.elev=obs.Ellip_Height;
end


end

%%%-----------------------------------------------------------------------
function ini=read_info(varargin)
% READ_INFO - read info parameter file
%
% file should be a list of property/values pairs
% separated by an equal sign, tab delimited

if nargin>0
    fname=varargin{1};
    if ~exist(fname,'file')
        error('File not found.')
    end
else
    [filename, pathname] = uigetfile( ...
        {'*.*', 'All Files (*.*)'},...
        'Select a File');
    if filename==0
        ini=[];
        return
    else
        fname=[pathname,filename];
    end
end



fid=fopen(fname);
names = textscan(fid,'%s%*[^\n]');
names=names{1};
frewind(fid);
data=textscan(fid,'%*s =%[^\n]');
data=data{1};

params={'project','%s';...
    'ops_area','%s';...
    'sci_pi','%s';...
    'cruise_id','%s';...
    'vessel_id','%s';...
    'vessel_operator','%s';...
    'receiver_type','%s';...
    'antenna_type','%s';...
    'echosounder_type','%s';...
    'datum','%s';...
    'map_projection','%s';...
    'zone','%s';...
    'geoid','%s';...
    'orig_speed_of_sound','%f';...
    'xducer_frequency','%f';...
    'xducer_blanking','%f';...
    'xducer_gain','%f';...
    'xducer_tx','%f';...
    'base_stn_id','%s';...
    'data_analyst','%s';...
    'note','%s';
    'xducer_pitch_applied_deg','%s'};

for i=1:length(names);
    ind=find(strcmpi(names{i},params(:,1)), 1);
    
    if ~isempty(ind)
        if ~isempty(data{i})
            data2=textscan(data{i},params{ind,2},'delimiter','\t');
            switch params{ind,2}
                case '%s'
                    ini.(names{i})=data2{:}{1};
                case '%f'
                    ini.(names{i})=data2{:};
            end
            
        else
            ini.(names{i})=[];
        end
    else
        errordlg(sprintf('Unknown keyword encountered: %s',...
            names{i}),'modal');
        error('ASTEVENS:read_ini:badString',...
            'Unknown keyword encountered: %s',...
            names{i})
        
    end
end

fclose(fid);
end
%%%%----------------------------------------------------------------------
function run_read_info(hfig,evnt) %#ok
gdata=guidata(hfig);

info=read_info;
if isfield(gdata,'xducer_settings')
    
    info.xducer_frequency=gdata.xducer_settings.Channel_1_frequency./10;
    info.xducer_blanking= gdata.xducer_settings.Blanking./10;
    info.xducer_gain= gdata.xducer_settings.Channel_1_Gain;
    info.xducer_tx=gdata.xducer_settings.Channel_1_TxPower;
end

if isfield(gdata,'info')
    if isfield(gdata.info,'ppk_filename')
        info.ppk_filename=gdata.info.ppk_filename;
    end
    if isfield(gdata.info,'gps_antenna_offset')
        info.gps_antenna_offset=gdata.info.gps_antenna_offset;
    end
end

if isfield(gdata,'sv')
    if ~isempty(gdata.sv)
        if gdata.sv.use_mean_sos
            info.applied_speed_of_sound=sprintf('%0.1f',...
                gdata.sv.mean_vel);
        else
            info.applied_speed_of_sound='Profile';
        end
    end
end

info.tv_ver=gdata.tv_ver;
gdata.info=info;
guidata(hfig,gdata)
show_meta_info(hfig);

set(gdata.showmeta,'visible','on')
end

%%%%%---------------------------------------------------------------------
function show_meta_info(hfig,evnt)%#ok

gdata=guidata(hfig);
h2=figure;
set(h2,'Name','Metadata','menubar','none','numbertitle','off');

if ~isempty(gdata.info)
    fields=fieldnames(gdata.info);
    cdata=struct2cell(gdata.info);
    
    data(1,1:length(fields))=fields;
    for i=1:length(fields)
        if isnumeric(cdata{i})
            cdata{i}=sprintf('%0.3f',cdata{i});
        end
    end
    data(2,1:length(fields))=cdata;
    
    
    mtable=uitable(h2,'Data',data','columnname',{'Property','Value'});
    set(mtable,'units','normalized','position',[0.1 0.1 0.8 0.8]);
    
    pix=get(h2,'position');
    cwidth=floor([0.3 0.45]*pix(3));
    set(mtable,'columnwidth',num2cell(cwidth))
end

end
%%%%%---------------------------------------------------------------------
function load_xducer(hfig,evnt) %#ok
gdata=guidata(hfig);

[filename, pathname] = uigetfile( ...
    {'*.xml', 'XML Files (*.xml)'},...
    'Select a .xml file',gdata.filepath);
if filename==0
    return
end
fname=fullfile(pathname,filename);
gdata.xducer_settings=read_echart_settings(fname);

set(gdata.xdmenu,'visible','on');
gdata.info.orig_speed_of_sound=sprintf('%0.1f',gdata.xducer_settings.Velocity);

gdata.info.xducer_frequency=gdata.xducer_settings.Channel_1_frequency./10;
gdata.info.xducer_blanking= gdata.xducer_settings.Blanking./10;
gdata.info.xducer_gain= gdata.xducer_settings.Channel_1_Gain;
gdata.info.xducer_tx=gdata.xducer_settings.Channel_1_TxPower;


guidata(hfig,gdata);

end

%%%-----------------------------------------------------------------------
function showXducerInfo(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);
h2=figure;
set(h2,'Name','Transducer Settings','menubar','none','numbertitle','off');


fields=fieldnames(gdata.xducer_settings);
cdata=struct2cell(gdata.xducer_settings);

data(1,1:length(fields))=fields;
data(2,1:length(fields))=cellfun(@(x)(sprintf('%0.0f',x)),cdata,'un',0);


mtable=uitable(h2,'Data',data','columnname',{'Property','Value'});
set(mtable,'units','normalized','position',[0.1 0.1 0.8 0.8]);

pix=get(h2,'position');
cwidth=floor([0.3 0.45]*pix(3));
set(mtable,'columnwidth',num2cell(cwidth))


end

%%%%----------------------------------------------------------------------
function showRawData(hfig,evnt) %#ok
gdata=guidata(hfig);

gdata.rawflag=1;
set(gdata.hp,'Title','Raw Data Controls');
set(gdata.list1,'visible','off');
set(gdata.text2,'visible','off');
set(gdata.edit2,'visible','off');
set(gdata.text3,'visible','off');
set(gdata.edit3,'visible','off');
set(gdata.text5,'visible','off');
set(gdata.edit4,'visible','off');
set(gdata.text4,'visible','off');
set(gdata.push5,'callback',@cleanProfileRaw);
set(gdata.menu11,'callback',@cleanProfileRaw);
set(gdata.menu8,'visible','on')
set(gdata.menu9,'visible','off')
set(gdata.rawe1,'visible','on')
set(gdata.rawt1,'visible','on')
set(gdata.rawe2,'visible','on')
set(gdata.rawt2,'visible','on')
set(gdata.rawe3,'visible','on')
set(gdata.rawt3,'visible','on')
set(gdata.rawp1,'visible','on')
set(gdata.menu15,'visible','off')
set(gdata.menu14,'visible','off')
set(gdata.menu16,'visible','on')


if isempty(gdata.bindata);
    try
        
        rawfile=char(strcat(gdata.filepath,strtok(gdata.fnames(...
            gdata.newInd),'.'),'.bin'));
        
        set(gdata.text1,'string','Loading .BIN file...')
        drawnow
        gdata.bindata=readBin(rawfile,gdata.bdata);
        
        
        
        gdata.bindata.range=-gdata.bindata.range.*gdata.invert;
        txtstr= ['Displaying transect ',...
            num2str(gdata.pointer),' of ',...
            num2str(length(gdata.fnames)),' - '...
            gdata.fnames{gdata.newInd}];
        
        set(gdata.text1,'string',txtstr,'foregroundcolor','k');
        
        %make sure bindata and digitized data have same dimensions
        [junk,ia,ib] = intersect(gdata.bdata.vtime,...
            gdata.bindata.vtime);%#ok
        
        gdata.bdata.mtime=gdata.bdata.mtime(ia);
        gdata.bdata.vtime=gdata.bdata.vtime(ia);
        gdata.bdata.x=gdata.bdata.x(ia);
        gdata.bdata.y=gdata.bdata.y(ia);
        gdata.bdata.distance=gdata.bdata.distance(ia);
        if ~isempty(gdata.bdata.adist)
            gdata.bdata.adist=gdata.bdata.adist(ia);
            gdata.bdata.offline=gdata.bdata.offline(ia);
        end
        if isfield(gdata.bdata,'lon')
            gdata.bdata.lon=gdata.bdata.lon(ia);
            gdata.bdata.lat=gdata.bdata.lat(ia);
        end
        gdata.bdata.zraw=gdata.bdata.zraw(ia);
        gdata.bdata.tide=gdata.bdata.tide(ia);
        gdata.bdata.zc=gdata.bdata.zc(ia);
        
        if gdata.numedits~=0
            gdata.eraw=gdata.eraw(ia);
            gdata.ezc=gdata.ezc(ia);
            
            minz=find(ia==1,1,'first');
            if ~isempty(minz)
                for i=1:length(gdata.numedits)
                    gdata.edits{i}=gdata.edits{i}-(minz-1);
                end
            end
        end
        gdata.bindata.vtime=gdata.bindata.vtime(ib);
        gdata.bindata.vals=gdata.bindata.vals(:,ib);
        
        
        
    catch %#ok
        cla
        txtstr= ['Problems loading raw data file:',...
            gdata.fnames{gdata.newInd}];
        
        set(gdata.text1,'string',txtstr,'foregroundcolor','r');
        gdata.rawflag=0;
    end
end


%limit bin data to max val defined in gui
ind=find(gdata.bindata.range>=gdata.maxrawz);
gdata.bindata.range=gdata.bindata.range(ind);
gdata.bindata.vals=gdata.bindata.vals(ind,:);


set(gcf,'renderer','zbuffer')
hold off
gdata.rim=imagesc((1:numel(gdata.bindata.vtime)),...
    gdata.bindata.range,...
    gdata.bindata.vals);

if isfinite(gdata.rclims)
    caxis(gdata.rclims);
end


colormap(flipud(gray));
c1=colorbar;
set(get(c1,'ylabel'),'string','Backscatter Intensity',...
    'fontsize',12)

xlabel('Ping Number','fontsize',12)
ylabel('Depth (m)','fontsize',12);
set(gca,'ydir','normal')
hold on

gdata.l1=plot(1:numel(gdata.bdata.zraw),gdata.bdata.zraw,'r-');
gdata.xlimr=get(gca,'xlim');
gdata.ylimr=get(gca,'ylim');


    if isfield(gdata,'gg')
        
       gdata=rmfield(gdata,'gg');
    end

guidata(hfig,gdata)


end

%%%%-----------------------------------------------------------------------
function showCleanData(hfig,evnt) %#ok
gdata=guidata(hfig);


set(gdata.hp,'Title','Filter Controls');
set(gdata.list1,'visible','on');
set(gdata.text2,'visible','on');
set(gdata.edit2,'visible','on');
set(gdata.text3,'visible','on');
set(gdata.edit3,'visible','on');
set(gdata.text5,'visible','on');
set(gdata.edit4,'visible','on');
set(gdata.text4,'visible','on');
set(gdata.push5,'callback',@cleanProfile);
set(gdata.menu8,'visible','off')
set(gdata.menu9,'visible','on')
set(gdata.menu11,'callback',@cleanProfile);
set(gdata.rawe1,'visible','off')
set(gdata.rawt1,'visible','off')
set(gdata.rawe2,'visible','off')
set(gdata.rawt2,'visible','off')
set(gdata.rawe3,'visible','off')
set(gdata.rawt3,'visible','off')
set(gdata.rawp1,'visible','off')
set(gdata.menu16,'visible','off')

if gdata.alongflag==1
    
    axes(gdata.ax1);
    hold off
    gdata.l1=plot(gdata.bdata.adist,gdata.bdata.zc);
    set(gdata.menu15,'visible','on')
else
    axes(gdata.ax1);
    hold off
    gdata.l1=plot(gdata.bdata.distance,gdata.bdata.zc);
    if ~isempty(gdata.bdata.adist)
        set(gdata.menu14,'visible','on')
    end
end



if ~isempty(gdata.topo) && gdata.alongflag==1
    hold on
    gdata.tlineh=plot(gdata.topo.dist,gdata.topo.z,'r-',...
        'linewidth',2);
    gdata.xlimo=[min([gdata.bdata.adist;gdata.topo.dist]),...
        max([gdata.bdata.adist;gdata.topo.dist])];
    gdata.ylimo=[min([gdata.bdata.zc;gdata.topo.z]),...
        max([gdata.bdata.zc;gdata.topo.z])];
    
    
end

set(gca,'xlim',gdata.xlimo,'ylim',gdata.ylimo)
xlabel('Distance (m)')
ylabel('Elevation (m)')

gdata.rawflag=0;
guidata(hfig,gdata);
end

%%%%-----------------------------------------------------------------------
function cpfcn(hfig,evnt)
gdata=guidata(hfig);
pan off
switch evnt.Character
    case 'e'
        if gdata.rawflag==1
            cleanProfileRaw(hfig);
        else
            cleanProfile(hfig);
        end
    case 'd'
        if gdata.rawflag==1
            handdigitize(hfig);
        end
    case 'v'
        if gdata.rawflag==1
            showCleanData(hfig)
        else
            if isempty(gdata.bindata)
                rawfile=char(strcat(gdata.filepath,strtok(gdata.fnames(...
                    gdata.newInd),'.'),'.bin'));
                if exist(rawfile,'file');
                    showRawData(hfig)
                end
            else
                showRawData(hfig)
            end
        end
    case 'g'
        showgps(hfig)
    case 'z'
        zoomto(hfig)
    case 'f'
        fullExtents(hfig)
    case 'o'
        showOffline(hfig)
    case 'w'
        calcHs(hfig)
    case 'u'
        if gdata.numedits>0
            undo(hfig)
        end
    case 'l'
        applylf(hfig)
    case '-'
        zoomout(hfig)
    case '='
        zoomin(hfig)
    case 'm'
        measure_dist(hfig)
    case ''
        if strcmp(evnt.Key,'control')
            disp('yes')
        end
    case 'b'
        exportAll(hfig)
    case 'a'
        if gdata.alongflag
            viewAbsolute(hfig)
        else
            viewAlong(hfig)
        end
end
end
%%%%%----------------------------------------------------------------------
function measure_dist(hfig,evnt) %#ok

gdata=guidata(hfig);

gdata.told=get(gdata.text1,'string');
set(gdata.text1,'string',...
    'Click and drag to measure distance',...
    'visible','on');


waitforbuttonpress;



cc = get(gdata.ax,'CurrentPoint');
gdata.xy1 = cc(1,1:2);

hold on
guidata(hfig,gdata)

set(hfig,'WindowButtonMotionFcn',@lmotion);
set(hfig,'WindowButtonUpFcn',@ldone);


uiwait
end
%%%%--------------------------------------------------------------

function lmotion(hfig,evnt) %#ok
%LMOTION- callback for buttondown
gdata=guidata(hfig);



mousenew = get(gdata.ax,'CurrentPoint');
gdata.xy2 = mousenew(1,1:2);

xd = gdata.xy2(1) - gdata.xy1(1);
yd = gdata.xy2(2) - gdata.xy1(2);
dist=sqrt(xd.*xd +yd.*yd)';



if isfield(gdata,'m1')
    set(gdata.m1,'xdata',[gdata.xy1(1) gdata.xy2(1)],...
        'ydata',[gdata.xy1(2) gdata.xy2(2)]);
    
    
else
    gdata.m1=line([gdata.xy1(1) gdata.xy2(1)],...
        [gdata.xy1(2) gdata.xy2(2)]);
    set(gdata.m1,'color','r')
    
    
end

set(gdata.text1,'string',sprintf(['dx = %.2f, ',...
    'dy = %0.2f ,distance = %0.2f'],...
    xd,yd,dist));

guidata(hfig,gdata)
end
%----------------------------------------------------------

function ldone(hfig,evnt) %#ok
%LDONE- callback for buttonup
gdata=guidata(hfig);

delete(gdata.m1);

set(hfig,'WindowButtonMotionFcn',[]);
set(hfig,'WindowButtonUpFcn',[]);

set(gdata.text1,'string',gdata.told)
set(hfig,'pointer','arrow')
gdata=rmfield(gdata,{'m1','told'});
guidata(hfig,gdata)
uiresume
end
%%%%%--------------------------------------------------------------------
function scrollfcn(hfig,evnt) %#ok

gdata=guidata(hfig);



mousenew = get(gdata.ax1,'CurrentPoint');
xy=mousenew(1,1:2);
xl=get(gdata.ax1,'xlim');
yl=get(gdata.ax1,'ylim');

xrange=diff(xl);
yrange=diff(yl);
gdata.xlims=[xy(1)-(xrange/2) xy(1)+(xrange/2)];
gdata.ylims=[xy(2)-(yrange/2) xy(2)+(yrange/2)];



if gdata.rawflag
    xdata=1:numel(gdata.bindata.vtime);
    ydata=gdata.bindata.range;
    
    xind=floor(gdata.xlims(1)):ceil(gdata.xlims(2));
    yind=find(ydata>=gdata.ylims(1) & ydata<=gdata.ylims(2));
    
    ia=intersect(xind,xdata);
    
    if all(~isempty([ia(:);yind(:)]));
        assignin('base','bz',gdata.binzoom)
        gdata.binzoom.x=xdata(ia);
        gdata.binzoom.y=ydata(yind);
        gdata.binzoom.cdata=gdata.bindata.vals(yind,ia);
        set(gdata.rim,'xdata',gdata.binzoom.x,...
            'ydata',gdata.binzoom.y,...
            'cdata',gdata.binzoom.cdata);
        guidata(hfig,gdata)
    end
    
    
end


set(gdata.ax1,'xlim',gdata.xlims,...
    'ylim',gdata.ylims)



guidata(hfig,gdata);




end
%%%%%----------------------------------------------------------------------
function zoomout(hfig,evnt) %#ok

gdata=guidata(hfig);
mousenew = get(gdata.ax1,'CurrentPoint');
xy=mousenew(1,1:2);

xl=get(gdata.ax1,'xlim');
xrange=diff(xl)*1.1;
xd=get(gdata.l1,'xdata');
yd=get(gdata.l1,'ydata');

gdata.xlims=[xy(1)-(xrange/2) xy(1)+(xrange/2)];

if gdata.rawflag~=1
    if (xy(1)-(xrange/2))<gdata.xlimo(1);
        gdata.xlims=[gdata.xlimo(1) xy(1)+(xrange/2)+...
            (gdata.xlimo(1)-(xy(1)-(xrange/2)))];
    end
    if (xy(1)+(xrange/2))>gdata.xlimo(2);
        gdata.xlims=[xy(1)-(xrange/2)-(xy(1)+(xrange/2))>gdata.xlimo(2),...
            gdata.xlimo(2)];
    end
else
    if (xy(1)-(xrange/2))<gdata.xlimr(1);
        gdata.xlims=[gdata.xlimr(1) xy(1)+(xrange/2)+...
            (gdata.xlimr(1)-(xl(1)-(xrange/2)))];
    end
    if (xy(1)+(xrange/2))>gdata.xlimr(2);
        gdata.xlims=[xy(1)-(xrange/2)-(xy(1)+(xrange/2))>gdata.xlimr(2),...
            gdata.xlimr(2)];
    end
end


yind=find(xd>=gdata.xlims(1) & xd<=gdata.xlims(2));
gdata.ylims=[min(yd(yind)) max(yd(yind))];

if gdata.rawflag
    xdata=1:numel(gdata.bindata.vtime);
    ydata=gdata.bindata.range;
    
    xind=floor(gdata.xlims(1)):ceil(gdata.xlims(2));
    yind=find(ydata>=gdata.ylims(1) & ydata<=gdata.ylims(2));
    
    ia=intersect(xind,xdata);
    
    if all(~isempty([ia';yind]));
        assignin('base','bz',gdata.binzoom)
        gdata.binzoom.x=xdata(ia);
        gdata.binzoom.y=ydata(yind);
        gdata.binzoom.cdata=gdata.bindata.vals(yind,ia);
        set(gdata.rim,'xdata',gdata.binzoom.x,...
            'ydata',gdata.binzoom.y,...
            'cdata',gdata.binzoom.cdata);
        guidata(hfig,gdata)
    end
    
    
end

set(gdata.ax1,'xlim',gdata.xlims,...
    'ylim',gdata.ylims)

guidata(hfig,gdata);



end
%%%%%----------------------------------------------------------------------
function zoomin(hfig,evnt) %#ok

gdata=guidata(hfig);
mousenew = get(gdata.ax1,'CurrentPoint');
xy=mousenew(1,1:2);

xl=get(gdata.ax1,'xlim');
xrange=diff(xl)*0.9;
xd=get(gdata.l1,'xdata');
yd=get(gdata.l1,'ydata');

gdata.xlims=[xy(1)-(xrange/2) xy(1)+(xrange/2)];

if gdata.rawflag~=1
    if (xy(1)-(xrange/2))<gdata.xlimo(1);
        gdata.xlims=[gdata.xlimo(1) xy(1)+(xrange/2)+...
            (gdata.xlimo(1)-(xy(1)-(xrange/2)))];
    end
    if (xy(1)+(xrange/2))>gdata.xlimo(2);
        gdata.xlims=[xy(1)-(xrange/2)-(xy(1)+(xrange/2))>gdata.xlimo(2),...
            gdata.xlimo(2)];
    end
else
    if (xy(1)-(xrange/2))<gdata.xlimr(1);
        gdata.xlims=[gdata.xlimr(1) xy(1)+(xrange/2)+...
            (gdata.xlimr(1)-(xl(1)-(xrange/2)))];
    end
    if (xy(1)+(xrange/2))>gdata.xlimr(2);
        gdata.xlims=[xy(1)-(xrange/2)-(xy(1)+(xrange/2))>gdata.xlimr(2),...
            gdata.xlimr(2)];
    end
end

yind=find(xd>=gdata.xlims(1) & xd<=gdata.xlims(2));
gdata.ylims=[min(yd(yind)) max(yd(yind))];

if gdata.rawflag
    xdata=1:numel(gdata.bindata.vtime);
    ydata=gdata.bindata.range;
    
    xind=floor(gdata.xlims(1)):ceil(gdata.xlims(2));
    yind=find(ydata>=gdata.ylims(1) & ydata<=gdata.ylims(2));
    
    ia=intersect(xind,xdata);
    
    if all(~isempty([ia';yind]));
        assignin('base','bz',gdata.binzoom)
        gdata.binzoom.x=xdata(ia);
        gdata.binzoom.y=ydata(yind);
        gdata.binzoom.cdata=gdata.bindata.vals(yind,ia);
        set(gdata.rim,'xdata',gdata.binzoom.x,...
            'ydata',gdata.binzoom.y,...
            'cdata',gdata.binzoom.cdata);
        guidata(hfig,gdata)
    end
    
    
end

set(gdata.ax1,'xlim',gdata.xlims,...
    'ylim',gdata.ylims)

guidata(hfig,gdata);



end

%%%%%----------------------------------------------------------------------

function movePhoto(hfig,evnt)%#ok

gdata=guidata(hfig);
back=get(gdata.push1,'value');
forw=get(gdata.push2,'value');
set(gdata.text1,'string','Loading file...');
pause(0.000001)

pNum=str2double(get(gdata.edit1,'string'));


%Do nothing if invalid number is supplied to edit box 1
if isnan(pNum)
    pNum=gdata.pointer;
end


% move pointer forward or back accoring to pushbuttons
if forw==1
    set(gdata.push1,'enable','on')
    if gdata.fInd(gdata.pointer)~=gdata.fInd(end)
        gdata.pointer=gdata.pointer+1;
        gdata.newInd=gdata.fInd(gdata.pointer);
        
        if gdata.pointer==gdata.fInd(end);
            set(gdata.push2,'enable','off');
        end
        
    else
        txtstr= ['Displaying transect ',...
            num2str(gdata.pointer),' of ',...
            num2str(length(gdata.fnames)),' - '...
            gdata.fnames{gdata.newInd}];
        set(gdata.text1,'string',txtstr,'foregroundcolor','k');
        return
    end
elseif back==1;
    set(gdata.push2,'enable','on');
    if gdata.fInd(gdata.pointer)~=gdata.fInd(1)
        gdata.pointer=gdata.pointer-1;
        gdata.newInd=gdata.fInd(gdata.pointer);
        
        if gdata.pointer==1;
            set(gdata.push1,'enable','off')
        end
    else
        txtstr= ['Displaying transect ',...
            num2str(gdata.pointer),' of ',...
            num2str(length(gdata.fnames)),' - '...
            gdata.fnames{gdata.newInd}];
        set(gdata.text1,'string',txtstr,'foregroundcolor','k');
        return
    end
else
    if pNum>=length(gdata.fInd);
        pNum=length(gdata.fInd);
        set(gdata.edit1,'string',num2str(pNum));
        set(gdata.push2,'enable','off')
        set(gdata.push1,'enable','on')
    elseif pNum<=1
        pNum=1;
        set(gdata.edit1,'string',num2str(1));
        set(gdata.push1,'enable','off');
        set(gdata.push2,'enable','on');
    end
    gdata.pointer=pNum;
    gdata.newInd=gdata.fInd(gdata.pointer);
    
end




%update view

set(gdata.edit1,'string',num2str(gdata.pointer));

%
% try

if gdata.invert==1
    [gdata.bdata,gdata.hdr,gdata.rdata]=readRAW(...
        'filename',[gdata.filepath,gdata.fnames{gdata.newInd}],...
        'invert','minz',gdata.minz,'maxz',gdata.maxz,...
        'maxoffset',gdata.maxoffset);
else
    [gdata.bdata,gdata.hdr,gdata.rdata]=readRAW(...
        'filename',[gdata.filepath,gdata.fnames{gdata.newInd}],...
        'minz',gdata.minz,'maxz',gdata.maxz,...
        'maxoffset',gdata.maxoffset);
end

%deal with tide options
gdata.tideopt.badtide=[];
guidata(hfig,gdata);
applytidecorr(hfig);
gdata=guidata(hfig);

hold off
if isfield(gdata.bdata,'zc')==0
    gdata.bdata.zc=gdata.bdata.zraw;
    gdata.bdata.tide=zeros(size(gdata.bdata.zraw));
    gdata.rdata.tide=zeros(size(gdata.rdata.antennaH));
end

if isfield(gdata,'sv');
    if gdata.invert==1
        zraw=-gdata.bdata.zraw;
    end
    zsos=apply_sos_prof(zraw,gdata.sv.sos_orig,...
        [gdata.sv.depth gdata.sv.sos],gdata.sv.use_mean_sos,...
        gdata.sv.mean_vel)-gdata.manoff;
    if gdata.invert==1
        zsos=-zsos;
    end
    
    gdata.bdata.zc=zsos-gdata.bdata.tide;
else
    gdata.bdata.zc=(gdata.bdata.zraw+gdata.manoff)-...
        gdata.bdata.tide;
end

%what about using ppk GPS
if gdata.applyppk==2
    if isfield(gdata.info,'ppk_filename')
        gdata.info=rmfield(gdata.info,'ppk_filename');
    end
    guidata(hfig,gdata)
    applyppk(hfig)
    gdata=guidata(hfig);
else
    if isfield(gdata.info,'ppk_filename')
        gdata.info=rmfield(gdata.info,'ppk_filename');
    end
    gdata.applyppk=0;
end



if isfield(gdata.bdata,'tide')
    set(gdata.wv,'enable','on')
    
else
    set(gdata.wv,'enable','off')
    
end

if gdata.echofix==1
    gdata.bdata=echofix(gdata.fixlen,gdata.bdata);
end


if isfield(gdata,'adistc');
    gdata=rmfield(gdata,'adistc');
end

figure(gdata.fh)
hold off
gdata.l1=plot(gdata.bdata.distance,gdata.bdata.zc);

fun=@(x)([nanmin(x)-(0.05*(nanmax(x)-nanmin(x))),...
    nanmax(x)+(0.05*(nanmax(x)-nanmin(x)))]);
gdata.xlimo=fun(gdata.bdata.distance);
gdata.ylimo=fun(gdata.bdata.zc);
gdata.xlims=fun(gdata.bdata.distance);
gdata.ylims=fun(gdata.bdata.zc);


set(gca,'xlim',gdata.xlimo,...
    'ylim',gdata.ylimo)

txtstr= ['Displaying transect ',...
    num2str(gdata.pointer),' of ',...
    num2str(length(gdata.fnames)),' - '...
    gdata.fnames{gdata.newInd}];

set(gdata.text1,'string',txtstr,'foregroundcolor','k');
guidata(hfig,gdata)

% catch %#ok
%     cla
%     txtstr= ['Problems loading transect ',...
%         num2str(gdata.pointer),' of ',...
%         num2str(length(gdata.fnames)),' - '...
%         gdata.fnames{gdata.newInd}];
%
%     set(gdata.text1,'string',txtstr,'foregroundcolor','r');
% end

rawfile=char(strcat(gdata.filepath,strtok(gdata.fnames(...
    gdata.newInd),'.'),'.bin'));
if exist(rawfile,'file');
    
    if ~isfield(gdata,'menu7')
        gdata.menu7=uimenu('label','View');
        gdata.menu8=uimenu(gdata.menu7,'label','Digitized Data',...
            'callback',@showCleanData,'visible','off');
        gdata.menu9=uimenu(gdata.menu7,'label','Raw Data',...
            'callback',@showRawData);
    else
        set(gdata.menu7,'visible','on');
        set(gdata.menu8,'visible','off');
        set(gdata.menu9,'visible','on');
    end
    set(gdata.chrawz,'visible','on')
else
    set(gdata.chrawz,'visible','off')
    if isfield(gdata,'menu7')
        set(gdata.menu9,'visible','off');
    end
    
end

if isfield(gdata,'eraw')
    gdata=rmfield(gdata,{'eraw';'ezc'});
end

gdata.alongflag=0;
gdata.bindata=[];
gdata.rawflag=0;
gdata.flag=[];
set(gdata.hp,'Title','Filter Controls');
set(gdata.list1,'visible','on');
set(gdata.text2,'visible','on');
set(gdata.edit2,'visible','on');
set(gdata.text3,'visible','on');
set(gdata.edit3,'visible','on');
set(gdata.text5,'visible','on');
set(gdata.edit4,'visible','on');
set(gdata.text4,'visible','on');
set(gdata.push5,'callback',@cleanProfile);
set(gdata.menu11,'callback',@cleanProfile);

set(gdata.rawe1,'visible','off')
set(gdata.rawt1,'visible','off')
set(gdata.rawe2,'visible','off')
set(gdata.rawt2,'visible','off')
set(gdata.rawe3,'visible','off')
set(gdata.rawt3,'visible','off')
set(gdata.rawp1,'visible','off')

set(gdata.menu15,'visible','off');
if ~isempty(gdata.bdata.adist)
    set(gdata.menu14,'visible','on')
end

gdata.numedits=0;
gdata.edits=[];
set(gdata.menu12,'visible','off')
set(gdata.menu13,'visible','off')
gdata.topo=[];
gdata.cdata=[];
guidata(hfig,gdata)
setFocus(hfig);


end

%%%%%---------------------------------------------------------------------
function xyzFilt(hfig, eventdata, handles) %#ok

gdata=guidata(hfig);
gdata.fType=get(gdata.list1,'value');
gdata.fDist=str2double(get(gdata.edit2,'string'));
gdata.fStrength=str2double(get(gdata.edit3,'string'));
gdata.wLen=str2double(get(gdata.edit4,'string'));
gdata.pan=get(gdata.toggle1,'value');
gdata.xlims=get(gca,'xlim');
gdata.ylims=get(gca,'ylim');
if gdata.pan==1;
    pan off
    set(gdata.toggle1,'value',0)
end

if isfield(gdata.bdata,'lon')
    geo=1;
else
    geo=0;
end

% ind=find(isfinite(gdata.bdata.zc));
% if length(ind)~=length(gdata.bdata.zc);
%     gdata.bdata.zc=gdata.bdata.zc(min(ind):max(ind));
%     gdata.bdata.x=gdata.bdata.x(min(ind):max(ind));
%     gdata.bdata.y=gdata.bdata.y(min(ind):max(ind));
%     gdata.bdata.distance=gdata.bdata.distance(min(ind):max(ind));
%     if ~isempty(gdata.bdata.adist)
%         gdata.bdata.adist=gdata.bdata.adist(min(ind):max(ind));
%     end
%     if isfield(gdata.bdata,'offline');
%         if ~isempty(gdata.bdata.offline)
%             gdata.bdata.offline=gdata.bdata.offline(min(ind):max(ind));
%         end
%     end
%     gdata.bdata.mtime=gdata.bdata.mtime(min(ind):max(ind));
%     if geo
%         gdata.bdata.lon=gdata.bdata.lon(min(ind):max(ind));
%         gdata.bdata.lat=gdata.bdata.lat(min(ind):max(ind));
%     end
%     if isfield(gdata.bdata,'tide')
%         gdata.bdata.tide=gdata.bdata.tide(min(ind):max(ind));
%     end
%     gdata.bdata.zraw=gdata.bdata.zraw(min(ind):max(ind));
%     gdata.xlimo=[min(gdata.bdata.distance),...
%         max(gdata.bdata.distance)];
%     gdata.ylimo=[min(gdata.bdata.zc),...
%         max(gdata.bdata.zc)];
% end



dist1=diff(gdata.bdata.distance);

%pre-filter data
dLen=length(gdata.bdata.zc);
fLen=round(gdata.fDist/nanmean(dist1(:)));
remLen=rem(dLen,fLen);

%group soundings into groups that are roughly "fDist" in width
zPF=reshape(gdata.bdata.zc(1:end-remLen),...
    [fLen floor(dLen/fLen)]);
% xPF=reshape(gdata.bdata.x(1:end-remLen),...
%     [fLen floor(dLen/fLen)]);
% yPF=reshape(gdata.bdata.y(1:end-remLen),...
%     [fLen floor(dLen/fLen)]);
% dnPF=reshape(gdata.bdata.mtime(1:end-remLen),...
%     [fLen floor(dLen/fLen)]);
% cDistPF=reshape(gdata.bdata.distance(1:end-remLen),...
%     [fLen floor(dLen/fLen)]);
% if ~isempty(gdata.bdata.adist)
%     caDistPF=reshape(gdata.bdata.adist(1:end-remLen),...
%         [fLen floor(dLen/fLen)]);
% end
% if geo
%     latPF=reshape(gdata.bdata.lat(1:end-remLen),...
%         [fLen floor(dLen/fLen)]);
%     lonPF=reshape(gdata.bdata.lon(1:end-remLen),...
%         [fLen floor(dLen/fLen)]);
% end

%detrend soundings in each group
zdt=detrend(zPF);
zdtm=nanmean(zdt);
zdtstd=nanstd(zdt);
[m,n]=size(zPF);


%loop through and remove outliers from each group

gdata.flag=zeros(m,n);
for i=1:n;
    ind=find(abs(zdt(:,i)-zdtm(i))>gdata.fStrength*zdtstd(i));
    zPF(ind,i)=NaN;
%     xPF(ind,i)=NaN;
%     yPF(ind,i)=NaN;
%     if geo
%         latPF(ind,i)=NaN;
%         lonPF(ind,i)=NaN;
%     end
%     dnPF(ind,i)=NaN;
end

% deal with the end of the profile
if remLen~=0
    ze=gdata.bdata.zc(end-remLen+1:end);
%     xe=gdata.bdata.x(end-remLen+1:end);
%     ye=gdata.bdata.y(end-remLen+1:end);
%     dc=gdata.bdata.distance(end-remLen+1:end);
%     if ~isempty(gdata.bdata.adist)
%         da=gdata.bdata.adist(end-remLen+1:end);
%         caDistPF=[caDistPF(:);da];
%     end
%     if geo
%         lone=gdata.bdata.lon(end-remLen+1:end);
%         late=gdata.bdata.lat(end-remLen+1:end);
%     end
%     dne=gdata.bdata.mtime(end-remLen+1:end);
    
    zdt=detrend(ze);
    zdtm=nanmean(zdt);
    zdtstd=nanstd(zdt);
    ze(abs(zdt-zdtm)>gdata.fStrength*zdtstd)=NaN;
    
    
    zPF=[zPF(:);ze];
%     xPF=[xPF(:);xe];
%     yPF=[yPF(:);ye];
%     if geo
%         lonPF=[lonPF(:);lone];
%         latPF=[latPF(:);late];
%     end
%     cDistPF=[cDistPF(:);dc(:)];
%     dnPF=[dnPF(:);dne];
end

gdata.flag=isnan(zPF(:));

switch gdata.fType
    case 1
        zc=slidefun(@nanmean,gdata.wLen,zPF(:));
    case 2
        zc=slidefun(@nanmedian,gdata.wLen,zPF(:));
    case 3
        zc=slidefun(@min,gdata.wLen,zPF(:));
    case 4
        zc=slidefun(@max,gdata.wLen,zPF(:));
end

% 
% 
% if geo
%     lonc=slidefun(@nanmean,gdata.wLen,lonPF(:));
%     latc=slidefun(@nanmean,gdata.wLen,latPF(:));
% end
% xc=slidefun(@nanmean,gdata.wLen,xPF(:));
% yc=slidefun(@nanmean,gdata.wLen,yPF(:));
% gdata.distc=slidefun(@nanmean,gdata.wLen,cDistPF(:));
% if ~isempty(gdata.bdata.adist)
%     gdata.adistc=slidefun(@nanmean,gdata.wLen,caDistPF(:));
% end
% dn=slidefun(@nanmean,gdata.wLen,dnPF(:));

% ind=find(isfinite(zc));
% ind1=min(ind);
% ind2=max(ind);
% xc=xc(ind1:ind2);
% yc=yc(ind1:ind2);
% zc=zc(ind1:ind2);
% if geo
%     lonc=lonc(ind1:ind2);
%     latc=latc(ind1:ind2);
% end
% dn=dn(ind1:ind2);
% gdata.distc=gdata.distc(ind1:ind2);
% if isfield(gdata,'adistc')
%     if ~isempty(gdata.adistc)
%         gdata.adistc=gdata.adistc(ind1:ind2);
%     end
% end



%fill gaps less than max value
gaps=isnan(zc(1:end-1));
if isempty(gaps)~=1;
    [c,ind]=getchunks(gaps,'-full');
    bind=find(isnan(zc(ind))==1);
    bstart=ind(bind)-1;
    bend=bsxfun(@plus,bstart,c(bind)+1);
    bend(bstart==0)=[];
    bstart(bstart==0)=[];
    
    dist=cellfun(@(x,y)(gdata.bdata.distance(y)-...
        gdata.bdata.distance(x)),...
        num2cell(bstart),num2cell(bend));
    
    fillGaps=find(dist<gdata.maxGap);
    fillStart=bstart(fillGaps);
    fillEnd=bend(fillGaps);
    
    
    for i=1:numel(fillGaps)
        if find(isnan([zc(fillStart(i));zc(fillEnd(i))]));
            zc(fillStart(i):fillEnd(i))=NaN;
        else
            zc(fillStart(i):fillEnd(i))=...
                interp1([gdata.bdata.distance(fillStart(i));...
                gdata.bdata.distance(fillEnd(i))],...
                [zc(fillStart(i));zc(fillEnd(i))],...
                gdata.bdata.distance(fillStart(i):fillEnd(i)));
        end
    end
    
end







hold off
if gdata.alongflag==1
    gdata.l1=plot(gdata.bdata.adist,gdata.bdata.zc);
    hold on
    gdata.bb=plot(gdata.bdata.adist(gdata.flag==1),...
        gdata.bdata.zc(gdata.flag==1),'o',...
        'color',[0.6 0.6 0.6],'markersize',3,...
        'markerfacecolor',[0.6 0.6 0.6]);
    gdata.gg=plot(gdata.bdata.adist,zc,'r-','linewidth',2);
else
    
    gdata.l1=plot(gdata.bdata.distance,gdata.bdata.zc);
    hold on
    gdata.bb=plot(gdata.bdata.distance(gdata.flag==1),...
        gdata.bdata.zc(gdata.flag==1),'o',...
        'color',[0.6 0.6 0.6],'markersize',3,...
        'markerfacecolor',[0.6 0.6 0.6]);
    gdata.gg=plot(gdata.bdata.distance,zc,'r-','linewidth',2);
end

set(gca,'xlim',gdata.xlims,'ylim',gdata.ylims);

c1=legend('Raw Data','Outliers','Filtered Data');
set(c1,'box','off','location','best')

if ~isempty(gdata.topo);
    gdata.tlineh=plot(gdata.topo.dist,gdata.topo.z,'g-',...
        'linewidth',2);
end


gdata.cdata.filename=gdata.bdata.filename;
gdata.cdata.xc=gdata.bdata.x;
gdata.cdata.yc=gdata.bdata.y;
gdata.cdata.zc=zc;
if geo
    gdata.cdata.lat=gdata.bdata.lat;
    gdata.cdata.lon=gdata.bdata.lon;
end
gdata.cdata.mtime=gdata.bdata.mtime;

guidata(hfig,gdata)

setFocus(hfig)
end

%%%%-----------------------------------------------------------------------
function addMeta(hfig,eventdata,handles) %#ok
gdata=guidata(hfig);

[filename, pathname] = uigetfile( ...
    {'*.meta', 'META Files (*.meta)'},...
    'Select a Meta file',gdata.metapath);


if filename==0
    return
end

set(gdata.text1,'string',['Reading File: ',filename],...
    'foregroundcolor','r')
drawnow

%read the meta file
[mdata,meta]=read_meta([pathname, filename]);
metatime=cellstr(datestr([mdata.Year,mdata.Month,...
    mdata.Day,mdata.Hour,mdata.Minute,mdata.Second],31));
%
% if ~isfield(gdata,'sos_corr');
%     gdata.sos_corr.ratio=1;
% else
%     gdata.sos_corr.sos_old=str2double(meta.speed_of_sound);
% end
%
%
%
%
% gdata.info.vessel=meta.vessel_id;
% gdata.info.vessel_op=meta.vessel_operator;
% gdata.info.echosounder=meta.echosounder_type;
% if isfield(meta,'speed_of_sound')
%     gdata.info.speed_of_sound=meta.speed_of_sound;
% end
% if isfield(meta,'base_stn_id')
%     gdata.info.base_id=meta.base_stn_id;
% else
%     gdata.info.base_id='unknown';
% end
% gdata.info.analyst=meta.data_analyst;

% limit the bdata to the times when data from meta file are available
% based on time strings in respective files
rawtime=datestr(gdata.bdata.mtime,31);

if numel(gdata.bdata.zraw)~=numel(mdata.Depth)
    
    [ia,ib]=ismember(rawtime,metatime);
    ib(~ia)=[];
    if isempty(ib);
        errordlg({'Meta File does not match current .RAW file.';...
            'Residing data retained'})
        txtstr= ['Displaying transect ',...
            num2str(gdata.pointer),' of ',...
            num2str(length(gdata.fnames)),' - '...
            gdata.fnames{gdata.newInd}];
        set(gdata.text1,'string',txtstr,'foregroundcolor','k');
        
        return
    end
    
    
    [~,ib2]=ismember(metatime,rawtime);
    
    fp=(ib2(1)-ib(1))+1;
    idx=fp:length(metatime)+(fp-1);
    
    
    fields=fieldnames(gdata.bdata);
    for i = 2:length(fields);
        if ~isempty(gdata.bdata.(fields{i}))
            gdata.bdata.(fields{i})=gdata.bdata.(fields{i})(idx);
        end
    end
end

% replace raw depth, and corrected depth in bdata;
gdata.bdata.zraw=mdata.Depth;

cla
% if gdata.alongflag==1;
%     gdata.l1=plot(gdata.bdata.adist,gdata.bdata.zc);
% else
%     gdata.l1=plot(gdata.bdata.distance,gdata.bdata.zc);
% end

gdata.metapath=pathname;

txtstr= ['Displaying transect ',...
    num2str(gdata.pointer),' of ',...
    num2str(length(gdata.fnames)),' - '...
    gdata.fnames{gdata.newInd}];
set(gdata.text1,'string',txtstr,'foregroundcolor','k');

guidata(hfig,gdata);


end
%%%%-----------------------------------------------------------------------
function addTopo(hfig,eventdata,handles) %#ok
gdata=guidata(hfig);

[filename, pathname] = uigetfile( ...
    {'*.xyz', 'XYZ Files (*.xyz)'},...
    'Select a Topo file',gdata.topopath);

if ~isequal (filename,0)
    % if it is a straight line followed, calculate offline distance
    xyz=load([pathname,filename]);
    
    if isfield(gdata.hdr,'linex');
        m=(gdata.hdr.liney(end)-gdata.hdr.liney(1))/...
            (gdata.hdr.linex(end)-gdata.hdr.linex(1));
        
        %determine along-line distance
        
        theta=atan(m);
        xt = xyz(:,1).*cos(theta) + xyz(:,2).*sin(theta);
        yt = -xyz(:,1).*sin(theta) + xyz(:,2).*cos(theta);
        x0=gdata.hdr.linex.*cos(theta) + ...
            gdata.hdr.liney.*sin(theta);
        toff=detrend(yt,'constant');
        
        %calculate along line distance
        gdata.topo.dist=xt-x0(1);
        gdata.topo.offline=toff;
        gdata.topo.z=xyz(:,3);
        
        
        gdata.xlimo=[min([gdata.bdata.adist;gdata.topo.dist]),...
            max([gdata.bdata.adist;gdata.topo.dist])];
        gdata.ylimo=[min([gdata.bdata.zc;gdata.topo.z]),...
            max([gdata.bdata.zc;gdata.topo.z])];
        gdata.xlims=[min([gdata.bdata.adist;gdata.topo.dist]),...
            max([gdata.bdata.adist;gdata.topo.dist])];
        gdata.ylims=[min([gdata.bdata.zc;gdata.topo.z]),...
            max([gdata.bdata.zc;gdata.topo.z])];
        
        
    end
    
    viewAlong(hfig);
    
    hold on
    gdata.tlineh=plot(gdata.topo.dist,gdata.topo.z,...
        'g-','linewidth',2);
    set(gca,'xlim',gdata.xlims,'ylim',gdata.ylims)
    
    set(get(gca,'ylabel'),'string','Elevation (m)',...
        'fontsize',12)
    
    gdata.alongflag=1;
    gdata.topopath=pathname;
    guidata(hfig,gdata)
end

end

%%%%-----------------------------------------------------------------------
function setFocus(hfig,eventdata,handles) %#ok
gdata=guidata(hfig);

pos = get(gcf, 'Position');  % User might have moved figure.
pointpos = get(0,'pointerlocation');  % Current pointer location.
set(0, 'PointerLocation', [pos(1)+(pos(3)/2),pos(2) + (pos(4)/2)]);
% Now we simulate a mouseclick on the figure using the JAVA.
gdata.ja.mousePress(java.awt.event.InputEvent.BUTTON1_MASK);  % Click down
gdata.ja.mouseRelease(java.awt.event.InputEvent.BUTTON1_MASK); % Let up.
set(0,'pointerlocation',pointpos);  % Put the pointer back.
pause(.025)   % drawnow does NOT work here.

end
%%%%-----------------------------------------------------------------------
function zoomto(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);
set(gdata.push3,'backgroundcolor','g');

gdata.pan=get(gdata.toggle1,'value');
if gdata.pan==1;
    set(gdata.toggle1,'value',0);
    pan off
end

k=waitforbuttonpress; %#ok
point1 = get(gca,'CurrentPoint');
finalRect = rbbox; %#ok
point2 = get(gca,'CurrentPoint');
point1 = point1(1,1:2);
point2 = point2(1,1:2);

xlimsN=sort([point1(1),point2(1)]);
ylimsN=sort([point1(2),point2(2)]);
set(gdata.ax1,'xlim',xlimsN,'ylim',ylimsN);

gdata.xlims=xlimsN;
gdata.ylims=ylimsN;
if xlimsN(1)>=1;
    gdata.pno=round(xlimsN(1));
else
    gdata.pno=1;
end

if gdata.rawflag
    xdata=1:numel(gdata.bindata.vtime);
    ydata=gdata.bindata.range;
    
    xind=floor(gdata.xlims(1)):ceil(gdata.xlims(2));
    yind=find(ydata>=gdata.ylims(1) & ydata<=gdata.ylims(2));
    
    ia=intersect(xind,xdata);
    
    if all(~isempty([ia(:);yind(:)]));
        
        gdata.binzoom.x=xdata(ia);
        gdata.binzoom.y=ydata(yind);
        gdata.binzoom.cdata=gdata.bindata.vals(yind,ia);
        set(gdata.rim,'xdata',gdata.binzoom.x,...
            'ydata',gdata.binzoom.y,...
            'cdata',gdata.binzoom.cdata);
        guidata(hfig,gdata)
    end
    
    
end


set(gdata.push3,'backgroundcolor',[0.9255    0.9137    0.8471]);
guidata(hfig,gdata);

setFocus(hfig)

end

%%%%%----------------------------------------------------------------------
function [] = batch_add_meta_gui(hf,evnt); %#ok
% The basic layout of this GUI was made with the help of guidegetter,
% available on the File Exchange at Mathworks.com

gd=guidata(hf);
gd2.metapath=gd.filepath;
gd2.ncpath=gd.filepath;

hf2 = figure('units','normalized',...
    'position',[0.332 0.639 0.153 0.28],...
    'menubar','none','name','Add/modify Metadata',...
    'numbertitle','off','color',[0.94 0.94 0.94]);

uicontrol(hf2,'style','pushbutton',...
    'units','normalized',...
    'position',[0.0806 0.76 0.253 0.0975],...
    'string','Open Meta Info',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'callback',@select_meta_txt);
gd2.mtext = uicontrol(hf2,'style','text',...
    'units','normalized',...
    'position',[0.392 0.783 0.509 0.0613],...
    'string','No File Selected',...
    'backgroundcolor',[0.94 0.94 0.94]);

gd2.nclist = uicontrol(hf2,'style','listbox',...
    'units','normalized',...
    'position',[0.396 0.203 0.52 0.515],...
    'string','No files selected',...
    'backgroundcolor',[1 1 1]);
uicontrol(hf2,'style','pushbutton',...
    'units','normalized',...
    'position',[0.0842 0.621 0.253 0.0975],...
    'string','Open .nc files',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'callback',@select_nc_files);


uicontrol(hf2,'style','pushbutton',...
    'units','normalized',...
    'position',[0.363 0.0223 0.253 0.0919],...
    'string','Cancel',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'callback',@(h,e)(close(hf2)));


gd2.done = uicontrol(hf2,'style','pushbutton',...
    'units','normalized',...
    'position',[0.663 0.0223 0.253 0.0919],...
    'string','Done',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'enable','off',...
    'callback',@add_meta_nc);

gd2.goflag=[0 0];
guidata(hf2,gd2)
end

function select_meta_txt(hf2,evnt) %#ok

gd2=guidata(hf2);
[filename, pathname] = uigetfile( ...
    {'*.*', 'All Files (*.*)'},...
    'Select a File',gd2.metapath);
if filename==0
    return
    
end


gd2.metafile=filename;
gd2.metapath=pathname;
gd2.goflag(1)=1;

set(gd2.mtext,'string',gd2.metafile)

if all(gd2.goflag==1)
    set(gd2.done,'enable','on')
end

guidata(hf2,gd2)
end

function select_nc_files(hf2,evnt) %#ok

gd2=guidata(hf2);
[filename, pathname] = uigetfile( ...
    {'*.nc', 'NC Files (*.nc)'},...
    'Select a File(s)','multiselect','on',...
    gd2.ncpath);
if pathname==0
    return
    
end

if ~iscell(filename)
    filename={filename};
end;
gd2.ncfiles=filename;
gd2.ncpath=pathname;
gd2.goflag(2)=1;


set(gd2.nclist,'string',gd2.ncfiles)

if all(gd2.goflag==1)
    set(gd2.done,'enable','on')
end

guidata(hf2,gd2)
end

function add_meta_nc(hf2,evnt) %#ok

gd2=guidata(hf2);

%read the metadata text file
ini=read_info([gd2.metapath,gd2.metafile]);
natts=fieldnames(ini); %new atts

for i=1:length(gd2.ncfiles)
    for j=1:length(natts)
        ncwriteatt([gd2.ncpath,gd2.ncfiles{i}],'/',...
            natts{j},ini.(natts{j}))
    end
    
end

set(gd2.nclist,'string','No files selected')
gd2.goflag(2)=0;
set(gd2.done,'enable','off')

gd2.ncpath=[];
gd2.ncfiles=[];
guidata(hf2,gd2);
end
%%%%-----------------------------------------------------------------------
function [] = batch_add_offset(hf,evnt); %#ok
% The basic layout of this GUI was made with the help of guidegetter,
% available on the File Exchange at Mathworks.com

gd=guidata(hf);
gd2.metapath=gd.filepath;
gd2.ncpath=gd.filepath;

hf2 = figure('units','normalized',...
    'position',[0.332 0.639 0.153 0.28],...
    'menubar','none','name','Batch Add Manual Offset',...
    'numbertitle','off','color',[0.94 0.94 0.94]);

uicontrol(hf2,'style','Text',...
    'units','normalized',...
    'position',[0.0806 0.86 0.253 0.0975],...
    'string','Manual Offset',...
    'backgroundcolor',[0.94 0.94 0.94]);
gd2.mtext = uicontrol(hf2,'style','edit',...
    'units','normalized',...
    'position',[0.392 0.883 0.509 0.0613],...
    'string',gd.manoff,...
    'backgroundcolor',[0.94 0.94 0.94],...
    'callback',@checkoffset);
gd2.checkbox=uicontrol(hf2,'style','Checkbox',...
    'units','normalized',...
    'position',[0.0806 0.76 0.53 0.0975],...
    'string','Write XYZ Files',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'value',0);
gd2.nclist = uicontrol(hf2,'style','listbox',...
    'units','normalized',...
    'position',[0.396 0.203 0.52 0.515],...
    'string','No files selected',...
    'backgroundcolor',[1 1 1]);
uicontrol(hf2,'style','pushbutton',...
    'units','normalized',...
    'position',[0.0842 0.621 0.253 0.0975],...
    'string','Open .nc files',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'callback',@select_nc_files_offset);

gd2.textstr=uicontrol(hf2,'style','text',...
        'units','normalized',...
    'position',[0.0806 0.1 0.453 0.0919],...
    'string','Working Please Wait...',...
    'foregroundcolor','r',...
    'visible','off');

uicontrol(hf2,'style','pushbutton',...
    'units','normalized',...
    'position',[0.363 0.0223 0.253 0.0919],...
    'string','Cancel',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'callback',@(h,e)(close(hf2)));


gd2.done = uicontrol(hf2,'style','pushbutton',...
    'units','normalized',...
    'position',[0.663 0.0223 0.253 0.0919],...
    'string','Done',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'enable','off',...
    'callback',@add_offset_nc);

gd2.goflag=[0 0];
guidata(hf2,gd2)
end

function checkoffset(hf2,evnt) %#ok

gd2=guidata(hf2);
offset=str2double(get(gd2.mtext,'string'));
if isfinite(offset)
    gd2.goflag(1)=1;
else
    gd2.goflag(1)=0;
end

if all(gd2.goflag==1)
    set(gd2.done,'enable','on')
else 
    set(gd2.done,'enable','off')
end

guidata(hf2,gd2)
end

function select_nc_files_offset(hf2,evnt) %#ok

gd2=guidata(hf2);
[filename, pathname] = uigetfile( ...
    {'*.nc', 'NC Files (*.nc)'},...
    'Select a File(s)','multiselect','on',...
    gd2.ncpath);
if pathname==0
    return
    
end

if ~iscell(filename)
    filename={filename};
end;
gd2.ncfiles=filename;
gd2.ncpath=pathname;
gd2.goflag(2)=1;


set(gd2.nclist,'string',gd2.ncfiles)

if all(gd2.goflag==1)
    set(gd2.done,'enable','on')
end

guidata(hf2,gd2)
end


function add_offset_nc(hf2,evnt) %#ok
%add_manual_offset_nc -add a vertical offset to netcdf file
%   applies offset to 'zc' and 'zf' (corrected and smoothed fields in the
%   bathy group. offset is added (positive makes profile shallower).


gd2=guidata(hf2);
set(gd2.textstr,'visible','on')
drawnow
offset=str2double(get(gd2.mtext,'string'));
writexyz=get(gd2.checkbox,'value');

for i=1:length(gd2.ncfiles)
    info=ncinfo([gd2.ncpath,gd2.ncfiles{i}]);
    groups=arrayfun(@(x)(x.Name),info.Groups,'un',0);
    gidx=find(strcmpi('bathy',groups));
    if ~isempty(gidx)
        bfields=arrayfun(@(x)(x.Name),info.Groups(gidx).Variables,'un',0);
        [cfields,ci]=intersect({'zc';'zf'},bfields);
        if ~isempty(cfields)
            %first see if there already is a manual offset applied
            gnames=arrayfun(@(x)(x.Name),info.Attributes,'un',0)';
            
            if any(strcmpi('manual_offset',gnames))
                orig_off=ncreadatt([gd2.ncpath,gd2.ncfiles{i}],...
                    '/','manual_offset');
                fin_off=str2double(orig_off)+offset;
            else
                fin_off=offset;
            end
            
            ncwriteatt([gd2.ncpath,gd2.ncfiles{i}],'/',...
                'manual_offset',fin_off);
            str=sprintf('%0.3f',fin_off);
            ncwriteatt([gd2.ncpath,gd2.ncfiles{i}],'/',...
                'manual_offset',str);
            
            
            for j = 1:length(ci)
                zn=ncread([gd2.ncpath,gd2.ncfiles{i}],...
                    ['bathy/',cfields{j}]);
                ncwrite([gd2.ncpath,gd2.ncfiles{i}],...
                    ['bathy/',cfields{j}],zn+offset);
            end
            
            if writexyz
                
                x=ncread([gd2.ncpath,gd2.ncfiles{i}],...
                    'bathy/x');
                y=ncread([gd2.ncpath,gd2.ncfiles{i}],...
                    'bathy/y');
                dlmwrite([gd2.ncpath,strtok(gd2.ncfiles{i},'.'),...
                    '.xyz'],[x y zn+offset],...
                    'delimiter',' ' ,'precision','%0.3f');
            end
        end
    end
end
  

set(gd2.textstr,'visible','off')
set(gd2.nclist,'string','No files selected')
gd2.goflag(2)=0;
set(gd2.done,'enable','off')

gd2.ncpath=[];
gd2.ncfiles=[];
guidata(hf2,gd2);

end
%%%%%----------------------------------------------------------------------
function panIm(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);


if gdata.rawflag
    
    set(gdata.rim,'xdata',(1:numel(gdata.bindata.vtime)),...
        'ydata',gdata.bindata.range,...
        'cdata',gdata.bindata.vals);
end


gdata.pan=get(gdata.toggle1,'value');

if gdata.pan==1;
    set(gdata.toggle1,'backgroundcolor','g')
    pan on
else
    set(gdata.toggle1,'backgroundcolor',[ 0.9255    0.9137    0.8471])
    pan off
end




guidata(hfig,gdata);
setFocus(hfig)
end

%%%%----------------------------------------------
function fullExtents(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);
gdata.pan=get(gdata.toggle1,'value');

if gdata.pan==1;
    set(gdata.toggle1,'value',0);
    pan off
end

if gdata.rawflag~=1
    set(gca,'xlim',gdata.xlimo)
    set(gca,'ylim',gdata.ylimo)
else
    set(gca,'xlim',gdata.xlimr)
    set(gca,'ylim',gdata.ylimr)
    set(gdata.rim,'xdata',(1:numel(gdata.bindata.vtime)),...
        'ydata',gdata.bindata.range,...
        'cdata',gdata.bindata.vals);
end

guidata(hfig,gdata);

setFocus(hfig);
end

%%%%----------------------------------------------
function export1(hfig, eventdata, handles)%#ok

gdata=guidata(hfig);

assignin('base','bdata',gdata.bdata);
assignin('base','rdata',gdata.rdata);
assignin('base','hdr',gdata.hdr);

if ~isempty(gdata.bindata)
    assignin('base','bindata',gdata.bindata)
end


if ~isempty(gdata.info)
    assignin('base','info',gdata.info)
end

if isempty(gdata.cdata)~=1
    assignin('base','fdata',gdata.cdata)
end

if isfield(gdata,'ppk')
    assignin('base','ppk_data',gdata.ppk)
end


if isfield(gdata,'wvs')
    if isempty(gdata.wvs)~=1
        assignin('base','wvsData',gdata.wvs)
    end
end

end

%%%%----------------------------------------------
function export2(hfig, eventdata, handles)%#ok

gdata=guidata(hfig);

namer=strtok(gdata.bdata.filename,'.');

[filename, pathname] = uiputfile( ...
    {'*.xyz', 'XYZ Files'}, ...
    'Save as',[gdata.outpath,namer,'.xyz']);

if filename==0
    return
end

if isempty(gdata.cdata)==1
    ind=find(isfinite(gdata.bdata.zc));
    ind1=min(ind);
    ind2=max(ind);
    xr=gdata.bdata.x(ind1:ind2);
    yr=gdata.bdata.y(ind1:ind2);
    zr=gdata.bdata.zc(ind1:ind2);
    dlmwrite([pathname,filename],[xr yr zr],'precision','%.3f',...
        'delimiter','\t')
else
    ind=find(isfinite(gdata.cdata.zc));
    ind1=min(ind);
    ind2=max(ind);
    xr=gdata.cdata.xc(ind1:ind2);
    yr=gdata.cdata.yc(ind1:ind2);
    zr=gdata.cdata.zc(ind1:ind2);
    dlmwrite([pathname,filename],[xr yr zr],'precision','%.3f',...
        'delimiter','\t')
end

gdata.outpath=pathname;
if gdata.outpath(end)~=filesep;
    gdata.outpath=[gdata.outpath,filesep];
end

guidata(hfig,gdata)
end

%%%%-----------------------------------------------------------------------
function gesetup(hfig,eventdata,handles) %#ok
gdata=guidata(hfig);

prompt={'Thin:';'Color Minimum:';'Color Maximum:';...
    'Symbol Scale (0-1):';'Colormap:'};
name='GE Export Setup';
numlines=1;
format={'%.0f';'%.2f';'%.2f';'%.1f';'%s'};

defaultanswer=cellfun(@(x,y)(sprintf(x,y)),...
    format,struct2cell(gdata.ge),'uni',0);

answer=inputdlg(prompt,name,numlines,defaultanswer);

gdata.ge.thin=str2double(answer{1});
gdata.ge.cmin=str2double(answer{2});
gdata.ge.cmax=str2double(answer{3});
gdata.ge.scale=str2double(answer{4});
gdata.ge.cmap=answer{5};
guidata(hfig,gdata);
end

%%%%%%---------------------------------------------------------------------
function metainput(hfig,event,handles) %#ok
% now defunct-06/26/2014, remove soon
gdata=guidata(hfig);

prompt={'Vessel ID:';'Vessel Operator';....
    'Echosounder Type:';'Speed of Sound (m/s):';...
    'Base Station ID:';'Analyst:'};
name='Metadata input';
numlines=1;
format={'%s';'%s';'%s';'%s';'%s';'%s'};

defaultanswer=cellfun(@(x,y)(sprintf(x,y)),...
    format,struct2cell(gdata.info),'uni',0);

answer=inputdlg(prompt,name,numlines,defaultanswer);

if isempty(answer)
    return
end


gdata.info.vessel=answer{1};
gdata.info.vessel_op=answer{2};
gdata.info.echosounder=answer{3};
gdata.info.speed_of_sound=answer{4};
gdata.info.base_id=answer{5};
gdata.info.analyst=answer{6};
guidata(hfig,gdata);
end

%%%%%--------------------------------------------------------------------
function export3(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);
if isempty(gdata.cdata)==1
    gxyz=gdata.bdata;
else
    gxyz=gdata.cdata;
end

namer=strtok(gdata.bdata.filename,'.');
[filename, pathname] = uiputfile('*.kml', 'Save profile as',...
    [gdata.outpath1,namer]);

if filename==0
    return
end

ind=find(isfinite(gxyz.zc) & isfinite(gxyz.lon) & ...
    isfinite(gxyz.lat) & isfinite(gxyz.mtime));
lon=gxyz.lon(ind(1:gdata.ge.thin:end));
lat=gxyz.lat(ind(1:gdata.ge.thin:end));
z=gxyz.zc(ind(1:gdata.ge.thin:end));
mtime=gxyz.mtime(ind(1:gdata.ge.thin:end));

gescatter([pathname,filename],lon,lat,z,'time',mtime,...
    'clims',[gdata.ge.cmin gdata.ge.cmax],'colormap',gdata.ge.cmap,...
    'scale',gdata.ge.scale)

gdata.outpath1=pathname;
if gdata.outpath1(end)~=filesep;
    gdata.outpath1=[gdata.outpath1,filesep];
end

guidata(hfig,gdata)
end

%%%%%%---------------------------------------------------------------------
function export4(hfig, eventdata, handles)%#ok
%more defunct output options
gdata=guidata(hfig);

namer=strtok(gdata.bdata.filename,'.');

[filename, pathname] = uiputfile( ...
    {'*.txt', 'TXT Files'}, ...
    'Save as',[gdata.outpath,namer,'.txt']);

if filename==0
    return
end

h = waitbar(0,'Creating Infobank file, Please wait...');
set(h,'name','Creating Infobank File');

if isempty(gdata.cdata)==1
    ind=find(isfinite(gdata.bdata.zc)==1 & ...
        isfinite(gdata.bdata.mtime));
    xr=gdata.bdata.lon(ind);
    yr=gdata.bdata.lat(ind);
    zr=gdata.bdata.zc(ind);
    dn=gdata.bdata.mtime(ind);
else
    ind=find(isfinite(gdata.cdata.zc)==1 & ...
        isfinite(gdata.cdata.mtime));
    xr=gdata.cdata.lon(ind);
    yr=gdata.cdata.lat(ind);
    zr=gdata.cdata.zc(ind);
    dn=gdata.cdata.mtime(ind);
    
end

% dn1=dn(1);
% diffdn=nanmean(diff(dn));
% dnvec=dn1+((1:numel(dn))*diffdn)';

dni=dn(1):(1/86400):max(dn);
[dnu,ui]=unique(dn);
xc=interp1(dnu,xr(ui),dni);
yc=interp1(dnu,yr(ui),dni);
zc=interp1(dnu,zr(ui),dni);


[year,month,day,hr,mini,sec]=datevecfix(dni,'precision',2);
doy=(datenum(year,month,day)-datenum(year,1,1))+1;
secr=floor(sec);
tsec=floor((sec-secr).*10);


fid=fopen([pathname,filename],'wt');
for i=1:length(xc)
    fprintf(fid,'%d%0.3d%0.2d%0.2d%0.2d%d\t',...
        year(i),doy(i),hr(i),mini(i),secr(i),tsec(i));
    fprintf(fid,'%.6f\t',yc(i));
    fprintf(fid,'%.6f\t',xc(i));
    fprintf(fid,'%.2f\n',zc(i));
    %     waitbar(i/length(xc),h,sprintf('%d%% complete...',...
    %         round((i/length(xc))*100)));
end
fclose(fid);
close(h);

gdata.outpath=pathname;
if gdata.outpath(end)~=filesep;
    gdata.outpath=[gdata.outpath,filesep];
end
guidata(hfig,gdata)

end

%%%%%----------------------------------------------------------------------
function export5(hfig, eventdata, handles,varargin)%#ok
%defunct output format, no longer available to write
gdata=guidata(hfig);

if isempty(varargin)
    namer=strtok(gdata.bdata.filename,'.');
    [filename, pathname] = uiputfile( ...
        {'*.meta', 'META Files'}, ...
        'Save as',[gdata.outpath,namer,'.txt']);
    
    if filename==0
        return
    end
else
    pathname=gdata.outpath;
    filename=varargin{1};
end

h = waitbar(0,'Creating Meta file, Please wait...');
set(h,'name','Creating Meta File');

fid=fopen([pathname,filename],'wt');
%write the file header
fprintf(fid,'Meta File Version:\t%0.2f\n',1.20);
fprintf(fid,'Raw Filename:\t%s\n',gdata.bdata.filename);
fprintf(fid,'Vessel ID:\t%s\n',gdata.info.vessel);
fprintf(fid,'Vessel Operator:\t%s\n',gdata.info.vessel_op);
fprintf(fid,'Echosounder Type:\t%s\n',gdata.info.echosounder);
if isfield(gdata,'xducer_settings');
    fprintf(fid,'Transducer Frequency (kHz):\t%0.0f\n',...
        gdata.xducer_settings.Channel_1_frequency./10);
    fprintf(fid,'Transducer Blanking (m):\t%0.1f\n',...
        gdata.xducer_settings.Blanking./10);
    fprintf(fid,'Transducer Gain:\t%0.0f\n',...
        gdata.xducer_settings.Channel_1_Gain);
    fprintf(fid,'Transducer Transmit Power:\t%0.0f\n',...
        gdata.xducer_settings.Channel_1_TxPower);
end

fprintf(fid,'Speed of Sound (m/s):\t%s\n',gdata.info.speed_of_sound);
if isfield(gdata.rdata,'antennaH')
    fprintf(fid,'GPS Antenna Offset (m):\t%0.2f\n',...
        -mode(gdata.rdata.antennaH+gdata.rdata.tide));
    fprintf(fid,'Manual Elevation Offset (m):\t%0.2f\n',...
        gdata.manoff);
end

if gdata.applyppk~=0
    fprintf(fid,'Post-Process GPS Filename:\t%s\n',...
        gdata.ppk_info.filename);
end


fprintf(fid,'Base Station ID: \t%s\n',gdata.info.base_id);
fprintf(fid,'Data Analyst:\t%s\n',gdata.info.analyst);
fprintf(fid,'Transect Viewer Version:\t%0.2f\n\n',gdata.tv_ver);

%write column headers
fprintf(fid,'Date and Time (yyyy-mm-dd HH:MM:SS)\t');
fprintf(fid,'Easting (m)\tNorthing (m)\tDepth (m)');

if ~isempty(gdata.cdata)
    if numel(gdata.bdata.mtime)~=numel(gdata.cdata.zc)
        ind=find(gdata.bdata.mtime>=gdata.cdata.mtime(1) & ...
            gdata.bdata.mtime<=gdata.cdata.mtime(end));
    else
        ind=1:numel(gdata.bdata.mtime);
    end
else
    ind=1:numel(gdata.bdata.mtime);
end


%determine which types of data are available
dataout=cell(4,1);
dataout{1}=datestr(gdata.bdata.mtime(ind),31);
dataout{2}=gdata.bdata.x(ind);
dataout{3}=gdata.bdata.y(ind);
dataout{4}=gdata.bdata.zraw(ind);

dfmt={'%s\t','%0.2f\t','%0.2f\t','%0.2f\t'};


if isfield(gdata.bdata,'tide');
    fprintf(fid,'\tTide Correction (m)\tElevation (m)');
    nd=length(dataout);
    dataout=cat(1,dataout{:},cell(2,1));
    dataout{nd+1}=gdata.bdata.tide(ind);
    dataout{nd+2}=gdata.bdata.zc(ind);
    dfmt=cat(1,dfmt{:},{'%0.2f\t';'%0.2f\t'});
    
end

if ~isempty(gdata.cdata)
    fprintf(fid,'\tSmoothed Elevation (m)');
    nd=length(dataout);
    dataout{nd+1}=gdata.cdata.zc;
    dfmt=cat(1,dfmt{:},{'%0.2f\t'});
    
end

if isfield(gdata.bdata,'adist')
    if ~isempty(gdata.bdata.adist)
        fprintf(fid,'\tOffline (m)\tDistance (m)');
        nd=length(dataout);
        dataout=cat(1,dataout{:},cell(2,1));
        dataout{nd+1}=gdata.bdata.offline(ind);
        dataout{nd+2}=gdata.bdata.adist(ind);
        dfmt=cat(1,dfmt{:},{'%0.2f\t';'%0.2f\t'});
    end
else
    fprintf(fid,'\tDistance (m)');
    nd=length(dataout);
    dataout=cat(1,dataout{:},cell(1,1));
    dataout{nd+1}=gdata.bdata.distance(ind);
    dfmt=cat(1,dfmt{:},{'%0.2f\t'});
end
dfmt{end}=[dfmt{end},'\n'];
fprintf(fid,'\n');


for i =1:length(ind)
    for j=1:length(dataout)
        fprintf(fid,dfmt{j},dataout{j}(i,:));
    end
end

waitbar(1,h,'Done!');
fclose(fid);
close(h);

gdata.outpath=pathname;
if gdata.outpath(end)~=filesep;
    gdata.outpath=[gdata.outpath,filesep];
end
guidata(hfig,gdata)

end

%%%%-----------------------------------------------------------------------
function export6(hfig,evnt,varargin) %#ok
gdata=guidata(hfig);

if isempty(varargin)
    
    namer=strtok(gdata.bdata.filename,'.');
    [filename, pathname] = uiputfile( ...
        {'*.nc', 'netCDF Files'}, ...
        'Save as',[gdata.outpath,namer,'.nc']);
    
    if filename==0
        return
    else
        gdata.outpath=pathname;
        guidata(hfig,gdata);
        opt.outfile=[pathname,filename];
    end
else
    opt.outfile=varargin{1};
end


opt.fpath=gdata.filepath;
opt.bdata=gdata.bdata;
opt.rdata=gdata.rdata;
opt.hdr=gdata.hdr;
if isfield(gdata,'info')
    if ~isempty(gdata.info)
        opt.info=gdata.info;
    end
end
if isfield(gdata,'bindata')
    if ~isempty(gdata.bindata)
        opt.bindata=gdata.bindata;
        if all(sign(gdata.bindata.range))
            opt.bindata.range=-gdata.bindata.range;
        end
    end
end

%smoothed data present
if ~isempty(gdata.cdata)==1
    type=get(gdata.list1,'value');
    switch type
        case 1
            opt.fdata.ftype='mean';
        case 2
            opt.fdata.ftype='median';
        case 3
            opt.fdata.ftype='min';
        case 4
            opt.fdata.ftype='max';
    end
    opt.fdata.fdist=str2double(get(gdata.edit2,'string'));
    opt.fdata.fstrength=str2double(get(gdata.edit3,'string'));
    opt.fdata.wlen=str2double(get(gdata.edit4,'string'));
    opt.fdata.maxgap=gdata.maxGap;
    opt.fdata.zf=gdata.cdata.zc;
end

if isfield(gdata,'sv')
    opt.sv=gdata.sv;
end

raw2nc(opt);




end
%%%%%----------------------------------------------------------------------
function batch_opt = tv_batch_gui(hfig,evnt)%#ok

gd=guidata(hfig);

hf = figure('units','normalized','position',[0.253 0.387 0.134 0.238],...
    'menubar','none','name','Configure Batch Output',...
    'numbertitle','off','color',[0.925 0.914 0.847]);

bod.check5 = uicontrol(hf,'style','checkbox','units','normalized',...
    'position',[0.178 0.743 0.615 0.0888],'string','XYZ File (.xyz)',...
    'backgroundcolor',[0.925 0.914 0.847],'value',gd.batch.out_xyz);
bod.check4 = uicontrol(hf,'style','checkbox','units','normalized',...
    'position',[0.178 0.612 0.615 0.0888],...
    'string','NetCDF File (.nc)',...
    'backgroundcolor',[0.925 0.914 0.847],'value',gd.batch.out_nc);
bod.check3 = uicontrol(hf,'style','checkbox','units','normalized',...
    'position',[0.178 0.48 0.615 0.0888],...
    'string','Google Earth File (.kml)',...
    'backgroundcolor',[0.925 0.914 0.847],'value',gd.batch.out_kml);
bod.check2 = uicontrol(hf,'style','checkbox','units','normalized',...
    'position',[0.178 0.349 0.615 0.0888],'string','Shape FIle (.shp)',...
    'backgroundcolor',[0.925 0.914 0.847],'value',gd.batch.out_shp);
uicontrol(hf,'style','pushbutton','units','normalized',...
    'position',[0.615 0.0493 0.295 0.0822],'string','Done',...
    'backgroundcolor',[0.925 0.914 0.847],'callback',@bodclose);

guidata(hf,bod);

uiwait
bod=guidata(hf);
gd.batch=bod.batch;
guidata(hfig,gd);

close(hf)
end


function bodclose(hf,evnt) %#ok

bod=guidata(hf);
bod.batch.out_xyz=get(bod.check5,'value');
bod.batch.out_nc=get(bod.check4,'value');
bod.batch.out_kml=get(bod.check3,'value');
bod.batch.out_shp=get(bod.check2,'value');


guidata(hf,bod);
uiresume
end

% %%%%---------------------------------------------------------------------
function exportAll(hfig, eventdata, handles)%#ok

gdata=guidata(hfig);

namer=strtok(gdata.bdata.filename,'.');

[filename, pathname] = uiputfile( ...
    {'*.*',  'All Files (*.*)'},...
    'Select File Name',[gdata.outpath,namer]);

if filename==0
    return
end


%export xyz
if gdata.batch.out_xyz
    h = waitbar(0,'Creating .XYZ File, Please wait...');
    set(h,'name','Creating .XYZ File');
    
    if isempty(gdata.cdata)==1
        dlmwrite([pathname,filename,'.xyz'],[gdata.bdata.x,...
            gdata.bdata.y,gdata.bdata.zc],...
            'precision','%.3f','delimiter','\t')
    else
        dlmwrite([pathname,filename,'.xyz'],[gdata.cdata.xc,...
            gdata.cdata.yc,gdata.cdata.zc],...
            'precision','%.3f','delimiter','\t')
    end
    waitbar(1,h,'Done!');
    close(h);
end

if gdata.batch.out_kml
if isempty(gdata.cdata)==1
    gxyz=gdata.bdata;
    
    ind=find(isfinite(gxyz.zc) & ...
        isfinite(gxyz.x) & ...
        isfinite(gxyz.y) & ...
        isfinite(gxyz.lat) & ...
        isfinite(gxyz.lon) & ...
        isfinite(gxyz.mtime));
    

    zr=gxyz.zc(ind);
    lat=gxyz.lat(ind);
    lon=gxyz.lon(ind);
    dn=gxyz.mtime(ind);
else
    gxyz=gdata.cdata;
    ind=find(isfinite(gxyz.zc) & ...
        isfinite(gxyz.xc) & ...
        isfinite(gxyz.yc) & ...
        isfinite(gxyz.lat) & ...
        isfinite(gxyz.lon) & ...
        isfinite(gxyz.mtime));
    
    zr=gxyz.zc(ind);
    lat=gxyz.lat(ind);
    lon=gxyz.lon(ind);
    dn=gxyz.mtime(ind);
end


%export kml
gescatter([pathname,filename,'.kml'],lon(1:gdata.ge.thin:end),...
    lat(1:gdata.ge.thin:end),zr(1:gdata.ge.thin:end),...
    'time',dn(1:gdata.ge.thin:end),...
    'clims',[gdata.ge.cmin gdata.ge.cmax],'colormap',gdata.ge.cmap,...
    'scale',gdata.ge.scale);
end




% %export meta
% export5(hfig,[],[],[filename,'.meta']);

%export netcdf
if gdata.batch.out_nc
     h = waitbar(0,'Creating .NC File, Please wait...');
    export6(hfig,[],[pathname,filename,'.nc']);
    waitbar(1,h,'Done!')
    close(h);
end

if gdata.batch.out_shp
    
    if isempty(gdata.cdata)==1
        shpdata.X=gdata.bdata.x;
        shpdata.Y=gdata.bdata.y;
        shpdata.Z=gdata.bdata.zc;
    else
        shpdata.X=gdata.cdata.xc;
        shpdata.Y=gdata.cdata.yc;
        shpdata.Z=gdata.cdata.zc;
    end
    
    pidx=all(isfinite(cell2mat(struct2cell(shpdata)')),2);
    shpdata=structfun(@(x)(x(pidx)),shpdata,'un',0);
    
    shpdata.mtime=gdata.bdata.mtime(pidx);
    shpdata.Line_Number=gdata.hdr.lineNum;
    
    raw2shp([pathname,filename],shpdata)
end


gdata.outpath=pathname;
if gdata.outpath(end)~=filesep;
    gdata.outpath=[gdata.outpath,filesep];
end

guidata(hfig,gdata)

end

function export_shp(hfig,evnt) %#ok
gdata=guidata(hfig);

namer=strtok(gdata.bdata.filename,'.');

[filename, pathname] = uiputfile( ...
    {'*.shp', 'SHP Files'}, ...
    'Save as',[gdata.outpath,namer,'.xyz']);

if filename==0
    return
end

if isempty(gdata.cdata)==1
    shpdata.X=gdata.bdata.x;
    shpdata.Y=gdata.bdata.y;
    shpdata.Z=gdata.bdata.zc;
else
    shpdata.X=gdata.cdata.xc;
    shpdata.Y=gdata.cdata.yc;
    shpdata.Z=gdata.cdata.zc;
end

pidx=all(isfinite(cell2mat(struct2cell(shpdata)')),2);
shpdata=structfun(@(x)(x(pidx)),shpdata,'un',0);

shpdata.mtime=gdata.bdata.mtime(pidx);
shpdata.Line_Number=gdata.hdr.lineNum;

raw2shp([pathname,filename],shpdata)


gdata.outpath=pathname;
if gdata.outpath(end)~=filesep;
    gdata.outpath=[gdata.outpath,filesep];
end

guidata(hfig,gdata)

end

function raw2shp(shpfile,shpdata)

h = waitbar(0,'Creating .SHP File, Please wait...');
np=length(shpdata.X);
dc=cell(np,1);
shp=struct('Geometry','Point',...
    'Line',dc,...
    'Time',dc,...
    'X',dc,...
    'Y',dc,...
    'Z',dc);

lnum=num2cell(repmat(shpdata.Line_Number,np,1));
[shp(:).Line]=deal(lnum{:});


for i=1:np;
    shp(i).Time=datestr(shpdata.mtime(i),...
        'yyyy-mm-dd HH:MM:SS.FFF');
    shp(i).X=shpdata.X(i);
    shp(i).Y=shpdata.Y(i);
    shp(i).Z=shpdata.Z(i);
end

waitbar(0,h,'Writing .SHP File, Please wait...');
shapewrite(shp,shpfile)
waitbar(1,h,'Done!')
close(h);



end

%%%%%----------------------------------------------------------------------
function cleanProfile(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);

if ~isfield(gdata,'eraw');
    gdata.eraw=gdata.bdata.zraw;
    gdata.ezc=gdata.bdata.zc;
end
if gdata.numedits==0
    set(gdata.menu12,'visible','on');
    set(gdata.menu13,'visible','on');
end
gdata.numedits=gdata.numedits+1;

handles=findobj(gdata.ax,'type','line');
set(gdata.push5,'backgroundcolor','g');


if gdata.pan==1;
    set(gdata.toggle1,'value',0);
    pan off
end


if length(handles)==1
    [pl,xs,ys] = selectdata('sel','lasso'); %#ok
    gdata.bdata.zc(pl)=NaN;
    gdata.bdata.zraw(pl)=NaN;
    
    set(handles,'ydata',gdata.bdata.zc)
    
    ind=find(isfinite(gdata.bdata.zc));
    
    if gdata.alongflag==0
        gdata.xlimo=[min(gdata.bdata.distance(ind)),...
            max(gdata.bdata.distance(ind))];
        gdata.ylimo=[min(gdata.bdata.zc(ind)),...
            max(gdata.bdata.zc(ind))];
    else
        if ~isempty(gdata.topo)
            gdata.xlimo=[min([gdata.bdata.adist;gdata.topo.dist]),...
                max([gdata.bdata.adist;gdata.topo.dist])];
            gdata.ylimo=[min([gdata.bdata.zc;gdata.topo.z]),...
                max([gdata.bdata.zc;gdata.topo.z])];
        else
            gdata.xlimo=[min(gdata.bdata.adist(ind)),...
                max(gdata.bdata.adist(ind))];
            gdata.ylimo=[min(gdata.bdata.zc(ind)),...
                max(gdata.bdata.zc(ind))];
        end
    end
    
    waitfor(pl)
    
    gdata.edits{gdata.numedits}=pl;
    
    guidata(hfig,gdata);
    setFocus(hfig)
    
    
else
    
    if ~isempty(gdata.flag)
        gdata.xlims=get(gca,'xlim');
        gdata.ylims=get(gca,'ylim');
        hold off
        if gdata.alongflag==1
            gdata.l1=plot(gdata.bdata.adist,gdata.bdata.zc);
            hold on
            gdata.bb=plot(gdata.bdata.adist(gdata.flag==1),...
                gdata.bdata.zc(gdata.flag==1),'o',...
                'color',[0.6 0.6 0.6],'markersize',3,...
                'markerfacecolor',[0.6 0.6 0.6]);
            gdata.gg=plot(gdata.adistc,gdata.bdata.zc,'r-','linewidth',2);
        else
            
            gdata.l1=plot(gdata.bdata.distance,gdata.bdata.zc);
            hold on
            gdata.bb=plot(gdata.bdata.distance(gdata.flag==1),...
                gdata.bdata.zc(gdata.flag==1),'o',...
                'color',[0.6 0.6 0.6],'markersize',3,...
                'markerfacecolor',[0.6 0.6 0.6]);
            gdata.gg=plot(gdata.bdata.distance,gdata.cdata.zc,...
                'r-','linewidth',2);
        end
        
        
        
        set(gdata.ax,'xlim',gdata.xlims,'ylim',gdata.ylims);
        
        c1=legend('Raw Data','Outliers','Filtered Data');
        set(c1,'box','off','location','best')
        
        if ~isempty(gdata.topo);
            gdata.tlineh=plot(gdata.topo.dist,gdata.topo.z,'r-',...
                'linewidth',2);
            [pl,xs,ys] = selectdata('sel','lasso',...
                'ignore',[gdata.gg gdata.bb gdata.tlineh]); %#ok
        else
            
            [pl,xs,ys] = selectdata('sel','lasso',...
                'ignore',[gdata.gg gdata.bb]); %#ok
        end
        
        gdata.bdata.zc(pl)=NaN;
        gdata.bdata.zraw(pl)=NaN;
        
        set(gdata.l1,'ydata',gdata.bdata.zc)
        gdata.edits{gdata.numedits}=pl;
        
        guidata(hfig,gdata);
        
        waitfor(pl)
        xyzFilt(hfig);
    else
        gdata.xlims=get(gca,'xlim');
        gdata.ylims=get(gca,'ylim');
        
        hold off
        gdata.l1=plot(gdata.bdata.adist,gdata.bdata.zc);
        hold on
        gdata.tlineh=plot(gdata.topo.dist,gdata.topo.z,'r-',...
            'linewidth',2);
        set(gdata.ax,'xlim',gdata.xlims,'ylim',gdata.ylims);
        
        [pl,xs,ys] = selectdata('sel','lasso',...
            'ignore',[gdata.tlineh]); %#ok
        
        gdata.bdata.zc(pl)=NaN;
        gdata.bdata.zraw(pl)=NaN;
        
        gdata.edits{gdata.numedits}=pl;
        set(gdata.l1,'ydata',gdata.bdata.zc)
        guidata(hfig,gdata);
    end
    
end

if ishandle(gdata.offFig)
    figure(gdata.offFig);
    
    set(gdata.offPlot,'xdata',...
        gdata.bdata.distance(isfinite(gdata.bdata.zc)),...
        'ydata',gdata.bdata.offline(isfinite(gdata.bdata.zc)))
end

set(gdata.push5,'backgroundcolor',[0.9255    0.9137    0.8471]);

setFocus(hfig)

end

%%%%----------------------------------------------------------------------
function undo(hfig,evnt) %#ok
gdata=guidata(hfig);

if gdata.numedits>0
    gdata.bdata.zraw(gdata.edits{gdata.numedits})=...
        gdata.eraw(gdata.edits{gdata.numedits});
    
    
    if isfield(gdata,'sv');
        if gdata.invert==1
            zraw=-gdata.bdata.zraw;
        end
        zsos=apply_sos_prof(zraw,gdata.sv.sos_orig,...
            [gdata.sv.depth gdata.sv.sos],gdata.sv.use_mean_sos,...
            gdata.sv.mean_vel)-gdata.manoff;
        if gdata.invert==1
            zsos=-zsos;
        end
        
        gdata.bdata.zc=zsos-gdata.bdata.tide;
        gdata.bdata.zc(gdata.edits{gdata.numedits})=...
            zsos(gdata.edits{gdata.numedits})-...
            gdata.bdata.tide(gdata.edits{gdata.numedits});
        
    else
        gdata.bdata.zc(gdata.edits{gdata.numedits})=...
            (gdata.bdata.zraw(gdata.edits{gdata.numedits})+...
            gdata.manoff)-...
            gdata.bdata.tide(gdata.edits{gdata.numedits});
    end
    
    
    gdata.edits(gdata.numedits)=[];
    gdata.numedits=gdata.numedits-1;
    
    
    if gdata.rawflag==0
        if numel(gdata.bdata.zc)~=numel(get(gdata.l1,'xdata'));
            delete(gdata.l1)
            gdata.l1=plot(gdata.bdata.distance,...
                gdata.bdata.zc,'b');
        else
            set(gdata.l1,'ydata',gdata.bdata.zc);
        end
    else
        if numel(gdata.bdata.zraw)~=numel(get(gdata.l1,'xdata'));
            delete(gdata.l1)
            gdata.l1=plot((1:numel(gdata.bindata.vtime)),...
                gdata.bdata.zraw,'r');
        else
            set(gdata.l1,'ydata',gdata.bdata.zraw);
        end
    end
    
    if gdata.numedits==0
        set(gdata.menu12,'visible','off')
        set(gdata.menu13,'visible','off')
    end
    guidata(hfig,gdata);
end
end

%%%%----------------------------------------------------------------------
function showLineFile(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);

ind=find(isfinite(gdata.bdata.zc));

if isfield(gdata,'lnw')~=1;
    
    gdata.lnw=readLNW;
    guidata(hfig,gdata);
end

minx=min(arrayfun(@(x)(min(x.x)),gdata.lnw));
maxx=max(arrayfun(@(x)(max(x.x)),gdata.lnw));
miny=min(arrayfun(@(x)(min(x.y)),gdata.lnw));
maxy=max(arrayfun(@(x)(max(x.y)),gdata.lnw));

xr=maxx-minx;
yr=maxy-miny;



h3=figure;
set(h3,'Name','Line File Viewer','numbertitle','off');

if xr>yr
    subplot(211)
    xs=0.025;
    ys=0.05;
else
    subplot(121)
    xs=0.05;
    ys=0.025;
end
hold on

%make ~50 labels per plot
numlines=length(gdata.lnw);
int=ceil(numlines/50);
for i=1:numlines;
    plot(gdata.lnw(i).x,gdata.lnw(i).y);
    if rem(i,int)==0
        text(gdata.lnw(i).x(1),gdata.lnw(i).y(1),...
            num2str(gdata.lnw(i).name),'fontsize',10)
    end
end
axis equal

if isempty(gdata.cdata)==1
    plot(gdata.bdata.x(ind),gdata.bdata.y(ind),'k.')
    plot(gdata.bdata.x(min(ind)),gdata.bdata.y(min(ind)),'rs')
    
    minx2=min(gdata.bdata.x(ind))-xs*xr;
    maxx2=max(gdata.bdata.x(ind))+xs*xr;
    miny2=min(gdata.bdata.y(ind))-ys*yr;
    maxy2=max(gdata.bdata.y(ind))+ys*yr;
    
else
    plot(gdata.cdata.xc(ind),gdata.cdata.yc(ind),'k.')
    plot(gdata.cdata.xc(min(ind)),gdata.cdata.yc(min(ind)),'rs')
    
    minx2=min(gdata.cdata.xc(ind))-xs*xr;
    maxx2=max(gdata.cdata.xc(ind))+xs*xr;
    miny2=min(gdata.cdata.yc(ind))-ys*yr;
    maxy2=max(gdata.cdata.yc(ind))+ys*yr;
    
end

p1(1)=line([minx2 maxx2],[miny2 miny2]);
p1(2)=line([minx2 maxx2],[maxy2 maxy2]);
p1(3)=line([minx2 minx2],[miny2 maxy2]);
p1(4)=line([maxx2 maxx2],[miny2 maxy2]);
set(p1(1:end),'color','r','linewidth',2)

if xr>yr
    subplot(212)
else
    subplot(122)
end
hold on

%make ~50 labels per plot
numlines=length(gdata.lnw);
for i=1:numlines;
    plot(gdata.lnw(i).x,gdata.lnw(i).y);
end
axis equal

if isempty(gdata.cdata)==1
    plot(gdata.bdata.x(ind),gdata.bdata.y(ind),'k.')
    plot(gdata.bdata.x(min(ind)),gdata.bdata.y(min(ind)),'rs')
    
else
    plot(gdata.cdata.xc(ind),gdata.cdata.yc(ind),'k.')
    plot(gdata.cdata.xc(min(ind)),gdata.cdata.yc(min(ind)),'rs')
    
end

set(gca,'xlim',[minx2 maxx2],'ylim',[miny2 maxy2]);
title(['Line ',sprintf('%d',gdata.hdr.lineNum)],'fontsize',14)


end


%%%%----------------------------------------------------------------------
function showLineInfo(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);
h2=figure;
set(h2,'Name','Line Info','menubar','none','numbertitle','off');
data{1,1}='filename:';
data{1,2}=gdata.bdata.filename;
data{2,1}='linefile:';
if isfield(gdata.hdr,'linefile')
    data{2,2}=gdata.hdr.linefile;
end
data{3,1}='line number:';
data{3,2}= num2str(gdata.hdr.lineNum);
data{4,1}='Date Collected:';
data{4,2}=datestr(gdata.hdr.mtime);
data{5,1}='Line length (m):';
data{5,2}=sprintf('%.2f',gdata.bdata.distance(end));
data{6,1}='Elapsed Time (min)';
data{6,2}=sprintf('%.2f',(gdata.bdata.mtime(end)-...
    gdata.bdata.mtime(1))*1440);
data{7,1}='Mean Speed (m/s)';
data{7,2}=sprintf('%.2f',gdata.bdata.distance(end)/...
    ((gdata.bdata.mtime(end)-gdata.bdata.mtime(1))*86400));

%calculate mean number of pings per m
pingrate=1/(mean(diff(gdata.bdata.mtime))*86400);
data{8,1}='Mean Pings per m';
data{8,2}=sprintf('%.2f',pingrate/str2double(data{7,2}));

if isfield(gdata.bdata,'tide')
    a=ver;
    tboxes={a.Name}';
    if any(strcmpi('signal processing toolbox',tboxes));
        ind=find(isfinite(gdata.bdata.tide));
        
        fs=1/(mean(diff(gdata.bdata.mtime(ind)))*86400);
        nfft=512;
        window = hanning(nfft);
        lf_cutoff=1/100;
        hf_cutoff=1/3;
        noverlap=nfft/2;
        
        
        amp=detrend(gdata.bdata.tide(ind));
        
        [Gpp, f] = p_welch(amp,nfft, fs, window, noverlap );
        df = f(3)-f(2);
        
        ff = find(f>=lf_cutoff,1,'first');
        lf = find(f<=hf_cutoff,1,'last');
        
        Hs = 4.004*sqrt(sum( Gpp(ff:lf)*df ) );
        
        [junk,ind2]=max(Gpp(ff:lf));%#ok
        f2=f(ff:lf);
        tp=1/f2(ind2);
        
        data{9,1}='Significant Wave Height (m)';
        data{9,2}=sprintf('%.2f',Hs);
        data{10,1}='Wave Period (s, relative)';
        data{10,2}=sprintf('%.2f',tp);
        
        gdata.wvs.fname=gdata.bdata.filename;
        gdata.wvs.mtime=gdata.bdata.mtime(ind);
        gdata.wvs.amp=amp;
        gdata.wvs.freq=f(ff:lf);
        gdata.wvs.famp=Gpp(ff:lf);
        gdata.wvs.hs=Hs;
        gdata.wvs.tp=tp;
        
        guidata(hfig,gdata);
        
    end
end

if isempty(gdata.tcorr)~=1
    len1=length(data);
    data{len1+1,1}='Tide correction (m)';
    data{len1+1,2}=sprintf('%0.2f',gdata.tcorr);
end

if isfield(gdata,'sv')
    len1=length(data);
    data{len1+1,1}='Original Speed of Sound (m/s)';
    if isfield(gdata.sv,'sos_orig')
        if ischar(gdata.sv.sos_orig)
            
            data{len1+1,2}=gdata.sv.sos_orig;
            
        else
            data{len1+1,2}=sprintf('%0.2f',gdata.sv.sos_orig);
        end
    else
        data{len1+1,2}='unknown';
    end
    
    data{len1+2,1}='Speed of Sound Applied (m/s)';
    if gdata.sv.use_mean_sos
        if ischar(gdata.sv.mean_vel)
            data{len1+2,2}=gdata.sv.mean_vel;
        else
            data{len1+2,2}=sprintf('%0.2f',gdata.sv.mean_vel);
        end
    else
        data{len1+2,2}='Profile';
    end
end


mtable=uitable(h2,'Data',data,'columnname',{'Property','Value'});
set(mtable,'units','normalized','position',[0.1 0.25 0.8 0.75]);

pix=get(h2,'position');
cwidth=floor([0.3 0.45]*pix(3));
set(mtable,'columnwidth',num2cell(cwidth))


end

%%%%-----------------------------------------------------------------------
function showLog(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);

h2=figure;
set(h2,'Name',['Log File Contents: ' ,gdata.logname],...
    'menubar','none','numbertitle','off');



mtable=uitable(h2,'data',[gdata.fnames,gdata.filesize],...
    'columnname',{'Filenames','File Size (kBytes)'});
set(mtable,'units','normalized','position',[0.2 0.1 0.6 0.8])
pix=get(h2,'position');
cwidth=floor(0.273*pix(3));
set(mtable,'columnwidth',{cwidth})

end
%%%%-----------------------------------------------------------------------
function calcHs(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);

ind=find(isfinite(gdata.bdata.tide));

fs=1/(mean(diff(gdata.bdata.mtime(ind)))*86400);
nfft=512;
window = hanning(nfft);
lf_cutoff=1/100;
hf_cutoff=1/3;
noverlap=nfft/2;


amp=detrend(gdata.bdata.tide(ind));

[Gpp, f] = p_welch(amp,nfft, fs, window, noverlap );
df = f(3)-f(2);

ff = find(f>=lf_cutoff,1,'first');
lf = find(f<=hf_cutoff,1,'last');

Hs = 4.004*sqrt(sum( Gpp(ff:lf)*df ) );

h3=figure;
set(h3,'Name','Wave height Estimate','menubar','none','numbertitle','off');
subplot(1,3,1:2);
plot(gdata.bdata.mtime(ind),amp);
set(gca,'xlim',[min(gdata.bdata.mtime(ind)) ,...
    max(gdata.bdata.mtime(ind))]);
pos=get(gca,'position');
pos(2)=pos(2)+0.2;
pos(4)=pos(4)-0.4;
set(gca,'position',pos)

xlims=get(gca,'xlim');
c1=line([xlims(1) xlims(2)],[0 0]);
set(c1,'color','r')

datetick('x',13,'keeplimits','keepticks')
xlabel('Time (HH:MM:SS)','fontsize',12)
ylabel('Relative antenna height (m)',...
    'fontsize',12)

subplot(1,3,3)
plot(f(ff:lf),Gpp(ff:lf));
xl=get(gca,'xlim');
yl=get(gca,'ylim');

xr=diff(xl);
yr=diff(yl);

text(xl(1)+0.5*xr,yl(1)+0.9*yr,...
    ['H_s = ',num2str(round(Hs*100)/100),' m']);

[junk,ind2]=max(Gpp(ff:lf)); %#ok
f2=f(ff:lf);
tp=1/f2(ind2);
text(xl(1)+0.5*xr,yl(1)+0.7*yr,...
    ['Tp = ',num2str(round(tp*10)/10),' s']);


pos=get(gca,'position');
pos(2)=pos(2)+0.2;
pos(4)=pos(4)-0.4;
set(gca,'position',pos)

xlabel('Frequency (Hz)','fontsize',12);
ylabel('Spectral Density (m^2/Hz)',...
    'fontsize',12)

gdata.wvs.fname=gdata.bdata.filename;
gdata.wvs.mtime=gdata.bdata.mtime(ind);
gdata.wvs.amp=amp;
gdata.wvs.freq=f(ff:lf);
gdata.wvs.famp=Gpp(ff:lf);
gdata.wvs.hs=Hs;
gdata.wvs.tp=tp;


guidata(hfig,gdata);

end
%%%%-----------------------------------------------------------------------
function dispHelp(hfig,eventdata,handles) %#ok
gdata=guidata(hfig);

str0=strjust(char(['Version',sprintf(' %0.2f',gdata.tv_ver)],...
    ' written by Andrew Stevens',...
    ['Last Modified ',gdata.modified],' '),'left');
str1=strjust(char('If you have questions',...
    'or want to report a bug, contact me:',' '),'left');
str2=strjust(char('Andrew Stevens',...
    'astevens@usgs.gov', 'tel: 831 460 7424'),'left');

str={str0;str1;str2};
h=helpdlg(str,'About CPS Processing Module');
g=get(h,'children');
set(g(1:2),'fontsize',18)

end
%%%%%---------------------------------------------------------------------
function dispKeys(hFig,eventdata,handles) %#ok

str={'The following is a list of keyboard shortcuts:';...
    ' ';...
    'e:  edit manually';...
    'u:  undo last edit';...
    'z:  zoom box';...
    '-:  zoom out';...
    '=:  zoom in';...
    'f:  show full extents';...
    'l:  local filter';...
    'o:  show offline plot';...
    'v:  toggle between raw (if available) and digitized data';...
    'g:  show gps plot';...
    'w:  calculate significant wave height';...
    'd:  hand digitize (raw data mode only)';...
    'm:  measure distance tool'};

h=helpdlg(str,'Keyboard Shortcuts');
g=get(h,'children');
set(g(1:2),'fontsize',18)


end

%%%%----------------------------------------------------------------------
function manualOffset(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);
prompt={'Manual offset (m):'};
name='Manual Offset';
numlines=1;
defaultanswer={num2str(gdata.manoff)};

answer=inputdlg(prompt,name,numlines,defaultanswer);
if isempty(answer)
    return
end
gdata.manoff=str2double(answer{:});

%apply manual offset

if isfield(gdata,'sv');
    if gdata.invert==1
        zraw=-gdata.bdata.zraw;
    end
    zsos=apply_sos_prof(zraw,gdata.sv.sos_orig,...
        [gdata.sv.depth gdata.sv.sos],gdata.sv.use_mean_sos,...
        gdata.sv.mean_vel)-gdata.manoff;
    if gdata.invert==1
        zsos=-zsos;
    end
    
    gdata.bdata.zc=zsos-gdata.bdata.tide;
    
else
    gdata.bdata.zc=(gdata.bdata.zraw+gdata.manoff)-...
        gdata.bdata.tide;
end

gdata.ylimo=gdata.ylimo-gdata.manoff;
gdata.ylimo=[min(gdata.bdata.zc),...
    max(gdata.bdata.zc)];

gdata.info.manual_offset=gdata.manoff;

if isempty(gdata.cdata)~=1
    guidata(hfig,gdata);
    xyzFilt(hfig)
else
    %update plot
    if gdata.rawflag~=1;
        set(gdata.l1,'ydata',gdata.bdata.zc)
        set(gdata.ax,'ylim',gdata.ylimo)
    end
    guidata(hfig,gdata);
    
end


end
%%%%----------------------------------------------------------------------
function interpDist(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);
prompt={'Maximum Distance to Interpolate (m):'};
name='Interp. Distance';
numlines=1;
defaultanswer={num2str(gdata.maxGap)};

answer=inputdlg(prompt,name,numlines,defaultanswer);
gdata.maxGap=str2double(answer{:});

guidata(hfig,gdata);

end
%%%%----------------------------------------------------------------------
function zlims(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);
prompt={'Min. Raw Sounding (m):';...
    'Max. Raw Sounding (m)'};

name='Sounder Limits';
numlines=1;
defaultanswer={num2str(gdata.minz);...
    num2str(gdata.maxz)};

answer=inputdlg(prompt,name,numlines,defaultanswer);
if ~isempty(answer)
    gdata.minz=str2double(answer{1});
    gdata.maxz=str2double(answer{2});
    guidata(hfig,gdata);
    
    movePhoto(hfig)
end
end

%%%%----------------------------------------------------------------------
function setmaxraw(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);
prompt={'Max. Raw Sounding (m):'};

name='Raw Data Limits';
numlines=1;
oldlimit=gdata.maxrawz;
defaultanswer={num2str(gdata.maxrawz)};

answer=inputdlg(prompt,name,numlines,defaultanswer);
if ~isempty(answer)
    gdata.maxrawz=str2double(answer{1});
    
    if gdata.maxrawz<oldlimit
        gdata.bindata=[];
    end
    
    guidata(hfig,gdata);
end
end


%%%%----------------------------------------------------------------------
function setrclims(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);
clims=get(gdata.ax1,'clim');

prompt={'Min.:';'Max.'};
name='Raw Data Color Limits:';

numlines=1;
defaultanswer=cellfun(@(x)(sprintf('%0.0f',x)),num2cell(clims),'un',0);

answer=inputdlg(prompt,name,numlines,defaultanswer');

if ~isempty(answer)
    gdata.rclims=sort(abs(str2double(answer)));
    set(gdata.ax1,'clim',gdata.rclims);
    guidata(hfig,gdata);
end
end


%%%%-----------------------------------------------------------------------
function run_sos_gui(hfig,evnt) %#ok

gdata=guidata(hfig);

if isfield(gdata,'sv')
    gdata.sv=sos_gui(gdata.sv);
else
    gdata.sv=sos_gui;
end

if ~isempty(gdata.sv)
    if gdata.invert==1
        zraw=-gdata.bdata.zraw;
    end
    
    zsos=apply_sos_prof(zraw,gdata.sv.sos_orig,...
        [gdata.sv.depth gdata.sv.sos],gdata.sv.use_mean_sos,...
        gdata.sv.mean_vel)-gdata.manoff;
    
    
    if gdata.invert==1
        zsos=-zsos;
    end
    
    gdata.bdata.zc=zsos-gdata.bdata.tide;
    gdata.info.orig_speed_of_sound=sprintf('%0.1f',...
        gdata.sv.sos_orig);
    if gdata.sv.use_mean_sos
        gdata.info.applied_speed_of_sound=sprintf('%0.1f',...
            gdata.sv.mean_vel);
    else
        gdata.info.applied_speed_of_sound='Profile';
    end
    
    gdata.ylimo=[min(gdata.bdata.zc),...
        max(gdata.bdata.zc)];
    
    
    if ~isempty(gdata.cdata)
        guidata(hfig,gdata);
        xyzFilt(hfig);
    else
        if gdata.rawflag==0
            set(gdata.l1,'ydata',gdata.bdata.zc);
            guidata(hfig,gdata);
        end
    end
    guidata(hfig,gdata);

end


end

function zsos=apply_sos_prof(zraw,osv,sosp,usemean,meanvel)
%zraw - raw depths
%osv - original sound velocity
%sosp - m x 2 matrix with col 1=depths, col2 = sos
%usemean - flag 0 = use profile, 1 = use single mean val
%meanvel - mean sos (m/s)

if ~usemean
    T=((2.*zraw)./osv)./2;
    ptime= diff([0;sosp(:,1)])./sosp(:,2);
    
    tedges=[0;cumsum(ptime);inf];
    dedges=[0;sosp(:,1)];
    sedges=[sosp(:,2);sosp(end,2)];
    [~,bin]=histc(T,tedges);
    zsos=nan(numel(zraw),1);
    zsos(bin~=0)=dedges(bin(bin~=0))+...
        ((T(bin~=0)-tedges(bin(bin~=0))).*sedges(bin(bin~=0)));
else
    ratio=meanvel/osv;
    zsos= zraw.*ratio;
end


end


function sv = sos_gui(varargin)
% The basic layout of this GUI was made with the help of guidegetter,
% available on the File Exchange at Mathworks.com


if nargin>0
    sv=varargin{1};
    
    if ~isempty(sv.depth)
        
        data=cell(length(sv.depth),2);
        
        data(:,1)=cellfun(@(x)(sprintf('%0.2f',x)),...
            num2cell(sv.depth),'un',0);
        data(:,2)=cellfun(@(x)(sprintf('%0.2f',x)),...
            num2cell(sv.sos),'un',0);
        list_enable='on';
        gd.sv=sv;
        
        if ~isfield(gd.sv,'sos_orig')
            gd.sv.sos_orig=NaN;
        end
        
    else
        if isempty(sv) || ~isfield(sv,'mean_vel')
            gd.sv.depth=[];
            gd.sv.sos=[];
            gd.sv.time=[];
            gd.sv.sos_orig=1500;
            gd.sv.use_mean_sos=0;
            gd.sv.use_prof=1;
            gd.sv.mean_vel=[];
            
            list_enable='off';
            data=cell(20,2);
        else
            gd.sv.depth=[];
            gd.sv.sos=[];
            gd.sv.time=[];
            if isfield(sv,'sos_orig')
                gd.sv.sos_orig=sv.sos_orig;
            else
                gd.sv.sos_orig=NaN;
            end
            gd.sv.use_mean_sos=1;
            gd.sv.use_prof=0;
            gd.sv.mean_vel=sv.mean_vel;
            
            list_enable='off';
            data=cell(20,2);
        end
        
    end
else
    gd.sv.depth=[];
    gd.sv.sos=[];
    gd.sv.time=[];
    gd.sv.sos_orig=1500;
    gd.sv.use_mean_sos=0;
    gd.sv.use_prof=1;
    gd.sv.mean_vel=[];
    
    list_enable='off';
    data=cell(20,2);
    
    
    
end



hf = figure('units','normalized',...
    'position',[0.354 0.491 0.293 0.414],...
    'menubar','none',...
    'name','sos_gui',...
    'numbertitle','off',...
    'color',[0.94 0.94 0.94]);

uicontrol(hf,'style','text',...
    'units','normalized',...
    'position',[0.05 0.895 0.18 0.06],...
    'string','Echosounder Sound Vel. (m/s)',...
    'backgroundcolor',[0.94 0.94 0.94]);
gd.edit0 = uicontrol(hf,'style','edit',...
    'units','normalized',...
    'position',[0.25 0.905 0.139 0.0492],...
    'string',sprintf('%0.2f',gd.sv.sos_orig),...
    'backgroundcolor',[1 1 1]);

gd.table1 = uitable(hf,'units','normalized',...
    'position',[0.05 0.215 0.35 0.673],...
    'ColumnName',{'Depth (m)','SV (m/s)'},...
    'data',data,...
    'ColumnEditable', [true true],...
    'CellEditCallback',@sos_edit_table,...
    'enable',list_enable);
gd.ax1 = axes('position',[0.55 0.215 0.4 0.665],...
    'ydir','rev',...
    'xaxislocation','top',...
    'nextplot','add');
gd.lh=plot(gd.sv.sos,gd.sv.depth,'ko-',...
    'markerfacecolor','k');

xlabel('\bf\itSound Velocity (m/s)')
ylabel('\bf\itDepth (m)')

uicontrol(hf,'style','text',...
    'units','normalized',...
    'position',[0.05 0.13 0.18 0.0425],...
    'string','Mean Velocity (m/s)',...
    'backgroundcolor',[0.94 0.94 0.94]);
gd.edit1 = uicontrol(hf,'style','edit',...
    'units','normalized',...
    'position',[0.25 0.125 0.139 0.0492],...
    'string',sprintf('%0.2f',gd.sv.mean_vel),...
    'backgroundcolor',[1 1 1],...
    'callback',@check_is_ready);

bg = uibuttongroup('units','normalized',...
    'position',[0.05 0.0125 0.35 0.1]);
gd.sv_prof = uicontrol(bg,'Style',...
    'radiobutton',...
    'units','normalized',...
    'String','Use Profile',...
    'Position',[0.1 0.5 0.7 0.35],...
    'HandleVisibility','off',...
    'callback',@check_is_ready);

gd.sv_mean = uicontrol(bg,'Style','radiobutton',...
    'String','Use Mean SV',...
    'units','normalized',...
    'Position',[0.1 0.1 0.7 0.35],...
    'HandleVisibility','off',...
    'value',gd.sv.use_mean_sos,...
    'callback',@check_is_ready);


uicontrol(hf,'style','pushbutton',...
    'units','normalized',...
    'position',[0.553 0.0134 0.157 0.0895],...
    'string','cancel',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'callback',@(h,e)(close(hf)));
gd.push2 = uicontrol(hf,'style','pushbutton',...
    'units','normalized',...
    'position',[0.733 0.0134 0.157 0.0895],...
    'string','Apply',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'enable',list_enable,...
    'callback',@(h,e)(uiresume(hf)));

menu1=uimenu('label','File');
uimenu(menu1,'label','Open',...
    'callback',@sos_gui_open)

menu2=uimenu('label','Utilities');
uimenu(menu2,'label','Avg. Profiles',...
    'callback',@sos_avg_prof);

gd.isready=[0 0];

guidata(hf,gd);

uiwait

if ishandle(hf) %if user pressed cancel
    gd=guidata(hf);
    
    sv=gd.sv;
    sv.mean_vel=str2double(get(gd.edit1,'string'));
    sv.sos_orig=str2double(get(gd.edit0,'string'));
    sv.use_mean_sos=get(gd.sv_mean,'value');
    sv.use_prof=get(gd.sv_prof,'value');
    close(hf)
else
    sv=[];
    return
    
end

end

function check_is_ready(hf,evnt) %#ok

gd=guidata(hf);
val=get(gd.sv_mean,'value');
if val==1
    gd.isready(1)=1;
else 
    if isempty(gd.sv.sos)
        gd.isready(1)=0;
    else 
        gd.isready(1)=1;
    end
end

sosnew=str2double(get(gd.edit1,'string'));
if isfinite(sosnew)
    gd.isready(2)=1;
end

if all(gd.isready==1)
    set(gd.push2,'enable','on')
else 
    set(gd.push2,'enable','off')
end

end


function sos_gui_open(hf,evnt) %#ok

gd=guidata(hf);

[filename, pathname, fidx] = uigetfile( ...
    {'*.mat', 'YSI CastAway Files (*.mat)';...
    '*.vel', 'Hypack SV files (*.vel)'},...
    'Select a file');

if filename==0
    return
end

switch fidx
    case 1
        sosp=load([pathname,filename]);
        sv.depth=sosp.Depth+...
            (diff([0;sosp.Depth])/2); %depth should be end of bin,...
        %not bin center
        sv.sos=sosp.Sound_velocity;
        sv.time=diff([0;sv.depth])./sv.sos;
        
    case 2
        
        fmt='%f %f';
        fid=fopen([pathname,filename],'r');
        svdata=textscan(fid,fmt,'headerlines',1);
        sv.depth=svdata{1};
        sv.sos=svdata{2};
        sv.time=diff([0;sv.depth])./sv.sos;
        fclose(fid);
end

set(gd.table1,'enable','on');
set(gd.push2,'enable','on');

gd.sv=sv;
guidata(hf,gd);
update_sos_gui(hf);

end
function update_sos_gui(hf,evnt) %#ok

gd=guidata(hf);

data=cell(length(gd.sv.depth),2);
data(:,1)=cellfun(@(x)(sprintf('%0.2f',x)),...
    num2cell(gd.sv.depth),'un',0);
data(:,2)=cellfun(@(x)(sprintf('%0.2f',x)),...
    num2cell(gd.sv.sos),'un',0);

set(gd.table1,'data',data)

set(gd.edit1,'string',sprintf('%0.2f',nanmean(gd.sv.sos)));

if ~isempty(gd.lh)
    set(gd.lh,'xdata',gd.sv.sos,'ydata',gd.sv.depth);
else
    gd.lh=plot(gd.sv.sos,gd.sv.depth,'ko-',...
        'markerfacecolor','k',...
        'parent',gd.ax1);
end
guidata(hf,gd);
end

function sos_avg_open(hf,evnt) %#ok

gd=guidata(hf);

[filename, pathname, fidx] = uigetfile( ...
    {'*.mat', 'YSI CastAway Files (*.mat)';...
    '*.vel', 'Hypack SV files (*.vel)'},...
    'Select a file','multiselect','on');

if pathname==0
    return
end

if ~iscell(filename)
    filename={filename};
end

switch fidx
    case 1
        sos=cellfun(@(x)(load([pathname,x])),filename);
        gd.sos_avg_depth=arrayfun(@(x)(x.Depth+...
           diff([0;x.Depth])./2),sos,'un',0);
        gd.sos_avg_p=arrayfun(@(x)(x.Sound_velocity),sos,'un',0);
        
        
    case 2
        
        fmt='%f %f';
        gd.sos_avg_depth=cell(1,length(filename));
        gd.sos_avg_p=cell(1,length(filename));
        for i =1:length(filename)
            fid=fopen([pathname,filename{i}],'r');
            svdata=textscan(fid,fmt,'headerlines',1);
            fclose(fid);
            gd.sos_avg_depth{i}=svdata{1};
            gd.sos_avg_p{i}=svdata{2};
        end
        
end

gd.sos_avg_ready=1;
gd.sos_path=pathname;
set(gd.list,'string',filename)

if all(gd.sos_avg_ready==1);
    set(gd.done,'enable','on')
end

guidata(hf,gd);
end


function [] = sos_avg_prof(hf,evnt) %#ok
% The basic layout of this GUI was made with the help of guidegetter,
% available on the File Exchange at Mathworks.com

gd=guidata(hf);

hf2 = figure('units','normalized',...
    'position',[0.354 0.651 0.177 0.254],...
    'menubar','none',...
    'name','sos_avg_prof',...
    'numbertitle','off',...
    'color',[0.94 0.94 0.94]);

uicontrol(hf2,'style',...
    'pushbutton','units',...
    'normalized','position',...
    [0.0914 0.752 0.204 0.128],...
    'string','Open Files',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'callback',@sos_avg_open);
gd2.list = uicontrol(hf2,'style','listbox',...
    'units','normalized',...
    'position',[0.327 0.628 0.558 0.248],...
    'string','No files selected',...
    'backgroundcolor',[1 1 1]);


uipanel1 = uipanel('parent',hf2,...
    'units','normalized',...
    'position',[0.0885 0.0803 0.537 0.314],...
    'title','');
gd2.binsize = uicontrol(uipanel1,'style','edit',...
    'units','normalized',...
    'position',[0.0562 0.22 0.433 0.432],...
    'string','0.5',...
    'backgroundcolor',[1 1 1]);
gd2.btxt = uicontrol(uipanel1,'style','text',...
    'units','normalized',...
    'position',[0.551 0.259 0.399 0.351],...
    'string','Bin size (m)',...
    'backgroundcolor',[0.94 0.94 0.94]);

gd2.done=uicontrol(hf2,'style','pushbutton',...
    'units','normalized',...
    'position',[0.678 0.0839 0.204 0.109],...
    'string','Done',...
    'enable','off',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'callback',@sos_do_avg);

uicontrol(hf2,'style','pushbutton',...
    'units','normalized',...
    'position',[0.678 0.23 0.204 0.109],...
    'string','Cancel',...
    'backgroundcolor',[0.94 0.94 0.94],...
    'callback',@(h,e)(close(hf2)));

gd2.sos_path=[pwd,filesep];
gd2.sos_avg_ready=0;
guidata(hf2,gd2)

uiwait
if ishandle(hf2) %if user pressed cancel
    gd2=guidata(hf2);
    close(hf2)
    
    set(gd.table1,'enable','on');
    set(gd.push2,'enable','on');
    
    sv.depth=gd2.dout(:,1);
    sv.sos=gd2.dout(:,2);
    sv.time=diff([0;sv.depth])./sv.sos;
    gd.sv=sv;
    
    guidata(hf,gd);
end
update_sos_gui(hf);
end

function sos_do_avg(hf,evnt) %#ok

gd=guidata(hf);

binsize=str2double(get(gd.binsize,'string'));

dall=cat(1,gd.sos_avg_depth{:});
zbins=[min(dall):binsize:max(dall),max(dall)];

[n,bin]=histc(dall,zbins);
bin(bin==numel(n))=numel(n)-1;
zcen=zbins(2:end)-(binsize/2);

sos_avg=accumarray(bin,cat(1,gd.sos_avg_p{:}),...
    [numel(n)-1 1],@mean,NaN);

%make sure no nans
dout=[zcen' sos_avg];
dout(any(isnan(dout),2),:)=[];

gd.dout=dout;
guidata(hf,gd);

uiresume
end

function sos_edit_table(hf,evnt)

gd=guidata(hf);
r = evnt.Indices(1);
c = evnt.Indices(2);
data=get(gd.table1,'Data');
edata=str2double(evnt.EditData);
if isnan(edata)
    data(r,:)=[];
else
    
    data{r,c}=sprintf('%0.2f',edata);
end

mdata=str2double(data);
gd.sv.depth=mdata(:,1);
gd.sv.sos=mdata(:,2);

guidata(hf,gd);
update_sos_gui(hf);

end


%%%%%----------------------------------------------------------------------
function [dstr,hdr] = read_meta(filen)
% READ_META - reads contents of a .meta file

if ~exist(filen,'file')
    error('File not found.')
end

fid=fopen(filen,'r');

%read the meta data
hdr_strings={'Meta File Version:';...
    'Raw Filename:';...
    'Vessel ID:';...
    'Vessel Operator:';...
    'Echosounder Type:';...
    'Transducer Frequency (kHz):';...
    'Transducer Blanking (m):';...
    'Transducer Gain:';...
    'Transducer Transmit Power:';...
    'Speed of Sound (m/s):';...
    'Post-Process GPS Filename:';...
    'GPS Antenna Offset (m):';...
    'Manual Elevation Offset (m):';...
    'Base Station ID: ';...
    'Data Analyst:';...
    'Transect Viewer Version:'};
hdr_par={'meta_ver';...
    'raw_filename';...
    'vessel_id';...
    'vessel_operator';...
    'echosounder_type';...
    'xducer_frequency';...
    'xducer_blanking';...
    'xducer_gain';...
    'xducer_tx';...
    'speed_of_sound';...
    'ppk_gps_filename';...
    'antenna_offset';...
    'manual_offset';...
    'base_stn_id';...
    'data_analyst';...
    'tv_ver'};

nline=1;
hline=1;
fmt='%[^\t]%[^\n]';
while ~isempty(hline);
    hline=fgetl(fid);
    if isempty(hline) || all(isspace(hline))
        break
    end
    hpar=textscan(hline,fmt);
    
    [~,ia] = intersect(hdr_strings,hpar{1});
    hdr.(hdr_par{ia})=hpar{2}{:};
    nline=nline+1;
end


%read the main data block
dhdrs=textscan(fgetl(fid),'%s',...
    'delimiter','\t');
data=textscan(fid,'%[^\n]');

time_fmt='%f-%f-%f %f:%f:%f';

%fix bug in past files where a tab delimtier was missing
if str2double(hdr.meta_ver)<1.02
    dhdrs{1}(1)={'Date and Time (yyyy-mm-dd HH:MM:SS)'};
    bad_hdr=strmatch('Offline (m)Distance (m)',dhdrs{1}); %#ok
    if ~isempty(bad_hdr);
        dhdrs{1}(bad_hdr)={'Offline (m)'};
        dhdrs{1}(bad_hdr+1)={'Distance (m)'};
    end
end
full_fmt=[time_fmt,repmat('\t%f',1,length(dhdrs{1})-1)];

dmat=num2cell(cell2mat(cellfun(@(x)(textscan(x,full_fmt,...
    'collectoutput',1)),data{1})),1);

%reformat data into struct
dstr_fields1={'Year';'Month';'Day';...
    'Hour';'Minute';'Second'};
dstr_fields2=strtok(dhdrs{1}(2:end),' ');
dfields=cat(1,dstr_fields1,dstr_fields2);

dstr=cell2struct(dmat',dfields);



end

%%%%%----------------------------------------------------------------------

function xducer=read_echart_settings(varargin)
% read data from echart presets file (.xml)

error(nargchk(0,1,nargin,'struct'));
if nargin<1;
    [filename, pathname] = uigetfile( ...
        {'*.xml', 'XML Files (*.xml)'},...
        'Select a .xml file');
    fname=fullfile(pathname,filename);
else
    fname=varargin{1};
    if ~exist(fname,'file')
        error('File not found.')
    end
end
xml=xml2struct(fname);

nparams=length(xml.ArrayOfParameter.Parameter);
parameters=cellfun(@(x)(xml.ArrayOfParameter.Parameter{x}.name.Text),...
    num2cell(1:nparams),'un',0);
values=cellfun(@(x)(xml.ArrayOfParameter.Parameter{x}.value.Text),...
    num2cell(1:nparams),'un',0);


[~,pidx]=unique(parameters);
params=parameters(sort(pidx));
p_space=cellfun(@(x)(isspace(x)),params,'un',0);

for i =1:length(p_space);
    params{i}(p_space{i})='_';
    params{i}(params{i}=='.')=[];
end


vals=num2cell(str2double(values(sort(pidx))));
xducer=cell2struct(vals',params');

end

%--------------------------------------------------------------------------
function [ s ] = xml2struct( file )
%Convert xml file into a MATLAB structure
% [ s ] = xml2struct( file )
%
% A file containing:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</Element>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Will produce:
% s.XMLname.Attributes.attrib1 = "Some value";
% s.XMLname.Element.Text = "Some text";
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
% s.XMLname.DifferentElement{1}.Text = "Some more text";
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
% s.XMLname.DifferentElement{2}.Text = "Even more text";
%
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012

if (nargin < 1)
    clc;
    help xml2struct
    return
end

if isa(file, 'org.apache.xerces.dom.DeferredDocumentImpl') || isa(file, 'org.apache.xerces.dom.DeferredElementImpl')
    % input is a java xml object
    xDoc = file;
else
    %check for existance
    if (exist(file,'file') == 0)
        %Perhaps the xml extension was omitted from the file name. Add the
        %extension and try again.
        if (isempty(strfind(file,'.xml')))
            file = [file '.xml'];
        end
        
        if (exist(file,'file') == 0)
            error(['The file ' file ' could not be found']);
        end
    end
    %read the xml file
    xDoc = xmlread(file);
end

%parse xDoc into a MATLAB structure
s = parseChildNodes(xDoc);

end

% ----- Subfunction parseChildNodes -----
function [children,ptext,textflag] = parseChildNodes(theNode)
% Recurse over node children.
children = struct;
ptext = struct; textflag = 'Text';
if hasChildNodes(theNode)
    childNodes = getChildNodes(theNode);
    numChildNodes = getLength(childNodes);
    
    for count = 1:numChildNodes
        theChild = item(childNodes,count-1);
        [text,name,attr,childs,textflag] = getNodeData(theChild);
        
        if (~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section'))
            %XML allows the same elements to be defined multiple times,
            %put each in a different cell
            if (isfield(children,name))
                if (~iscell(children.(name)))
                    %put existsing element into cell format
                    children.(name) = {children.(name)};
                end
                index = length(children.(name))+1;
                %add new element
                children.(name){index} = childs;
                if(~isempty(fieldnames(text)))
                    children.(name){index} = text;
                end
                if(~isempty(attr))
                    children.(name){index}.('Attributes') = attr;
                end
            else
                %add previously unknown (new) element to the structure
                children.(name) = childs;
                if(~isempty(text) && ~isempty(fieldnames(text)))
                    children.(name) = text;
                end
                if(~isempty(attr))
                    children.(name).('Attributes') = attr;
                end
            end
        else
            ptextflag = 'Text';
            if (strcmp(name, '#cdata_dash_section'))
                ptextflag = 'CDATA';
            elseif (strcmp(name, '#comment'))
                ptextflag = 'Comment';
            end
            
            %this is the text in an element (i.e., the parentNode)
            if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                    ptext.(ptextflag) = text.(textflag);
                else
                    %what to do when element data is as follows:
                    %<element>Text <!--Comment--> More text</element>
                    
                    %put the text in different cells:
                    % if (~iscell(ptext)) ptext = {ptext}; end
                    % ptext{length(ptext)+1} = text;
                    
                    %just append the text
                    ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                end
            end
        end
        
    end
end
end

% ----- Subfunction getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
% Create structure of node info.

%make sure name is allowed as structure name
name = toCharArray(getNodeName(theNode))';
name = strrep(name, '-', '_dash_');
name = strrep(name, ':', '_colon_');
name = strrep(name, '.', '_dot_');

attr = parseAttributes(theNode);
if (isempty(fieldnames(attr)))
    attr = [];
end

%parse child nodes
[childs,text,textflag] = parseChildNodes(theNode);

if (isempty(fieldnames(childs)) && isempty(fieldnames(text)))
    %get the data of any childless nodes
    % faster than if any(strcmp(methods(theNode), 'getData'))
    % no need to try-catch (?)
    % faster than text = char(getData(theNode));
    text.(textflag) = toCharArray(getTextContent(theNode))';
end

end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = struct;
if hasAttributes(theNode)
    theAttributes = getAttributes(theNode);
    numAttributes = getLength(theAttributes);
    
    for count = 1:numAttributes
        %attrib = item(theAttributes,count-1);
        %attr_name = regexprep(char(getName(attrib)),'[-:.]','_');
        %attributes.(attr_name) = char(getValue(attrib));
        
        %Suggestion of Adrian Wanner
        str = toCharArray(toString(item(theAttributes,count-1)))';
        k = strfind(str,'=');
        attr_name = str(1:(k(1)-1));
        attr_name = strrep(attr_name, '-', '_dash_');
        attr_name = strrep(attr_name, ':', '_colon_');
        attr_name = strrep(attr_name, '.', '_dot_');
        attributes.(attr_name) = str((k(1)+2):(end-1));
    end
end
end

%%-------------------------------------------------------------------------


function svel = sw_svel(S,T,P)

% SW_SVEL    Sound velocity of sea water
%=========================================================================
% SW_SVEL  $Revision: 1.3 $  $Date: 1994/10/10 05:53:00 $
%          Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE:  svel = sw_svel(S,T,P)
%
% DESCRIPTION:
%    Sound Velocity in sea water using UNESCO 1983 polynomial.
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78)]
%   T = temperature [degree C (IPTS-68)]
%   P = pressure    [db]
%       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT:
%   svel = sound velocity  [m/s]
%
% AUTHOR:  Phil Morgan 93-04-20  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%    Fofonoff, P. and Millard, R.C. Jr
%    Unesco 1983. Algorithms for computation of fundamental properties of
%    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%=========================================================================

% CALLER: general purpose
% CALLEE: none

% UNESCO 1983. eqn.33  p.46

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=3
    error('sw_svel.m: Must pass 3 parameters')
end %if

% CHECK S,T,P dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);
[mp,np] = size(P);


% CHECK THAT S & T HAVE SAME SHAPE
if (ms~=mt) | (ns~=nt)
    error('check_stp: S & T must have same dimensions')
end %if

% CHECK OPTIONAL SHAPES FOR P
if     mp==1  & np==1      % P is a scalar.  Fill to size of S
    P = P(1)*ones(ms,ns);
elseif np==ns & mp==1      % P is row vector with same cols as S
    P = P( ones(1,ms), : ); %   Copy down each column.
elseif mp==ms & np==1      % P is column vector
    P = P( :, ones(1,ns) ); %   Copy across each row
elseif mp==ms & np==ns     % PR is a matrix size(S)
    % shape ok
else
    error('check_stp: P has wrong dimensions')
end %if
[mp,np] = size(P);



% IF ALL ROW VECTORS ARE PASSED THEN LET US PRESERVE SHAPE ON RETURN.
Transpose = 0;
if mp == 1  % row vector
    P       =  P(:);
    T       =  T(:);
    S       =  S(:);
    
    Transpose = 1;
end %if
%***check_stp

%---------
% BEGIN
%--------

P = P/10;  % convert db to bars as used in UNESCO routines

%------------
% eqn 34 p.46
%------------
c00 = 1402.388;
c01 =    5.03711;
c02 =   -5.80852e-2;
c03 =    3.3420e-4;
c04 =   -1.47800e-6;
c05 =    3.1464e-9;

c10 =  0.153563;
c11 =  6.8982e-4;
c12 = -8.1788e-6;
c13 =  1.3621e-7;
c14 = -6.1185e-10;

c20 =  3.1260e-5;
c21 = -1.7107e-6;
c22 =  2.5974e-8;
c23 = -2.5335e-10;
c24 =  1.0405e-12;

c30 = -9.7729e-9;
c31 =  3.8504e-10;
c32 = -2.3643e-12;

Cw =    c00 + c01.*T + c02.*T.^2 + c03.*T.^3 + c04.*T.^4 + c05.*T.^5   ...
    + (c10 + c11.*T + c12.*T.^2 + c13.*T.^3 + c14.*T.^4).*P          ...
    + (c20 + c21.*T + c22.*T.^2 + c23.*T.^3 + c24.*T.^4).*P.^2        ...
    + (c30 + c31.*T + c32.*T.^2).*P.^3;

%-------------
% eqn 35. p.47
%-------------
a00 =  1.389;
a01 = -1.262e-2;
a02 =  7.164e-5;
a03 =  2.006e-6;
a04 = -3.21e-8;

a10 =  9.4742e-5;
a11 = -1.2580e-5;
a12 = -6.4885e-8;
a13 =  1.0507e-8;
a14 = -2.0122e-10;

a20 = -3.9064e-7;
a21 =  9.1041e-9;
a22 = -1.6002e-10;
a23 =  7.988e-12;

a30 =  1.100e-10;
a31 =  6.649e-12;
a32 = -3.389e-13;

A =     a00 + a01.*T + a02.*T.^2 + a03.*T.^3 + a04.*T.^4       ...
    + (a10 + a11.*T + a12.*T.^2 + a13.*T.^3 + a14.*T.^4).*P   ...
    + (a20 + a21.*T + a22.*T.^2 + a23.*T.^3).*P.^2            ...
    + (a30 + a31.*T + a32.*T.^2).*P.^3;


%------------
% eqn 36 p.47
%------------
b00 = -1.922e-2;
b01 = -4.42e-5;
b10 =  7.3637e-5;
b11 =  1.7945e-7;

B = b00 + b01.*T + (b10 + b11.*T).*P;

%------------
% eqn 37 p.47
%------------
d00 =  1.727e-3;
d10 = -7.9836e-6;

D = d00 + d10.*P;

%------------
% eqn 33 p.46
%------------
svel = Cw + A.*S + B.*S.*sqrt(S) + D.*S.^2;

if Transpose
    svel = svel';
end %if

end
%--------------------------------------------------------------------------
function applylf(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);


switch gdata.lfd.lftype
    case 1
        fun=@min;
    case 2
        fun=@max;
    case 3
        fun=@mean;
    case 4
        fun=@median;
end


if ~isfield(gdata,'eraw');
    gdata.eraw=gdata.bdata.zraw;
    gdata.ezc=gdata.bdata.zc;
end
if gdata.numedits==0
    set(gdata.menu12,'visible','on');
    set(gdata.menu13,'visible','on');
end
gdata.numedits=gdata.numedits+1;
if gdata.pan==1;
    set(gdata.toggle1,'value',0);
    pan off
end

if gdata.rawflag==1;
    pl = selectdata('sel','lasso','ignore',gdata.rim);
    
    
    
    if isfield(gdata,'sv')
        zr=slidefun(fun,gdata.lfd.lflen,...
            gdata.bdata.zraw(pl));
        if gdata.invert==1
            zr=-zr; %sv correction assumes depths postive
        end
      
        zsos=apply_sos_prof(zr,gdata.sv.sos_orig,...
            [gdata.sv.depth gdata.sv.sos],gdata.sv.use_mean_sos,...
            gdata.sv.mean_vel)-gdata.manoff;
        if gdata.invert==1
            zsos=-zsos;
            zr=-zr; %switch raw depths back to negative
        end
        
        gdata.bdata.zc(pl)=zsos-gdata.bdata.tide(pl);
        gdata.bdata.zraw(pl)=zr;
        
        
    else
        gdata.bdata.zraw(pl)=slidefun(fun,gdata.lfd.lflen,...
            gdata.bdata.zraw(pl));
        gdata.bdata.zc(pl)=((gdata.bdata.zraw(pl)-...
            gdata.bdata.tide(pl)).*gdata.invert)-gdata.manoff;
    end
    
    
    set(gdata.l1,'ydata',gdata.bdata.zraw);
    
    gdata.edits{gdata.numedits}=pl;
    
    
    guidata(hfig,gdata)
    return
end


handles=findobj(gdata.ax,'type','line');


if length(handles)==1
    [pl,xs,ys] = selectdata('sel','lasso'); %#ok
    %insert slidefun
    gdata.bdata.zc(pl)=slidefun(fun,gdata.lfd.lflen,...
        gdata.bdata.zc(pl));
    
    set(gdata.l1,'ydata',gdata.bdata.zc)
    
    ind=find(isfinite(gdata.bdata.zc));
    
    if gdata.alongflag==0
        gdata.xlimo=[min(gdata.bdata.distance(ind)),...
            max(gdata.bdata.distance(ind))];
        gdata.ylimo=[min(gdata.bdata.zc(ind)),...
            max(gdata.bdata.zc(ind))];
    else
        if ~isempty(gdata.topo)
            gdata.xlimo=[min([gdata.bdata.adist;gdata.topo.dist]),...
                max([gdata.bdata.adist;gdata.topo.dist])];
            gdata.ylimo=[min([gdata.bdata.zc;gdata.topo.z]),...
                max([gdata.bdata.zc;gdata.topo.z])];
        else
            gdata.xlimo=[min(gdata.bdata.adist(ind)),...
                max(gdata.bdata.adist(ind))];
            gdata.ylimo=[min(gdata.bdata.zc(ind)),...
                max(gdata.bdata.zc(ind))];
        end
    end
    
    waitfor(pl)
    
    gdata.edits{gdata.numedits}=pl;
    
    guidata(hfig,gdata);
    setFocus(hfig)
    
    
else
    
    if ~isempty(gdata.flag)
        gdata.xlims=get(gca,'xlim');
        gdata.ylims=get(gca,'ylim');
        hold off
        if gdata.alongflag==1
            gdata.l1=plot(gdata.bdata.adist,gdata.bdata.zc);
            hold on
            gdata.bb=plot(gdata.bdata.adist(gdata.flag==1),...
                gdata.bdata.zc(gdata.flag==1),'o',...
                'color',[0.6 0.6 0.6],'markersize',3,...
                'markerfacecolor',[0.6 0.6 0.6]);
            gdata.gg=plot(gdata.adistc,gdata.bdata.zc,'r-','linewidth',2);
        else
            
            gdata.l1=plot(gdata.bdata.distance,gdata.bdata.zc);
            hold on
            gdata.bb=plot(gdata.bdata.distance(gdata.flag==1),...
                gdata.bdata.zc(gdata.flag==1),'o',...
                'color',[0.6 0.6 0.6],'markersize',3,...
                'markerfacecolor',[0.6 0.6 0.6]);
            gdata.gg=plot(gdata.distc,gdata.cdata.zc,...
                'r-','linewidth',2);
        end
        
        
        
        set(gdata.ax,'xlim',gdata.xlims,'ylim',gdata.ylims);
        
        c1=legend('Raw Data','Outliers','Filtered Data');
        set(c1,'box','off','location','best')
        
        if ~isempty(gdata.topo);
            gdata.tlineh=plot(gdata.topo.dist,gdata.topo.z,'r-',...
                'linewidth',2);
            [pl,xs,ys] = selectdata('sel','lasso',...
                'ignore',[gdata.gg gdata.bb gdata.tlineh]); %#ok
        else
            
            [pl,xs,ys] = selectdata('sel','lasso',...
                'ignore',[gdata.gg gdata.bb]); %#ok
        end
        
        gdata.bdata.zc(pl)=slidefun(fun,gdata.lfd.lflen,...
            gdata.bdata.zc(pl));
        
        set(gdata.l1,'ydata',gdata.bdata.zc)
        gdata.edits{gdata.numedits}=pl;
        
        guidata(hfig,gdata);
        
        waitfor(pl)
        xyzFilt(hfig);
    else
        gdata.xlims=get(gca,'xlim');
        gdata.ylims=get(gca,'ylim');
        
        hold off
        gdata.l1=plot(gdata.bdata.adist,gdata.bdata.zc);
        hold on
        gdata.tlineh=plot(gdata.topo.dist,gdata.topo.z,'r-',...
            'linewidth',2);
        set(gdata.ax,'xlim',gdata.xlims,'ylim',gdata.ylims);
        
        [pl,xs,ys] = selectdata('sel','lasso',...
            'ignore',[gdata.tlineh]); %#ok
        
        gdata.bdata.zc(pl)=slidefun(fun,gdata.lfd.lflen,...
            gdata.bdata.zc(pl));
        
        gdata.edits{gdata.numedits}=pl;
        set(gdata.l1,'ydata',gdata.bdata.zc)
        guidata(hfig,gdata);
    end
    
end



setFocus(hfig)

end


%--------------------------------------------------------------------------
function filt_local(hfig,evnt) %#ok

gdata=guidata(hfig);

hf = figure('units','normalized','position',[0.253 0.521 0.145 0.104],...
    'menubar','none','name','Local Filter Options',...
    'numbertitle','off','color',[0.925 0.914 0.847]);
uicontrol(hf,'style','text','units',...
    'normalized','position',[0.0707 0.233 0.374 0.18],...
    'string','Filter Type','backgroundcolor',[0.925 0.914 0.847]);
uicontrol(hf,'style','text','units','normalized',...
    'position',[0.037 0.639 0.407 0.158],'string','Window Length',...
    'backgroundcolor',[0.925 0.914 0.847]);

lfd.edit1 = uicontrol(hf,'style','edit','units','normalized',...
    'position',[0.502 0.233 0.418 0.211],'string',...
    sprintf('%0.0f',gdata.lfd.lflen),...
    'backgroundcolor',[1 1 1]);

lfd.popupmenu1 = uicontrol(hf,'style','popupmenu',...
    'units','normalized','position',[0.508 0.617 0.407 0.203],...
    'string',{'Min';'Max';'Mean';'Median'},'backgroundcolor',[1 1 1],...
    'value',gdata.lfd.lftype);

uicontrol(hf,'style','pushbutton',...
    'units','normalized',...
    'position',[0.4 0.01 0.2 0.15],...
    'string','Done','callback',@lfclose)

guidata(hf,lfd);

uiwait
lfd=guidata(hf);
gdata.lfd=lfd;
guidata(hfig,gdata);
close(hf)
end
%%%%%---------------------------------------------------------------------
function lfclose(hf,evnt) %#ok

lfd=guidata(hf);

lfd.lftype=get(lfd.popupmenu1,'value');
lfd.lflen=str2double(get(lfd.edit1,'string'));

if rem(lfd.lflen,1)~=0
    set(lfd.edit1,'string','Enter Integer Value',...
        'backgroundcolor','r');
    return
else
    set(lfd.edit1,'backgroundcolor','w')
end

guidata(hf,lfd);
uiresume
end

%%%%-----------------------------------------------------------------------
function runEchoFix(hfig,eventdata,handles) %#ok
gdata=guidata(hfig);

if gdata.fixlen~=1
    %reload the original file
    if gdata.invert==1
        [gdata.bdata,gdata.hdr,gdata.rdata]=readRAW(...
            'filename',[gdata.filepath,gdata.fnames{gdata.newInd}],...
            'invert','minz',gdata.minz,'maxz',gdata.maxz);
    else
        [gdata.bdata,gdata.hdr,gdata.rdata]=readRAW(...
            'filename',[gdata.filepath,gdata.fnames{gdata.newInd}],...
            'minz',gdata.minz,'maxz',gdata.maxz);
    end
end

prompt={'Echofix window length:'};
name='Echofix Parameter';
numlines=1;
defaultanswer={num2str(gdata.fixlen)};

answer=inputdlg(prompt,name,numlines,defaultanswer);
gdata.fixlen=str2double(answer{1});

gdata.bdata=echofix(gdata.fixlen,gdata.bdata);
set(gdata.l1,'ydata',gdata.bdata.zc)


guidata(hfig,gdata);

end
%%%%%--------------------------------------------------------------------
function df=echofix(wlen,varargin)
%ECHOFIX-fix latency problem on hydrobox echosounder

%read in data file
if nargin<2
    [d,h,r]=readRAWn2; %#ok
elseif nargin~=2
    error('Unknown input options');
else
    d=varargin{1};
end

if wlen==1
    df=d;
    return
end

ind=find(~isnan(d.tide) & ~isnan(d.zraw));
pn=(1:numel(d.zraw))';
fz=ones(numel(d.zraw),1)*NaN;
fh=ones(numel(d.tide),1)*NaN;
rz=d.zraw(ind);
rh=d.tide(ind);


%break data into chunks of width= wlen
nchunks=floor(numel(rz)/wlen);
leftover=rem(numel(rz),wlen);

zmat=num2cell(reshape(rz(1:end-leftover),...
    [wlen nchunks]),1);
tmat=num2cell(reshape(rh(1:end-leftover),...
    [wlen nchunks]),1);
pmat=num2cell(reshape((1:numel(rz)-leftover),...
    [wlen nchunks]),1);

%deal with the leftover
if leftover~=0
    zmat(end+1)=num2cell(rz((end-leftover+1):end),1);
    tmat(end+1)=num2cell(rh((end-leftover+1):end),1);
    pmat(end+1)=num2cell([(numel(rz)-leftover+1):...
        (numel(rz))]',1);
end



%detrend data by removing mean from each chunk
zdt=cellfun(@(x)(detrend(x)),zmat,'uni',0);
tdt=cellfun(@(x)(detrend(x)),tmat,'uni',0);


%plot relative antenna height from gps
%and wiggles in echosounder data
% figure
% ax(1)=subplot(211);
% hold on
% h1=cellfun(@(x,y)(plot(x,y,'k')),pmat,zdt);
% h2=cellfun(@(x,y)(plot(x,y,'b')),pmat,tdt);

%find correlation for each chunk between the
%two data streams
maxlag=100;
[c,lags]=cellfun(@(x,y)(xcorr(x,y,maxlag)),...
    zdt,tdt,'uni',0);
c2=cellfun(@(x)(x(maxlag+1:end)),c,'uni',0);
lag2=cellfun(@(x)(x(maxlag+1:end)),lags,'uni',0);
lag=cellfun(@(x,y)(x(find(y==max(y),1,'first'))),...
    lag2,c2,'uni',0);

%apply calculated lag to echosounder data stream

zhc=cell(1,length(tmat));
zhc1={[ones(lag{1},1)*NaN;tmat{1}(1:end-lag{1})]};
zhc2=cellfun(@(x,y,z)([x(end-(z-1):end);y(1:end-z)]),...
    tmat(1:end-1),tmat(2:end),lag(2:end),'uni',0);
zhc(1)=zhc1;
zhc(2:end)=zhc2;


tc=cell(1,length(tmat));
tc1={[ones(lag{1},1)*NaN;tdt{1}(1:end-lag{1})]};
tc2=cellfun(@(x,y,z)([x(end-(z-1):end);y(1:end-z)]),...
    tdt(1:end-1),tdt(2:end),lag(2:end),'uni',0);
tc(1)=tc1;
tc(2:end)=tc2;



%convert cells back to vector
tcmat=cell2mat(tc');
zhcmat=cell2mat(zhc');
pnmat=cell2mat(pmat');
znmat=cell2mat(zmat');

%plot corrected antenna height
% h3=plot(pnmat,tcmat,'r');
% legend([h1(1);h2(1);h3],'Echosounder Wiggles',...
%     'Raw detrended gps','Corrected gps')

zclean=znmat-zhcmat;


%plot raw and corrected profile
% ax(2)=subplot(212);
% plot([1:numel(d.zc)],d.zc,'k')
% hold on
% plot(pn(ind),zclean,'r')
% legend('Raw','Corrected')

% linkaxes(ax,'x')

df=d;
fz(ind)=zclean;
df.zc=fz;
fh(ind)=zhcmat;
df.tide=fh;

if nargout>1
    varargout{1}=cell2mat(lag);
end

end


%%%%-----------------------------------------------------------------------
function [bdata,hdr,varargout]=readRAW(varargin)
%READRAW - reads hypack *.RAW files
%
%   [bdata,hdr] = readRAW reads in Hypack's native
%   raw data format and returns a structure with the data,
%   bdata, and pertinent header information, hdr.  This code
%   applies specifically to data files collected with the
%   USGS Coastal Profiling System (assumes a gps and echosounder).
%
%   NOTE: Because typically, there are more soundings than gps fixes,
%   the positions were interpolated to provide a position for each
%   sounding.  Additionally, only positions that were flagged as
%   being "Fix-RTK" were accepted.
%
%   INPUTS: None Required, but you can specify the data file
%   using the following sytax:
%
%       [bdata,hdr]=readRAW('filename','foo.raw') or
%       [bdata,hdr]=readRAW('filename','C:\bathy\foo.raw')
%
%   If no inputs are supplied, the user will be prompted
%   to supply file interactively. Additionally, if you want to see
%   raw (non-interpolated) positions and gps quality info,
%   supply a third output to the call to readHypackRaw:
%
%       [bdata,hdr,gps]=readRAW;
%
%   You can invert z values using the following syntax:
%
%       [bdata,hdr,gps]=readRAW('invert');
%
%   OUTPUTS:
%   bdata- as structure with the following fields:
%       mtime: matlab datenum for each sounding
%       x: x-position for each sounding (in coordinate system
%          supplied to Hypack.
%       y: y-position for each sounding (in coordinate system
%          supplied to Hypack.
%       distance: along-transect distance.
%       offline: if the line followed was straight (ie. planned line has
%           2 points) the offline distance is calculated.
%       lon: longitude
%       lat: latitude
%       tide: tide correction applied to raw soundings
%       zraw: raw sounding data (range, m)+ offset (see hdr).
%       zc: tide-corrected sounding data (m).
%
%   hdr- a structure with the following fields:
%       linefile: linefile used during data collection.
%       ellipsoid: ellipsoid supplied to Hypack.
%       projection: information about the projection system.
%       datumtrans: information about any datum transformations.
%       datestr: time and date when data collection began.
%       mtime: datenum when data collection began.
%       dev_(1...n): driver used for any attatched devices.
%       dev_(1...n)_offsets: offsets applied to any attached devices.
%       xstart: starting x-position (projected) of data collection.
%       ystart: starting y-position (projected) of data collection.
%       linenum: planned line number.
%       linex: x-coordinates of planned line.
%       liney: y-coordinates of planned line.
%
%   rdata- a structure with the following fields
%       ptime: datenum of gps fixes.
%       x: raw x-position (projected).
%       y: raw y-position (projected).
%       hdop: gps hdop for each fix.
%       numsats: number of satellites for each gps fix.
%       gpsmode: type of gps fix (0=invalid, 1= stand alone, 2=
%           differential, 3 and 4 = fix RTK, 5= RTK float).
%       lat: latitude.
%       lon: longitude.
%       antennaH: gps antenna height (m).
%       gpsTime: utc time.
%       tide: tide correction appled to raw soundings (m).
%
% A.Stevens @ USGS, 7/17/2007
% astevens@usgs.gov

%parse inputs
fpath=[];
invert=1;
minz=-inf;
maxz=inf;
maxoff=inf;
read_raw_msg=0;

if nargin>0
    [m1,n1]=size(varargin); %#ok
    opts={'filename','invert','minz','maxz','maxoffset','readmsg'};
    
    for i=1:n1;
        indi=strcmpi(varargin{i},opts);
        ind=find(indi==1);
        if isempty(ind)~=1
            switch ind
                case 1
                    fpath=varargin{i+1};
                    if exist(fpath,'file')==0
                        errordlg('Path not found. Try again');
                        return
                    end
                case 2
                    invert=-1;
                case 3
                    minz=varargin{i+1};
                case 4
                    maxz=varargin{i+1};
                case 5
                    maxoff=varargin{i+1};
                case 6
                    read_raw_msg=1;
            end
        else
        end
    end
end


if isempty(fpath)==1;
    [filename, pathname] = ...
        uigetfile('*.RAW', 'Pick an RAW-file');
    fpath=[pathname,filename];
end


%open file and check its length
fid=fopen(fpath);
%read header
numDev=0;
eoh = 0;
while eoh~=1
    l = fgetl(fid); %read line by line
    s1 = strread(l,'%s')';
    
    switch s1{1}
        case 'VER'
            hdr.hypack_ver=s1{2};
        case 'FIL'
            hdr.linefile=s1{3};
            
        case 'ELL' %read ellipsoid data
            hdr.ellipsoid = s1(2:end);
            
        case 'PRO' %projection information
            hdr.projection = s1(2:end);
            
        case 'GEO'
            s2=s1(2:end);
            for i=1:length(s2);
                if ~isempty(regexp(s2{i},'.geo','once'))
                    [junk,file,ext]=fileparts(s2{i}); %#ok
                    
                    hdr.geoid=[file,ext];
                end
            end
            hdr.orthometric_height_corr=s2{end};
            
            
        case 'DTM' %datum transformation
            hdr.datumtrans = s1(2:end);
            
        case 'TND' %time of data collection
            hdr.datestr = s1(2:end);
            [hr,mini,sec] = strread(s1{2},'%n:%n:%n');
            [mon,day,yr] = strread(s1{3},'%n/%n/%n');
            hdr.mtime=datenum(yr,mon,day,hr,mini,sec);
            
        case 'DEV' %read device and driver info
            numDev=numDev+1;
            namer=['dev_',num2str(numDev)];
            hdr.(namer)=s1(2:end);
            
        case 'OFF' %read offset info
            namer=['dev_',num2str(numDev),'_offsets'];
            hdr.(namer)=s1(2:end);
            
        case 'LIN' % planned line points
            numpoints=str2double(s1{2});
            xp=zeros(numpoints,1);
            yp=zeros(numpoints,1);
            
            for i=1:numpoints;
                lp=fgetl(fid);
                [junk,xp(i),yp(i)]=strread(lp,'%s %f %f'); %#ok
            end
            
        case 'LBP'
            hdr.xstart=str2double(s1{2});
            hdr.ystart=str2double(s1{3});
            
            
        case 'LNN' %planned line number
            lineNum=str2double(s1{2});
            hdr.lineNum=lineNum;
            hdr.linex=xp;
            hdr.liney=yp;
            
        case 'EOH'
            eoh=1;
            
    end
end

% read data
% split the file into the data code
% and the string of data

data = textscan(fid,'%s%[^\n]');

%start at first pos code string
firstP=find(strcmpi('POS',data{1})==1,1,'first');
data{1}=data{1}(firstP:end);
data{2}=data{2}(firstP:end);

%define codes and associated formats
codes={'EC2';'EC1';'POS';'QUA';'RAW';'TID';'KTC'};
formats={'%*d %f %f';...
    '%*d %f %f';...
    '%*d %f %f %f %*f';...
    ['%*d %*f %*d %*f %f %f %f',...
    '%*d %*d %*d %*f %*d'];...
    '%*d %*f %*d %f %f %f %f';...
    '%*d %f %f';...
    ['%*d',repmat(' %f',1,9)]};


for i=1:length(codes);
    
    %for each code, turn strings of data into
    %a matrix...kinda cool way to handle complex
    %files with many formats without loops
    data2={data{2}{strcmpi(codes{i},data{1})}}';
    if ~isempty(data2)
        r=cellfun(@(x) sscanf(x,formats{i})',...
            data2,'uni',false);
        nvals=cellfun(@numel,r);
        
        
        
        
%         data3=cell2mat(cellfun(@(x)(x(...
%             1:min(nvals(nvals~=0)))),r(nvals~=0),'un',0));
        data3=cell2mat(r(nvals==max(nvals)));
        

        %based on each of the codes, define variables
        %and format data into a structure
        switch codes{i}
            
            case 'EC1' %sounder data
                
                
                ztime=data3(:,1);
                sdata.vtime=ztime;
                sdata.ztime=...
                    datenum(yr,mon,day)+ztime./86400;
                sdata.zraw=data3(:,2);
                
                
            case 'EC2' %sounder data
                
                ztime=data3(:,1);
                sdata.vtime=ztime;
                sdata.ztime=...
                    datenum(yr,mon,day)+ztime./86400;
                sdata.zraw=data3(:,2);
                
                
            case 'POS'  %read position
                
                rdata.hygtime=data3(:,1);
                rdata.ptime=...
                    datenum(yr,mon,day)+...
                    rdata.hygtime./86400;
                rdata.x=data3(:,2);
                rdata.y=data3(:,3);
                
            case 'QUA' %gps quality
                
                rdata.hdop=data3(:,1);
                rdata.numsats=data3(:,2);
                rdata.gpsmode=data3(:,3);
                
            case 'RAW' %raw gps data
                
                lon=data3(:,2);
                lat=data3(:,1);
                rdata.antennaH=data3(:,3);
                rdata.gpsTime=data3(:,4);
                
                
                lat1=fix(lat/10000);
                lat2= (lat-(lat1*10000))/6000;
                rdata.lat=lat1+lat2;
                
                lon1=fix(lon/10000);
                lon2= (lon-(lon1*10000))/6000;
                rdata.lon=lon1+lon2;
                
            case 'TID'
                
                tdata.hyttime=data3(:,1);
                tdata.tide=data3(:,2);
                tdata.ttime=...
                    datenum(yr,mon,day)+...
                    tdata.hyttime./86400;
                
            case 'KTC'
                kdata.ktime=data3(:,1);
                kdata.ellipsoid_height_wgs84=...
                    data3(:,3);
                kdata.ellipsoid_height_local=...
                    data3(:,4);
                kdata.undulation=data3(:,5);
                kdata.k_value=data3(:,6);
                kdata.antenna_offset=data3(:,7);
                kdata.draft_correction=data3(:,8);
                
                
                
        end
    end
    
end

if read_raw_msg
    data2={data{2}{strcmpi('MSG',data{1})}}';
    fmt='%*d %f %[^\n]';
    gpsdata=cellfun(@(x)(textscan(x,fmt,...
        'delimiter',' ')),data2,'un',0);
    vtime=cell2mat(cellfun(@(x)(x(1)),gpsdata));
    gpsstr=cellfun(@(x)(x(2)),gpsdata);
    gpsstr=cat(1,gpsstr{:});
    
    
    msg=readmsg(gpsstr);
    msg.htime=vtime;
end




fclose(fid);

if ~exist('sdata','var') || ~exist('rdata','var')
    bdata=[];
    if nargout==3
        varargout={[]};
    end
    return
end

%done reading data, now do some processing
%interp values so there is a position for each sounding

% ind=find(rdata.gpsmode==3 | rdata.gpsmode==4); %use only fix-rtk
% positions..
% ind=find(rdata.gpsmode==3); %use only fix-rtk positions..
% if isempty(ind);
%     ind=1:length(rdata.gpsmode);
% end
% rdata=structfun(@(x)(x(ind)),rdata,'un',0);

[pathname,filename,ext]=fileparts(fpath);%#ok

bdata.filename=[filename,ext];
bdata.mtime=sdata.ztime;
bdata.vtime=sdata.vtime;
[btime,bind]=unique(rdata.ptime);


bdata.x=interp1(btime,rdata.x(bind),sdata.ztime);
bdata.y=interp1(btime,rdata.y(bind),sdata.ztime);

%calculate cumulative distance along line
x1=bdata.x(1:end-1);
x2=bdata.x(2:end);
y1=bdata.y(1:end-1);
y2=bdata.y(2:end);

xd = x2 - x1;
yd = y2 - y1;
dist=sqrt(xd.*xd +yd.*yd);
dist(isnan(dist))=0;
bdata.distance=[0,cumsum(dist)'];


% if it is a straight line followed, calculate offline distance
if isfield(hdr,'linex');
    if length(hdr.linex)==2
        m=(hdr.liney(end)-hdr.liney(1))/...
            (hdr.linex(end)-hdr.linex(1));
        b=hdr.liney(1)-(m*hdr.linex(1));
        yi=m.*bdata.x+b;
        
        %determine distance from line to points
        h=bdata.y-yi;
        ang=(pi/2)-abs(atan(m));
        d=h*sin(ang);
        bdata.offline=d;
        
        
        %determine along-line distance
        
        theta=atan(m);
        xt = bdata.x.*cos(theta) + bdata.y.*sin(theta);
        yt = -bdata.x.*sin(theta) + bdata.y.*cos(theta);
        x0=hdr.linex.*cos(theta) + ...
            hdr.liney.*sin(theta);
        
        %calculate along line distance
        bdata.adist=xt-x0(1);
    else
        bdata.offline=[];
        bdata.adist=[];
    end
    
end

if isfield(rdata,'lon');
    bdata.lon=interp1(btime,rdata.lon(bind),sdata.ztime);
    bdata.lat=interp1(btime,rdata.lat(bind),sdata.ztime);
end

sdata.zraw(sdata.zraw<minz | sdata.zraw>maxz)=NaN;
if ~isempty(bdata.offline);
    sdata.zraw(abs(bdata.offline)>maxoff)=NaN;
end;

bdata.zraw=(sdata.zraw).*invert; %invert z if specified

if exist('tdata','var')
    
    
    if exist('kdata','var')
        fields=fieldnames(kdata);
        [junk,ia,ib]=intersect(rdata.hygtime,kdata.ktime); %#ok
        for i=1:length(fields)
            
            rdata.(fields{i})=nan(numel(rdata.ptime),1);
            rdata.(fields{i})(ia)=kdata.(fields{i})(ib);
        end
    end
    
    [junk,ia,ib]=intersect(rdata.hygtime,tdata.hyttime); %#ok
    rdata.tide=nan(numel(rdata.ptime),1);
    rdata.tide(ia)=tdata.tide(ib);
    
    %only use values with valid tide data
    warning('off','MATLAB:interp1:NaNinY')
    try
        bdata.tide=interp1(rdata.ptime,rdata.tide,sdata.ztime);
        bdata.zc=(sdata.zraw+bdata.tide).*invert; %apply tide correction
    catch %#ok
        [ttime,tind]=unique(tdata.ttime);
        tide=rdata.tide(tind);
        bdata.tide=interp1(ttime,tide,sdata.ztime);
        bdata.zc=(sdata.zraw+bdata.tide).*invert; %apply tide correction
    end
end



%if called output raw position and gps quality data
if nargout==3
    varargout{1}=rdata;
end
if nargout==4;
    varargout{1}=rdata;
    varargout{2}=msg;
end

end
%%------------------------------------------------------------------------------------------------

function msg= readmsg(str)

data = cellfun(@(x)(textscan(x,'%s%[^\n]',...
    'delimiter',',')),str,'un',0);

id=cellfun(@(x)(x(1)),data);
id=cat(1,id{:});
data2=cellfun(@(x)(x(2)),data);
data2=cat(1,data2{:});

%define codes and associated formats
codes={'$GPRMC';'$SDDPT';'$PTNL'};
formats={['%s %*s %f %*s %f %*s',...
    '%f %f %s %*[^\n]'];...
    '%f %*[^\n]';...
    ['%*s %s %s %f %*s %f',...
    '%*s %f %f %f %3s%f %*[^\n]']};


for i=1:length(codes);
    
    data3=data2(strcmpi(codes{i},id));
    if ~isempty(data3)
        r=cellfun(@(x)(textscan(x,formats{i},...
            'delimiter',',')),...
            data3,'un',0);
        
        switch codes{i}
            case '$GPRMC'
                %time
                times=cellfun(@(x)(textscan(char(x{1}),...
                    '%2.0f%2.0f%2.0f')),r,'un',0);
                dates=cellfun(@(x)(textscan(char(x{6}),...
                    '%2.0f%2.0f%2.0f')),r,'un',0);
                tmat=cell2mat(cat(2,fliplr(cat(1,dates{:})),...
                    cat(1,times{:})));
                dn=datenum(tmat(:,1)+2000,tmat(:,2),tmat(:,3),...
                    tmat(:,4),tmat(:,5),tmat(:,6));
                %lat
                lat=cellfun(@(x)(x{2}),r);
                lat1=fix(lat/100);
                lat2= (lat-(lat1*100))/60;
                latd=lat1+lat2;
                
                %lon
                lon=cellfun(@(x)(x{3}),r);
                lon1=fix(lon/100);
                lon2= (lon-(lon1*100))/60;
                lond=-(lon1+lon2);
            case '$SDDPT'
                depth=cell2mat(cat(1,r{:}));
                
            case '$PTNL'
                times=cellfun(@(x)(textscan(char(x{1}),...
                    '%2.0f%2.0f%2.0f.%2.0f',1)),r,'un',0);
                dates=cellfun(@(x)(textscan(char(x{2}),...
                    '%2.0f%2.0f%2.0f')),r,'un',0);
                tmat=cell2mat(cat(2,cat(1,dates{:}),...
                    cat(1,times{:})));
                dn=datenum(tmat(:,3)+2000,tmat(:,1),tmat(:,2),...
                    tmat(:,4),tmat(:,5),tmat(:,6)+(tmat(:,7)/100));
                %lat
                lat=cellfun(@(x)(x{3}),r);
                lat1=fix(lat/100);
                lat2= (lat-(lat1*100))/60;
                latd=lat1+lat2;
                
                %lon
                lon=cellfun(@(x)(x{4}),r);
                lon1=fix(lon/100);
                lon2= (lon-(lon1*100))/60;
                lond=-(lon1+lon2);
                
                elev=cellfun(@(x)(x{9}),r);
        end
    end
end


%collect output
msg.mtime=dn;
msg.lon=lond;
msg.lat=latd;
if exist('elev','var')
    msg.elev=elev;
end
if exist('depth','var');
    msg.depth=depth;
end




end
%%%%----------------------------------------------------------------------
function [lnw,varargout] = readLNW(varargin)
%READLNW -read Hypack line file(.lnw)
%
%   lnw = readLNW reads the formatted text file that
%   HYPACK uses for planning lines.  The function
%   requires no inputs.
%
%   OUTPUT:
%       line data are read into a structure
%       with fields name (line number), x and y.
%
%   USAGE OPTIONS:
%       lnw=readLNW('filename',example.lnw)- reads
%           in specified file rather than prompting
%           the user.
%
%       lnw=readLNW('plotit')- plots lines with
%           labels.
%
%       [lnw, sumstats]=readLNW - outputs a structure
%           summarizing line lengths, total distance
%           for planning purposes.
%
% A. Stevens 3/22/2007
% astevens@usgs.gov

%check for input file
fpath=[];
if nargin>=1
    ind=strcmpi(varargin,'filename');
    if any(ind)==1
        fpath=varargin{ind+1};
        if exist(fpath,'file')==0;
            erstr=['File does not exist in MATLAB path. ',...
                'Check your file and path and try again...'];
            error(erstr)
        end
    end
end

if isempty(fpath)==1;
    [filename, pathname] = ...
        uigetfile('*.lnw', 'Pick an LNW-file');
    fpath=[pathname,filename];
end

%read file
fid=fopen(fpath,'r');
l1=fgetl(fid);
[junk,numlines]=strread(l1,'%s %n'); %#ok

lind=1;
while lind<=numlines;
    ln=fgetl(fid);
    [nmr,numpoints]=strread(ln,'%s %n');
    if strcmpi(nmr,'LIN')==1;
        for i=1:numpoints;
            lp=fgetl(fid);
            [junk,x(i),y(i)]=strread(lp,'%s %f %f'); %#ok
        end
        lend=fgetl(fid);
        [junk,lname]=strread(lend,'%s %s'); %#ok
        
        
        lnw(lind).name=char(lname);
        lnw(lind).x=x;
        lnw(lind).y=y;
        
        lind=lind+1;
        clear x y lname
    end
end

%make a simple plot if specified
if any(strcmpi(varargin,'plotit'))==1;
    figure
    hold on
    
    %make ~50 labels per plot
    int=ceil(numlines/50);
    
    for i=1:numlines;
        plot(lnw(i).x,lnw(i).y);
        if rem(i,int)==0
            text(lnw(i).x(1),lnw(i).y(1),...
                num2str(lnw(i).name),'fontsize',10)
        end
    end
    axis equal
end

if nargout>1
    
    stats.numlines=numlines;
    for i=1:numlines
        x1=lnw(i).x(1:end-1);
        x2=lnw(i).x(2:end);
        y1=lnw(i).y(1:end-1);
        y2=lnw(i).y(2:end);
        
        xd = x2 - x1;
        yd = y2 - y1;
        distr(i)=sum(sqrt(xd.*xd +yd.*yd))';
    end
    stats.totalLength=sum(distr);
    stats.lineLens=distance;
    varargout{1}=stats;
end

end

%%%%----------------------------------------------------------------------
function R = slidefun (FUN, W, V, windowmode, varargin)
% SLIDEFUN - apply function to a moving window over a vector
%
%   R = SLIDEFUN(FUN, W, V) evaluates the function FUN to a moving window
%   of W consecutive elements of the vector V. The function FUN is
%   specified by a function handle or a function name and should return a
%   scalar for a vector input. W specifies the size of the window and
%   should be a positive integer. At the two edges of the vector less than
%   W elements will be used. R will have the same size as V.
%
%   Effectively SLIDEFUN applies the function FUN to a moving window of
%   consecutive elemens V(x0:(x0+W)) using FEVAL, and returns the result in
%   R(x).
%
%   The window [x0:(x0+W)] is positioned relative to the current point x in V.
%   R = SLIDEFUN(FUN, W, V, WINDOWMODE) denotes the type of windowing being
%   used. WINDOWMODE can be (the first letters of) one the following:
%       - 'central', or '' (default): the window is centered around each point,
%          so that x0 equals x - floor(W/2);
%       - 'backward': window is using W points before the current point, so
%          that x0 equals x-W+1
%       - 'forward': window is using W points following the current point,
%          so that x0 equals x
%
%  R = SLIDEFUN(FUN,W,V,WINDOWMODE, P1,P2,...) provides for aditional
%  parameters which are passed to the function FUN.
%
%   Example 1) Sliding max filter - return the maximum of every three
%   consecutive elements:
%      V =  [1  2  3  9  4  2  1  1  5  6] ;
%      R = slidefun(@max, 3, V)
%      % -> [2  3  9  9  9  4  2  5  6  6]
%      % So R(i) = max(R(i-1:i+1)) ;
%      % and R(1) = max(V(1:2))
%
%   Example 2) Sum every four consecutive elements
%      V =  [1  2  3  4  3  2  1] ;
%      R = slidefun('sum',4, V)
%      % -> [3  6 10 12 12 10  6]
%      % So R(i) = sum(R(i-2:i+1)) ;
%      %    R(1) = sum([1 2]) ;
%      %    R(2) = sum([1 2 3]) ;
%
%   Example 3) Range of every three consecutive elements
%      myfun = inline('max(x) - min(x)','x') ; % (Matlab R13)
%      V =  [1  4  3  3  3  2  9  8] ;
%      R = slidefun(myfun,3, V)
%      % -> [3  3  1  0  1  7  7  1]
%
%   Example 4) Mimick cumsum
%      V = 1:10 ;
%      R = slidefun(@sum, numel(V), V, 'backward')
%      isequal(R,cumsum(V))
%
%   Example 5) Inverse cumprod ignoring zeros
%      V = [1:3 0 5:8] ;
%      myfun = inline('prod(x(x~=0))','x') ;
%      R = slidefun(myfun, numel(V), V, 'forward')
%
%   Note that for some specific functions (e.g., MEAN) filter can do the
%   same job faster. See FEVAL for more information about passing
%   functions.
%
%   See also FEVAL, INLINE, FILTER, BOOTSTRP
%   and Function handles, Anonymous functions.

% Written and tested in Matlab R13
% version 3.0 (oct 2006)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% 1.0 (sep 2006). This file was inspired by a post on CSSM in sep 2006.
% 2.0 (oct 2006). Use for-loop instead of large matrices
% 3.0 (oct 2006). Added windowmode option (after File 9428 by John
%                 D'Errico)

% check input arguments,expected
% <function name>, <window size>, <vector>, <windowmode>, <optional arguments ...>
narginchk(3,Inf) ;

if nargin==3 || isempty(windowmode),
    windowmode = 'central' ;
end

% based on code by John D'Errico
if ~ischar(windowmode),
    error('WindowMode should be a character array') ;
else
    validmodes = {'central','backward','forward'} ;
    windowmode = strmatch(lower(windowmode), validmodes) ;
    if isempty(windowmode),
        error('Invalid window mode') ;
    end
end
% windowmode will 1, 2, or 3

if (numel(W) ~= 1) || (fix(W) ~= W) || (W < 1),
    error('Window size W must be a positive integer scalar.') ;
end

nV = numel(V) ;

if isa(FUN,'function_handle')
    FUN = func2str(FUN);
end
if ~exist(FUN),
    error('The function "%s" can not be found.',FUN) ;
end

if nV==0,
    % trivial case
    R = V ;
    return
end

% can the function be applied succesfully?
try
    R = feval(FUN,V(1:min(W,nV)),varargin{:}) ;
    % feval did ok. Now check for scalar output
    if numel(R) ~= 1,
        error('Function "%s" does not return a scalar output for a vector input.', FUN) ;
    end
catch
    % Rewrite the error, likely to be caused by feval
    % For instance, function expects more arguments, ...
    
    ERR = lasterror ;
    ERR.message = strrep(ERR.message,['feval' 10],sprintf('%s%c(FEVAL error): ',mfilename,10)) ;
    if numel(varargin)>0,
        ERR.message = sprintf('%s\r(This could be caused by the additional arguments.)',ERR.message) ;
    end
    rethrow(ERR) ;
end % try-catch

% where is the first relative element
switch windowmode
    case 1 % central
        x0 = -floor(W/2) ;
    case 2 % backward
        x0 = -(W-1) ;
    case 3 % forward
        x0 = 0 ;
end
x1 = x0+W-1 ; % last relative element
x = x0:x1 ; % window vector (has W elements)

R = zeros(size(V)) ; % pre-allocation !!

% The engine: seperation in three sections is faster than using a single
% loop with calls to min and max.

% 1. leading elements
for i=1:-x0,
    R(i) = feval(FUN,V(1:i+x1),varargin{:}) ;
end
% 2. main portion of V
for i=(-x0+1):(nV-x1),
    R(i) = feval(FUN,V(x+i),varargin{:}) ;
end
% 3. trailing elements
for i=(nV-x1+1):nV,
    R(i) = feval(FUN,V((i+x0):nV),varargin{:}) ;
end

% Some timings:
%   V = rand(10000,1) ; N = 3 ; tic ; R = slidefun(@max,N, V) ; toc ;
%   % -> elapsed_time = 0.1870
%   V = rand(10000,1) ; N = 100 ; tic ; R = slidefun(@max,N, V) ; toc ;
%   % -> elapsed_time = 0.3440
end

%%%%----------------------------------------------------------------------
function y = nanmean(x)
%NANMEAN Average or mean ignoring NaNs.
%   NANMEAN(X) returns the average treating NaNs as missing values.
%   For vectors, NANMEAN(X) is the mean value of the non-NaN
%   elements in X.  For matrices, NANMEAN(X) is a row vector
%   containing the mean value of each column, ignoring NaNs.
%
%   See also NANMEDIAN, NANSTD, NANMIN, NANMAX, NANSUM.

%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.8 $  $Date: 1997/11/29 01:45:53 $

if isempty(x) % Check for empty input.
    y = NaN;
    return
end

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));

if min(size(x))==1,
    count = length(x)-sum(nans);
else
    count = size(x,1)-sum(nans);
end

% Protect against a column of all NaNs
i = find(count==0);
count(i) = ones(size(i));
y = sum(x)./count;
y(i) = i + NaN;

end
%%%%-----------------------------------------------------------------------
function y = nanmedian(x)
%NANMEDIAN NaN protected median value.
%   NANMEDIAN(X) returns the median treating NaNs as missing values.
%   For vectors, NANMEDIAN(X) is the median value of the non-NaN
%   elements in X.  For matrices, NANMEDIAN(X) is a row vector
%   containing the median value of each column, ignoring NaNs.
%
%   See also NANMEAN, NANSTD, NANMIN, NANMAX, NANSUM.

%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.8 $  $Date: 1998/07/21 15:04:39 $

[m,n] = size(x);
x = sort(x); % NaNs are forced to the bottom of each column

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));
if min(size(x))==1,
    n = length(x)-sum(nans);
    if n == 0
        y = NaN;
    else
        if rem(n,2)     % n is odd
            y = x((n+1)/2);
        else            % n is even
            y = (x(n/2) + x(n/2+1))/2;
        end
    end
else
    n = size(x,1)-sum(nans);
    y = zeros(size(n));
    
    % Odd columns
    odd = find(rem(n,2)==1 & n>0);
    idx =(n(odd)+1)/2 + (odd-1)*m;
    y(odd) = x(idx);
    
    % Even columns
    even = find(rem(n,2)==0 & n>0);
    idx1 = n(even)/2 + (even-1)*m;
    idx2 = n(even)/2+1 + (even-1)*m;
    y(even) = (x(idx1)+x(idx2))/2;
    
    % All NaN columns
    i = find(n==0);
    y(i) = i + nan;
end
end
%%%%----------------------------------------------------------------------
function y = nanstd(x)
%NANSTD Standard deviation ignoring NaNs.
%   NANSTD(X) returns the same standard deviation treating NaNs
%   as missing values.
%
%   For vectors, NANSTD(X) is the standard deviation of the
%   non-NaN elements in X.  For matrices, NANSTD(X) is a row
%   vector containing the standard deviation of each column,
%   ignoring NaNs.
%
%   See also NANMEAN, NANMEDIAN, NANMIN, NANMAX, NANSUM.

%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.7 $  $Date: 1997/11/29 01:45:55 $

nans = isnan(x);
i = find(nans);

% Find mean
avg = nanmean(x);

if min(size(x))==1,
    count = length(x)-sum(nans);
    x = x - avg;
else
    count = size(x,1)-sum(nans);
    x = x - avg(ones(size(x,1),1),:);
end

% Replace NaNs with zeros.
x(i) = zeros(size(i));

% Protect against a column of all NaNs
i = find(count==0);
count(i) = ones(size(i));
y = sqrt(sum(x.*x)./max(count-1,1));
y(i) = i + NaN;

end

%%%%----------------------------------------------------------------------
function [pointslist,xselect,yselect] = selectdata(varargin)
% selectdata: graphical selection of data points on a plot using the mouse
% usage: pointslist = selectdata         % uses all default options
% usage: pointslist = selectdata(prop1,val1,prop2,val2,...)
%
% SELECTDATA allows the user to select data points on a given plot
% using the mouse, in a variety of modes. 'Lasso' mode allows the
% selection by a user directed lasso around the points. 'Brush' mode
% selects points as you brush over them with the mouse. 'Rect' mode
% draws a rectangle, selecting any points inside the rectangle.
% 'Closest' mode looks for a aingle mouse click, finding the closest
% point to the mouse.
%
% Returned is a list of the point indexes selected, additionally you
% can specify that the points be deleted from the plot.
%
% Arguments: (input)
%  The input arguments to SELECTDATA are all property/value pairs.
%  (See PARSE_PV_PAIRS for more details.)
%  http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=9082&objectType=FILE
%
%  Property names and character values can be shortened as long as
%  the shortening is unambiguous. Capitalization is ignored.
%
%  Valid Properties: 'Action', 'Axes', 'BrushSize', 'Identify',
%           'Ignore' , 'Pointer', 'SelectionMode'
%
%  'Action' - {'list', 'delete'}
%           'delete' causes the selected points to be deleted from
%           the figure after their selection is final.
%
%           DEFAULT VALUE: 'List'
%
%  'Axes' - denotes the axes from which to select points
%
%           DEFAULT VALUE: gca
%
%  'BrushShape' - {'rect', 'circle'}
%
%           DEFAULT VALUE: 'circle'
%
%           Sets the shape of the brush to be used. Both brush
%           shapes are relative to the figure axes, so a nominally
%           "circular" brush is actually elliptical if the axis
%           units/lengths are unequal.
%
%           Only used when SelectionMode is 'brush'.
%
%  'BrushSize' - Only used when SelectionMode is 'brush'.
%
%           DEFAULT VALUE: 0.05
%
%           The default value will specify a rectangular brush
%           that has dimensions of 5% of the current axes. Note
%           that on a set of square axes, the brush will always
%           look square, even if the axes have very different
%           units.
%
%           0 < brushsize <= 0.25
%
%  'Identify' - {'on', 'off'}
%           Causes the selected points to be temporarily over-plotted
%           with a new filled red marker ('o') as they are selected.
%
%           Points selected with a lasso, or rect may be deselected
%           as the tool is modified. Brush selections are cumulative.
%
%           DEFAULT VALUE: 'on'
%
%  'Ignore' - a data handle, or []
%
%           A list of data handles to be ignored. This allows you to
%           use selectdata on only some sets of points, while others
%           in the same figure are ignored. This is a useful option
%           when you may have plotted some data points but also a
%           curve fit through your data. You can then cause the plotted
%           curve to be ignored by selectdata.
%
%           DEFAULT VALUE: []
%
%  'Pointer' - {'crosshair' | 'fullcrosshair' | 'arrow' | 'ibeam' |
%           'watch' | 'topl' | 'topr' | 'botl' | 'botr' | 'left' |
%           'top' | 'right' | 'bottom' | 'circle' | 'cross' | 'fleur' |
%           'hand' }
%
%           Changes the cursor pointer while selection is active.
%           After selection terminates, the figure pointer is
%           restored to its old setting.
%
%           DEFAULT VALUE: 'crosshair'
%
%  'Return' - {'selected' | 'unselected' }
%
%           Selection of data points, perhaps if used to indicate
%           outliers in data, would normally return that set of
%           points selected. But some users may prefer to see the
%           list of points returned to be those points which were
%           NOT selected.
%
%           DEFAULT VALUE: 'selected'
%
%  'SelectionMode' - {'Lasso', 'Brush', 'Rect', 'Closest'}
%
%           DEFAULT VALUE: 'Lasso'
%
%           If 'Brush' is used, then the brush will be a rectangular
%           brush, with a size as defined by the 'BrushSize' property.
%           The brush will be centered at the current mouse coordinates.
%           Brush size will be a fraction of the axis limits for the
%           current axes. Click inside the axes, then drag the "brush"
%           around. Any points the brush crosses over will be selected.
%
%           If 'Lasso' is chosen, then click inside the axes to define
%           one end of the lasso, then drag with the mouse still down,
%           causing the mouse to define a general curvilinear region.
%           The polygon will close automatically when the mouse button
%           is released. A point is "inside" the lasso if inpolygon
%           identifies it as so. BEWARE of convoluted lassos that
%           intersect themselves.
%
%           If 'Rect' is chosen, then click inside the axes to define
%           one corner of the rect, dragging to specify the opposite
%           corner, just as rbbox would do.
%
%           If 'Closest' was chosen, then a single mouse click in the
%           figure window is used, then that point which is closest in
%           in Euclidean distance (in window units) is chosen. You
%           can move the cursor around (don't release the mouse until
%           you are done) and the currently selected point will be
%           highlighted.
%
%  'Verify' - { 'off' | 'on' }
%
%           If set to 'on', this causes a dialog box to pop up after
%           selection. The user can then acccept the selection, redo
%           it, or cancel out, causing no points to be selected.
%
%           Note, if cancel is chosen from the dialog, and 'return'
%           was specify to return those points 'unselected', then ALL
%           the points will actually be returned.
%
%           DEFAULT VALUE: 'off'
%
%
% Note: other properties are available for use, but I've chosen to
% leave them semi-hidden because they don't seem terribly useful to
% most users. These properties allow you to specify the colors of the
% selection tool itself, the colors of the selected point markers,
% the transparency of the selection tool, the marker itself, etc.
% Default values for these parameters are:
%
%           FlagMarker    = 'o'
%           FlagColor     = 'r'
%           Fill          = 'on'
%           FillColor     = 'y'
%           FillEdgeColor = 'b'
%           FillTrans     = 0.5
%           MaxBrush      = 0.25
%           RemoveTool    = 'on'
%           RemoveFlagged = 'on'
%
% Further documentation on these parameters can be found by editting
% selectdata.m itself.
%
%
% Arguments: (output)
%  pointslist - list of points that were selected. If only one
%           dataset was found, then points list is a simple vector.
%           If multiple sets of points were found in the axes, then
%           points list will be a cell array.
%
%           NOTE: Each set of points is peeled off the stach in the
%           matlab stores it.
%
%  xselect - array (or cell array) containing the x coordinates of
%           those points identified in the selection
%
%  yselect - array (or cell array) containing the y coordinates of
%           those points identified in the selection
%
%
% Example:
%  Plot a set of points, then select some with the mouse
%  using a rectangular brush, return the indices of those
%  points selected. (Note: user interaction on the plot
%  will be necessary.)
%
%    x = 0:.1:1;
%    y = x.^2;
%    plot(x,y,'o')
%    pl = selectdata('selectionmode','brush');
%
% Example:
%  Select a single point with the mouse, then delete that
%  selected point from the plot.
%
%    pl = selectdata('selectionmode','closest','action','delete');
%
% Example:
%  Select some points using a rect(rbbox) tool, also return
%  the (x,y) coordinates from multiple curves plotted. Use
%  shortened versions of the properties and values.
%
%    plot(rand(5,2),rand(5,2),'o')
%    [pl,xs,ys] = selectdata('sel','r');
%
% Example:
%  Plot a curve and some data points on one plot, select
%  some points from the data plotted, but ignore the
%  smooth curve, even if the lasso passes over it.
%    x = 0:.01:1;
%    y = exp(x);
%    ynoisy = y + randn(size(y))/2;
%    h1 = plot(x,y,'-');
%    hold on
%    h2 = plot(x,ynoisy,'o');
%    [pl,xs,ys] = selectdata('sel','lasso','ignore',h1);
%
% See also:
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 2.0
% Release date: 2/7/07

% defaults for the parameters
params.Axes = gca;
params.SelectionMode = 'lasso';
params.Action = 'list';
params.BrushShape = 'circle';
params.BrushSize = .05;
params.Identify = 'on';
params.Ignore = [];
params.Pointer = 'cross';
params.Return = 'selected';
params.Verify = 'off';

% Undocumented options, also unchecked for validity.
% These parameters control the marker used to identify
% those points which are currently selected.
% FlagMarker must be a valid plot marker, FlagColor
% must be a valid color spec. The default values are...
params.FlagMarker = 'o';
params.FlagColor = 'r';

% More (unchecked) options that are yours to fiddle with
% (or not.) These control the fill color to be applied
% to the interior of the lasso, rect, and brush. Also
% controlled are the degree of transparency to be applied.
params.Fill = 'on';
params.FillColor = 'y';
params.FillEdgeColor = 'b';
params.FillTrans = 0.5; % must be in the interval [0,1]

% The maximum relative brushsize allowed is also
% controlled here, just in case I ever wanted to allow
% the brush to be a bit larger.
params.MaxBrush = 0.25;

% After selection has been accomplished, the flagged points
% and the selection tool itself are normally deleted. But
% in some circumstances it may be useful to not delete
% those objects. 'off' will cause these objects to remain.
params.RemoveTool = 'on';
params.RemoveFlagged = 'on';

% save this to restore later
axissize = axis;

% process any Property/value pairs.
if nargin>0
    params = parse_pv_pairs(params,varargin);
end

% check the supplied parameters for validity
params = check_params(params);

% bring the focus to the figure containing the
% designated axes
fighandle = get(params.Axes,'parent');
figure(fighandle)

% get the current figure pointer, so we can
% restore it later on
oldpointer = get(fighandle,'Pointer');

% extract xdata and ydata from the specified axes
% get the children of the axes
hc = findobj(params.Axes,'type','line');


% are any of the data handles to be ignored?
if ~isempty(params.Ignore)
    hc = setdiff(hc,params.Ignore);
end

% strip out xdata and ydata
xdata = get(hc,'xdata');
ydata = get(hc,'ydata');

% if we must highlight the points as they are
% selected, then for efficiency we need to know
% how many we may expect.
flaghandle = [];
if ~iscell(xdata)
    xdata = xdata(:);
    ydata = ydata(:);
    % total number of data points
    npoints = length(xdata);
else
    for i = 1:length(xdata)
        xdata{i} = xdata{i}(:);
        ydata{i} = ydata{i}(:);
    end
    % total number of data points
    npoints = cellfun('length',xdata);
end

% set up a while loop in case we need to verify
% satisfaction
selectionflag = true;

while selectionflag
    % and the total number currently selected
    nsel = 0;
    
    % what SelectionMode was specified
    switch params.SelectionMode
        case 'closest'
            % Find the single closest point (in Euclidean distance)
            % to the mouse click
            
            % set the figure pointer
            if ~isempty(params.Pointer)
                set(fighandle,'Pointer',params.Pointer)
            end
            
            % mouse click?
            waitforbuttonpress;
            
            % Button Motion and Button Up
            set(fighandle,'WindowButtonMotionFcn',@CPmotion);
            set(fighandle,'WindowButtonUpFcn',@selectdone);
            
            % dx, dy to scale the distance
            dx = (axissize(2) - axissize(1));
            dy = (axissize(4) - axissize(3));
            
            % current closest point is
            cc = get(gca,'CurrentPoint');
            xy = cc(1,1:2);
            
            % what point was closest?
            [pointslist,xselect,yselect] = closestpoint(xy,xdata,ydata,dx,dy);
            nsel = 1;
            
            % identify any points?
            if strcmp(params.Identify,'on')
                flagpoints
            end
            
            % selecthandle is not needed for this mode op operation
            selecthandle = [];
            
            % wait until the mouse button was released
            uiwait
            
            % ....
            
            % resume.
            
            % all we need to do here is restore the figure pointer
            % if we changed it before
            if ~isempty(params.Pointer)
                set(fighandle,'Pointer',oldpointer)
            end
            
        case 'rect'
            % Selection rect as a polygon
            
            % mouse click?
            waitforbuttonpress;
            
            % set the figure pointer
            if ~isempty(params.Pointer)
                set(fighandle,'Pointer',params.Pointer)
            end
            
            % button down detected
            cc = get(gca,'CurrentPoint');
            rectxy1 = cc(1,1:2);
            rectxy2 = rectxy1 + eps(rectxy1);
            
            % make a polygon of the box, initially of nil area
            xv = [rectxy1(1), rectxy2(1), rectxy2(1), rectxy1(1), rectxy1(1)];
            yv = [rectxy1(2), rectxy1(2), rectxy2(2), rectxy2(2), rectxy1(2)];
            
            % no points should been selected
            [pointslist,xselect,yselect,nsel] = testpoly(xv,yv,xdata,ydata);
            
            % and plot it, filled or not
            hold on
            if strcmp(params.Fill,'on')
                % filled
                selecthandle = fill(xv,yv,params.FillColor);
                set(selecthandle,'facealpha',params.FillTrans, ...
                    'linestyle','--','edgecolor',params.FillEdgeColor)
            else
                % unfilled
                selecthandle = plot(xv,yv,'r:');
            end
            
            % we can undo the hold now
            hold off
            
            % Button Motion and Button Up
            set(fighandle,'WindowButtonMotionFcn',@rectmotion);
            set(fighandle,'WindowButtonUpFcn',@selectdone);
            
            % wait until the selection is done
            uiwait
            
            % ....
            
            % resume.
            
            % The rect already is a polygon, stored in (xv,yv)
            
        case 'lasso'
            % Selection lasso as a polygon
            
            % mouse click?
            waitforbuttonpress;
            
            % set the figure pointer
            if ~isempty(params.Pointer)
                set(fighandle,'Pointer',params.Pointer)
            end
            
            % button down detected
            cc = get(gca,'CurrentPoint');
            xlasso = cc(1,1);
            ylasso = cc(1,2);
            
            % form the polygon
            xv = xlasso;
            yv = ylasso;
            
            % and plot it, filled or not
            hold on
            if strcmp(params.Fill,'on')
                % filled
                selecthandle = fill(xv,yv,params.FillColor);
                set(selecthandle,'facealpha',params.FillTrans, ...
                    'linestyle','--','edgecolor',params.FillEdgeColor)
            else
                % unfilled
                selecthandle = plot(xv,yv,'r:');
            end
            
            % we can undo the hold now
            hold off
            
            % Button Motion and Button Up
            set(fighandle,'WindowButtonMotionFcn',@lassomotion);
            set(fighandle,'WindowButtonUpFcn',@selectdone);
            
            % wait until the selection is done
            uiwait
            
            % ....
            
            % resume.
            
            % The lasso already is a polygon, stored in (xv,yv)
            
        case 'brush'
            % paint over the data, with a rectangular brush
            
            % mouse click?
            waitforbuttonpress;
            
            % set the figure pointer
            if ~isempty(params.Pointer)
                set(fighandle,'Pointer',params.Pointer)
            end
            
            % button down detected
            bc = get(gca,'CurrentPoint');
            brushcenter = bc(1,1:2);
            
            % dx, dy for the brush
            bdx = params.BrushSize*(axissize(2) - axissize(1));
            bdy = params.BrushSize*(axissize(4) - axissize(3));
            
            if strcmpi(params.BrushShape,'rect')
                % make the brush polygon as a fixed size rectangle
                % that we can slide around
                xv = brushcenter(1) + [-1, 1, 1, -1, -1]*bdx/2;
                yv = brushcenter(2) + [-1, -1, 1, 1, -1]*bdy/2;
            else
                % a circle was specified
                theta = linspace(0,2*pi,100);
                xv = cos(theta)*bdx/2 + brushcenter(1);
                yv = sin(theta)*bdy/2 + brushcenter(2);
            end
            
            % draw the initial brush polygon, filled or not
            hold on
            if strcmp(params.Fill,'on')
                % filled
                selecthandle = fill(xv,yv,params.FillColor);
                set(selecthandle,'facealpha',params.FillTrans, ...
                    'linestyle','-','edgecolor',params.FillEdgeColor)
            else
                % unfilled
                selecthandle = plot(xv,yv,'r:');
            end
            hold off
            
            % have any points been selected?
            [pointslist,xselect,yselect,nsel] = testpoly(xv,yv,xdata,ydata);
            
            % identify any points?
            if strcmp(params.Identify,'on') && (nsel>0)
                flagpoints
            end
            
            % Button Motion and Button Up
            set(fighandle,'WindowButtonMotionFcn',@brushmotion);
            set(fighandle,'WindowButtonUpFcn',@selectdone);
            
            % wait until the selection is done
            uiwait
            
            % ....
            
            % resume.
            
    end
    
    % verify?
    if strcmpi(params.Verify,'on')
        % pop up a dialog
        ButtonName = questdlg( ...
            'Are you satisfied with the points selected?','???', ...
            'Yes','Redo Selection','Cancel Selection','Yes');
        
        switch ButtonName
            case 'Yes'
                % we can drop through
                selectionflag = false;
                
            case 'Cancel Selection'
                % drop out, with nothing selected
                if ~iscell(xdata)
                    pointslist = [];
                    xselect = [];
                    yselect = [];
                else
                    for i = 1:numel(xdata);
                        pointslist{i} = [];
                        xselect{i} = [];
                        yselect{i} = [];
                    end
                end
                
                % we can drop through
                selectionflag = false;
                
            case 'Redo Selection'
                % or try again. The while loop will cycle
                % until happy or canceled
                
        end
        
    else
        % no verification was requested, so we want to
        % terminate the while loop after only one pass through.
        selectionflag = false;
    end
    
end

% pointslist and xselect, yselect are already complete.
% Do we delete the selected points?
if strcmpi(params.Action,'delete')
    if ~iscell(xdata)
        % only one set, so xdata and ydata are not cell arrays
        
        % Which points from the data fall in the selection polygon?
        xdata(pointslist) = [];
        ydata(pointslist) = [];
        
        % drop those points from the plot
        set(hc,'xdata',xdata,'ydata',ydata)
    else
        % it was a cell array, so there were multiple sets.
        for i = 1:numel(xdata);
            
            xdata{i}(pointslist{i}) = [];
            ydata{i}(pointslist{i}) = [];
            
            % drop those points from the plot
            set(hc(i),'xdata',xdata{i},'ydata',ydata{i})
        end
    end
end

% was 'return' set to be the selected list or the unselected one?
% Do nothing if 'selected', we are already done.
if strcmpi(params.Return,'unselected')
    if ~iscell(xdata)
        % only one set, so xdata and ydata are not cell arrays
        pointslist = setdiff((1:npoints)',pointslist);
        xselect = xdata(pointslist);
        yselect = ydata(pointslist);
    else
        % it was a cell array, so there were multiple sets.
        for i = 1:numel(xdata);
            pointslist{i} = setdiff((1:npoints(i))',pointslist{i});
            
            xselect{i} = xdata{i}(pointslist{i});
            yselect{i} = ydata{i}(pointslist{i});
        end
    end
end



% =====================================================
%     begin nested functions
% =====================================================
    function brushmotion(src,evnt) %#ok
        % nested function for motion of the brush
        
        % get the new mouse position
        mousenew = get(params.Axes,'CurrentPoint');
        mousenew = mousenew(1,1:2);
        
        % make sure the axes are fixed
        axis(axissize)
        
        % how far did it move
        brushoffset = mousenew - brushcenter;
        brushcenter = mousenew;
        
        xv = xv + brushoffset(1);
        yv = yv + brushoffset(2);
        
        % did we brush over any new points
        [pl,xselect,yselect,nsel] = testpoly(xv,yv,xdata,ydata);
        
        % for a brush, we need to append any selected points to
        % the already selected list
        if ~iscell(xdata)
            % only one set, so xdata and ydata are not cell arrays
            pointslist = union(pointslist,pl);
            xselect = xdata(pointslist);
            yselect = ydata(pointslist);
            nsel = length(pointslist);
        else
            % it was a cell array, so there were multiple sets.
            for j = 1:numel(pointslist);
                pointslist{j} = union(pointslist{j},pl{j});
                pointslist{j} = pointslist{j}(:);
                
                if ~isempty(pointslist{j})
                    xselect{j} = xdata{j}(pointslist{j});
                    yselect{j} = ydata{j}(pointslist{j});
                end
            end
            
            % total of points selected
            nsel = sum(cellfun('length',pointslist));
        end
        
        % identify any points?
        if strcmp(params.Identify,'on')
            flagpoints
        end
        
        % replot the brush in its new position
        set(selecthandle,'xdata',xv,'ydata',yv)
        
    end
% =====================================================

% =====================================================
    function CPmotion(src,evnt) %#ok
        % nested function to select the closest point
        
        % get the new mouse position
        mousenew = get(params.Axes,'CurrentPoint');
        xy = mousenew(1,1:2);
        
        % make sure the axes stay fixed
        axis(axissize)
        
        % what point was closest?
        [pointslist,xselect,yselect] = closestpoint(xy,xdata,ydata,dx,dy);
        nsel = 1;
        
        % identify any points?
        if strcmp(params.Identify,'on')
            flagpoints
        end
        
    end
% =====================================================

% =====================================================
    function rectmotion(src,evnt) %#ok
        % nested function for expansion or contraction of the rect
        
        % get the new mouse position
        mousenew = get(params.Axes,'CurrentPoint');
        rectxy2 = mousenew(1,1:2);
        
        % make sure the axes are fixed
        axis(axissize)
        
        % update the rect polygon of the box, changing the second corner
        xv = [rectxy1(1), rectxy2(1), rectxy2(1), rectxy1(1), rectxy1(1)];
        yv = [rectxy1(2), rectxy1(2), rectxy2(2), rectxy2(2), rectxy1(2)];
        
        % did we brush over any new points
        [pointslist,xselect,yselect,nsel] = testpoly(xv,yv,xdata,ydata);
        
        % identify any points?
        if strcmp(params.Identify,'on')
            flagpoints
        end
        
        % replot the rect in its new position
        set(selecthandle,'xdata',xv,'ydata',yv)
        
    end
% =====================================================

% =====================================================
    function lassomotion(src,evnt) %#ok
        % nested function for expansion of the lasso
        
        % get the new mouse position
        mousenew = get(params.Axes,'CurrentPoint');
        mousenew = mousenew(1,1:2);
        
        % append the new point on the end of the last lasso
        xlasso = [xlasso,mousenew(1,1)];
        ylasso = [ylasso,mousenew(1,2)];
        
        % and close it to form the polygon
        xv = [xlasso,xlasso(1)];
        yv = [ylasso,ylasso(1)];
        
        % replot the newly extended lasso
        set(selecthandle,'xdata',xv,'ydata',yv)
        
        % did we enclose any new points?
        [pointslist,xselect,yselect,nsel] = testpoly(xv,yv,xdata,ydata);
        
        % identify any points?
        if strcmp(params.Identify,'on')
            flagpoints
        end
        
        % make sure the axes are maintained in size
        axis(axissize)
        
    end
% =====================================================

% =====================================================
    function selectdone(src,evnt) %#ok
        % nested function for mouse up
        
        % do we remove the selection tool?
        if strcmpi(params.RemoveTool,'on')
            % delete the selection object from the plot
            delete(selecthandle)
            selecthandle = [];
        end
        
        % do we remove the selection tool?
        if strcmpi(params.RemoveFlagged,'on')
            % also remove the flagged/plotted points
            if ~isempty(flaghandle)
                delete(flaghandle)
                flaghandle = [];
            end
        end
        
        % cancel the WindowButtonFcn's that we had set
        set(fighandle,'WindowButtonMotionFcn',[]);
        set(fighandle,'WindowButtonUpFcn',[]);
        
        % restore the figure pointer to its original setting
        if ~isempty(params.Pointer)
            set(fighandle,'Pointer',oldpointer)
        end
        
        % and resume execution, back in the mainline
        uiresume
    end
% =====================================================

% =====================================================
    function flagpoints
        % nested function for flagging the selected points
        
        % Are these the first points flagged? If so,
        % we need to plot them and set the marker, etc.
        if isempty(flaghandle) && (nsel > 0)
            % hold the figure, so we can add the flagged points
            hold on
            
            if ~iscell(xselect)
                flaghandle = plot(xselect,yselect,params.FlagMarker);
                set(flaghandle,'Color',params.FlagColor,'MarkerFaceColor',params.FlagColor)
            else
                flaghandle = plot(vertcat(xselect{:}),vertcat(yselect{:}),params.FlagMarker);
                set(flaghandle,'Color',params.FlagColor,'MarkerFaceColor',params.FlagColor)
            end
            
            % now release the hold
            hold off
        elseif ~isempty(flaghandle)
            % otherwise, we just need to update xdata and ydata
            
            if nsel == 0
                set(flaghandle,'xdata',[],'ydata',[]);
                
            elseif ~iscell(xselect)
                set(flaghandle,'xdata',xselect,'ydata',yselect);
            else
                set(flaghandle,'xdata',vertcat(xselect{:}),'ydata',vertcat(yselect{:}));
            end
        end
    end
% =====================================================

end % mainline end

% ================================================
%               end main function
% ================================================

% ================================================
%                  subfunctions
% ================================================
function [pl,xsel,ysel,nsel] = testpoly(xv,yv,xdata,ydata)
% checks which points are inside the given polygon

% was there more than one set of points found in the plot?
if ~iscell(xdata)
    % only one set, so xdata and ydata are not cell arrays
    
    % Which points from the data fall in the selection polygon?
    pl = find(inpolygon(xdata,ydata,xv,yv));
    nsel = length(pl);
    
    xsel = xdata(pl);
    ysel = ydata(pl);
else
    % it was a cell array, so there were multiple sets.
    pl = cell(size(xdata));
    xsel = pl;
    ysel = pl;
    nsel = 0;
    for i = 1:numel(xdata);
        pl{i} = find(inpolygon(xdata{i},ydata{i},xv,yv));
        nsel = nsel + length(pl{i});
        
        if ~isempty(pl{i})
            xsel{i} = xdata{i}(pl{i});
            ysel{i} = ydata{i}(pl{i});
        end
        
    end
end

end % subfunction end

% ================================================
%                  subfunction
% ================================================
function [pointslist,xselect,yselect] = closestpoint(xy,xdata,ydata,dx,dy)
% find the single closest point to xy, in scaled units
if ~iscell(xdata)
    % just one set of points to consider
    D = sqrt(((xdata - xy(1))/dx).^2 + ((ydata - xy(2))/dy).^2);
    [junk,pointslist] = min(D(:));
    xselect = xdata(pointslist);
    yselect = ydata(pointslist);
else
    % there is more than one set of points
    Dmin = inf;
    pointslist = cell(size(xdata));
    for i = 1:numel(xdata);
        D = sqrt(((xdata{i} - xy(1))/dx).^2 + ((ydata{i} - xy(2))/dy).^2);
        [mind,ind] = min(D(:));
        
        if mind < Dmin
            % searching for the closest point
            Dmin = mind;
            
            pointslist = cell(size(xdata));
            xselect = cell(size(xdata));
            yselect = cell(size(xdata));
            
            pointslist{i} = ind;
            xselect{i} = xdata{i}(ind);
            yselect{i} = ydata{i}(ind);
        end
    end
end

end % subfunction end

% ================================================
%                  subfunction
% ================================================

% ============================================
% subfunction - check_params
% ============================================
function par = check_params(par)
% check the parameters for acceptability
%
% Defaults
%  Axes = gca;
%  SelectionMode = 'lasso';
%  Action = 'list';
%  BrushSize = .05;

% Axes == gca by default
if isempty(par.Axes)
    par.Axes = gca;
else
    if ~ishandle(par.Axes)
        error 'Axes must be the handle to a valid set of axes.'
    end
end

% SelectionMode == 'brush' by default
if isempty(par.SelectionMode)
    par.SelectionMode = 'brush';
else
    valid = {'rect', 'brush', 'lasso', 'closest'};
    if ~ischar(par.SelectionMode)
        error 'Invalid Style: Must be character'
    end
    
    ind = strmatch(lower(par.SelectionMode),valid);
    if isempty(ind) || (length(ind)>1)
        error(['Invalid SelectionMode: ',par.SelectionMode])
    end
    par.SelectionMode = valid{ind};
end

% BrushShape == 'circle' by default
if isempty(par.BrushShape)
    par.BrushShape = 'circle';
else
    valid = {'rect', 'circle'};
    if ~ischar(par.BrushShape)
        error 'Invalid Style: Must be character'
    end
    
    ind = strmatch(lower(par.BrushShape),valid);
    if isempty(ind) || (length(ind)>1)
        error(['Invalid SelectionMode: ',par.BrushShape])
    end
    par.BrushShape = valid{ind};
end

% Action == 'list' by default
if isempty(par.Action)
    par.Action = 'list';
else
    valid = {'list', 'delete'};
    if ~ischar(par.Action)
        error 'Invalid Action: Must be character'
    end
    
    ind = strmatch(lower(par.Action),valid);
    if isempty(ind) || (length(ind)>1)
        error(['Invalid Action: ',par.Action])
    end
    par.Action = valid{ind};
end

% Pointer == 'crosshair' by default, but
% if empty, will not change the pointer.
if ~isempty(par.Pointer)
    valid = {'crosshair', 'fullcrosshair', 'arrow', 'ibeam', ...
        'watch', 'topl', 'topr', 'botl', 'botr', 'left', 'top', ...
        'right', 'bottom', 'circle', 'cross', 'fleur', ...
        'custom', 'hand'};
    
    if ~ischar(par.Pointer)
        error 'Invalid Pointer: Must be character'
    end
    
    ind = strmatch(lower(par.Pointer),valid,'exact');
    if isempty(ind)
        ind = strmatch(lower(par.Pointer),valid);
        if isempty(ind) || (length(ind)>1)
            error(['Invalid Pointer: ',par.Pointer])
        end
    end
    par.Pointer = valid{ind};
end

% Identify == 'on' by default
if isempty(par.Identify)
    par.Identify = 'on';
else
    valid = {'on', 'off'};
    if ~ischar(par.Identify)
        error 'Value for Identify is invalid: Must be character'
    end
    
    ind = strmatch(lower(par.Identify),valid);
    if isempty(ind) || (length(ind)>1)
        error(['Invalid Action: ',par.Identify])
    end
    par.Identify = valid{ind};
end

% Return == 'selected' by default
if isempty(par.Return)
    par.Return = 'selected';
else
    valid = {'selected', 'unselected'};
    if ~ischar(par.Return)
        error 'Value for Return is invalid: Must be character'
    end
    
    ind = strmatch(lower(par.Return),valid);
    if isempty(ind) || (length(ind)>1)
        error(['Invalid Action: ',par.Return])
    end
    par.Return = valid{ind};
end

% Verify == 'off' by default
if isempty(par.Verify)
    par.Verify = 'off';
else
    valid = {'on', 'off'};
    if ~ischar(par.Verify)
        error 'Value for Verify is invalid: Must be character'
    end
    
    ind = strmatch(lower(par.Verify),valid);
    if isempty(ind) || (length(ind)>1)
        error(['Invalid Action: ',par.Verify])
    end
    par.Verify = valid{ind};
end

% BrushSize == 0.05 by default
if isempty(par.BrushSize)
    par.BrushSize = 0.05;
else
    if (length(par.BrushSize)>1) || (par.BrushSize<=0) || (par.BrushSize>par.MaxBrush)
        error 'Brushsize must be scalar, and 0 < BrushSize <= 0.25'
    end
end

% Ignore == [] by default
if ~isempty(par.Ignore) && any(~ishandle(par.Ignore))
    error 'Ignore must be empty, or a data handle'
end

end % check_params


% ============================================
% Included subfunction - parse_pv_pairs
% ============================================
function params=parse_pv_pairs(params,pv_pairs)
% parse_pv_pairs: parses sets of property value pairs, allows defaults
% usage: params=parse_pv_pairs(default_params,pv_pairs)
%
% arguments: (input)
%  default_params - structure, with one field for every potential
%             property/value pair. Each field will contain the default
%             value for that property. If no default is supplied for a
%             given property, then that field must be empty.
%
%  pv_array - cell array of property/value pairs.
%             Case is ignored when comparing properties to the list
%             of field names. Also, any unambiguous shortening of a
%             field/property name is allowed.
%
% arguments: (output)
%  params   - parameter struct that reflects any updated property/value
%             pairs in the pv_array.
%
% Example usage:
% First, set default values for the parameters. Assume we
% have four parameters that we wish to use optionally in
% the function examplefun.
%
%  - 'viscosity', which will have a default value of 1
%  - 'volume', which will default to 1
%  - 'pie' - which will have default value 3.141592653589793
%  - 'description' - a text field, left empty by default
%
% The first argument to examplefun is one which will always be
% supplied.
%
%   function examplefun(dummyarg1,varargin)
%   params.Viscosity = 1;
%   params.Volume = 1;
%   params.Pie = 3.141592653589793
%
%   params.Description = '';
%   params=parse_pv_pairs(params,varargin);
%   params
%
% Use examplefun, overriding the defaults for 'pie', 'viscosity'
% and 'description'. The 'volume' parameter is left at its default.
%
%   examplefun(rand(10),'vis',10,'pie',3,'Description','Hello world')
%
% params =
%     Viscosity: 10
%        Volume: 1
%           Pie: 3
%   Description: 'Hello world'
%
% Note that capitalization was ignored, and the property 'viscosity'
% was truncated as supplied. Also note that the order the pairs were
% supplied was arbitrary.

npv = length(pv_pairs);
n = npv/2;

if n~=floor(n)
    error 'Property/value pairs must come in PAIRS.'
end
if n<=0
    % just return the defaults
    return
end

if ~isstruct(params)
    error 'No structure for defaults was supplied'
end

% there was at least one pv pair. process any supplied
propnames = fieldnames(params);
lpropnames = lower(propnames);
for i=1:n
    p_i = lower(pv_pairs{2*i-1});
    v_i = pv_pairs{2*i};
    
    ind = strmatch(p_i,lpropnames,'exact');
    if isempty(ind)
        ind = find(strncmp(p_i,lpropnames,length(p_i)));
        if isempty(ind)
            error(['No matching property found for: ',pv_pairs{2*i-1}])
        elseif length(ind)>1
            error(['Ambiguous property name: ',pv_pairs{2*i-1}])
        end
    end
    p_i = propnames{ind};
    
    % override the corresponding default in params
    params = setfield(params,p_i,v_i); %#ok
    
end

end % parse_pv_pairs
%%%%%----------------------------------------------------------------------
function [d, id] = getchunks(a, opt)

%GETCHUNKS Get the number of repetitions that occur in consecutive chunks.
%   C = GETCHUNKS(A) returns an array of n elements, where n is the number
%   of consecutive chunks (2 or more repetitions) in A, and each element is
%   the number of repetitions in each chunk. A can be LOGICAL, any
%   numeric vector, or CELL array of strings. It can also be a character
%   array (see below, for its special treatment).
%
%   Example:
%     A = [1 2 2 3 4 4 4 5 6 7 8 8 8 8 9];
%     getchunks(A)
%       ans =
%           2   3   4
%
%   [C, I] = GETCHUNKS(A) also returns the indeces of the beginnings of the
%   chunks.
%
%   GETCHUNKS(A, OPT) when OPT is '-full', includes single-element chunks
%   as well. Default OPT is '-reps' (meaning REPEATING chunks).
%
%   If A is a character array, then it finds words (consecutive
%   non-spaces), returning number of chararcters in each word and the
%   indeces to the beginnings of the words. In this case, OPT has no
%   effect.
%
%   Example:
%     A = 'This is a generic (simple) sentence';
%     [C, I] = getchunks(A)
%       C =
%            4     2     1     7     8     8
%
%       I =
%            1     6     9    11    19    28
%
%   See also HIST, HISTC.

%   Jiro Doke
%   Feb 16, 2006
%

%--------------------------------------------------------------------------
% Error checking
narginchk(1, 2);
if ndims(a) > 2 || min(size(a)) > 1
    error('Input must be a 2-D vector');
end

if nargin < 2
    opt = '-reps';
end

%--------------------------------------------------------------------------
% Process options
switch lower(opt)
    
    % Include single-element chunks
    case '-full'
        fullList = true;
        if ischar(a);
            fprintf('''-full'' option not applicable with CHAR arrays.\n');
        end
        
        % Only find 2 or more repeating blocks
    case '-reps'
        fullList = false;
        
    otherwise
        error('Unknown option. Allowed option: ''-full'' or ''-reps''');
end

%--------------------------------------------------------------------------
% Convert to row vector for STRFIND
a = a(:)';

%--------------------------------------------------------------------------
% Deal with differet classes
switch class(a)
    
    case 'double'
        % Leave as is
        
    case {'logical', 'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'single'}
        % Convert to DOUBLE
        a = double(a);
        
    case 'char'
        % Get non-space locations
        a = ~isspace(a);
        
    case 'cell'
        % Convert cell array of strings into unique numbers
        if all(cellfun('isclass', a, 'char'))
            [a, a, a] = unique(a);
        else
            error('Cell arrays must be array of strings.');
        end
        
    otherwise
        error('Invalid type. Allowed type: CHAR, LOGICAL, NUMERIC, and CELL arrays of strings.');
end

%--------------------------------------------------------------------------
% Character arrays (now LOGICAL) are dealt differently
if islogical(a)
    % Pad the array
    a  = [false, a, false];
    
    % Here's a very convoluted engine
    b  = diff(a);
    id = strfind(b, 1);
    d  = strfind(b, -1) - id;
    
    %--------------------------------------------------------------------------
    % Everything else (numeric arrays) are processed here
else
    % Pad the array
    a                 = [NaN, a, NaN];
    
    % Here's more convoluted code
    b                 = diff(a);
    b1                = b;  % to be used in fullList (below)
    ii                = true(size(b));
    ii(strfind(b, 0)) = false;
    b(ii)             = 1;
    c                 = diff(b);
    id                = strfind(c, -1);
    
    % Get single-element chunks also
    if fullList
        
        % And more convoluted code
        b1(id)          = 0;
        ii2             = find(b1(1:end-1));
        d               = [strfind(c, 1) - id + 1, ones(1, length(ii2))];
        id              = [id,ii2];
        [id,tmp]        = sort(id);
        d               = d(tmp);
        
    else
        d               = strfind(c, 1) - id + 1;
    end
end
end
%%%%%----------------------------------------------------------------------
function varargout = datevecfix(t,varargin)
%DATEVECFIX   Date components with rounded seconds to specified precision.
%
%   SYNTAX:
%                   V = DATEVEC(N);
%                   V = DATEVEC(S);
%                   V = DATEVEC(S,F);
%                   V = DATEVEC(S,P);
%                   V = DATEVEC(S,P,F);
%                   V = DATEVEC(S,F,P);
%                   V = DATEVEC(...,'Precision',K);
%     [Y,MO,D,H,MI,S] = DATEVEC(...);
%
%   INPUT:
%     N - Serial date (see DATENUM).
%     S - String date (see DATESTR).
%     F - Format date (see DATESTR).
%     P - Pivot date  (see DATESTR).
%     K - Seconds precision. Must be a positive integer or zero.
%         DEFAULT: 0 (rounds up to seconds)
%
%   OUTPUT:
%     V  - 6 colums matrix with [Y MO D H MI S] values.
%     Y  - Years.
%     MO - Months.
%     D  - Days.
%     H  - Hours.
%     MI - Minutes.
%     S  - Seconds.
%
%
%   DESCRIPTION:
%     This program is the same as DATEVEC but rounds the seconds up to the
%     specified precision and fixes the problem with the 60 seconds not
%     rounded value.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Option argument name 'Precision' may be a single 'P'.
%
%   EXAMPLE:
%     datevec(   366.9999999999999), % returns: [0 12 31 23 59 60]
%     60-ans(6),                     % returns: 1.1176e-008
%     datevecfix(366.9999999999999), % returns: [1  1  1  0  0  0]
%      0-ans(6),                     % returns: 0
%
%   SEE ALSO:
%     DATEVEC
%     and
%     TLABEL by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   datevecfix.m
%   VERSION: 2.0 (Jun 08, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>)
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Mar 04, 2008)
%   2.0      Added precision optional input. (Jun 08, 2009)

%   DISCLAIMER:
%   datevecfix.m is provided "as is" without warranty of any kind, under
%   the revised BSD license.

%   Copyright (c) 2008-2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Defaults:
K = 0;

% Checks number of inputs:
if     nargin<1
    error('CVARGAS:datevecfix:notEnoughInputs',...
        'At least 1 input is required.')
elseif nargout>7
    error('CVARGAS:datevecfix:tooManyOutputs',...
        'At most 6 outputs are allowed.')
end

% Checks precision:
if ((nargin>2) && (~isempty(varargin{end-1}) && ...
        strncmpi(varargin{end-1},'precision',min(length(varargin{end-1}),9))))
    K = varargin{end};
    if (~isempty(K) && isnumeric(K) && (numel(K)==1))
        K = round(abs(K));
    else
        error('CVARGAS:datevecfix:incorrectPrecisionInput',...
            'Precision must be a single integer.')
    end
    varargin(end-1:end) = [];
end


% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

% Find vector:
x = datevec(t,varargin{:});

% Rounds seconds up to precision:
x(:,6) = round(x(:,6)*10^K);

% Check seconds:
ibad      = (x(:,6)>=60*10^K);
x(ibad,6) = x(ibad,6)-60*10^K;
x(:,6)    = x(:,6)/10^K;
x(ibad,5) = x(ibad,5)+1;
% Check minutes:
ibad = (x(:,5)>=60);
x(ibad,5) = x(ibad,5)-60;
x(ibad,4) = x(ibad,4)+1;
% Check hours:
ibad = (x(:,4)>=24);
x(ibad,4) = x(ibad,4)-24;
x(ibad,3) = x(ibad,3)+1;
% Check days, months and years:
jbad = x(ibad,3)>=28; % only problematic days
if any(jbad)
    y = datevec(datenum(x(ibad(jbad),1:3)));
    x(ibad(jbad),1:3) = y(:,1:3);
end


% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

if nargout>1
    for k = 1:nargout
        varargout{k} = x(:,k);
    end
else
    varargout{1} = x;
end


% [EOF]   datevecfix.m
end
%%%%----------------------------------------------------------------------
function gescatter(filename,lon,lat,c,varargin)
% GESCATTER - create a scatter plot in Google Earth
%
%   GESCATTER(FILENAME,LON,LAT,C) - creates a .kml file that
%   displays colored circles at the locations specified by the
%   vectors LON and LAT similar to ML's builtin function, SCATTER.
%   The color of the circles is scaled relative to the
%   values provided in third input, C.
%
%   OPTIONS AND SYNTAX - Optional inputs are entered as
%   property/value pairs.  Valid optional properties are:
%
%   GESCATTER(...,'colormap','hot') - uses Matlabs 'hot'
%   colormap instead of the default (jet). Also accepts function
%   handles (@hot), or custom colormaps (m x 3 matrices).
%
%   GESCATTER(...,'clims',[low high]) - limit the color
%   values to a specified values (similar to CAXIS on a ML
%   figure). Clims should be supplied as a 2-element array.
%
%   GESCATTER(...,'time',timevector) - assigns a time to each
%   point. The length of the timevector array should be the same
%   as LAT, LON, and C.
%
%   GESCATTER(...,'scale',size) - scales the size of the dots in the
%   Google Earth file.  Default value is 0.4.
%
%   EXAMPLE
%
%   %generate some data
%   x=(0:0.05:6*pi);
%   lon = -122.170087 + cos(x)*0.01;
%   lat = 37.455697 + x*0.001;
%
%   %color the points according to their latitude
%   gescatter('foo.kml',lon,lat,lat)
%
% SEE ALSO scatter

% A. Stevens @ USGS 3/04/2009
% astevens@usgs.gov

%default values
clims=[min(c) max(c)];
cmap=fliplr(jet);
t=[];
scale=0.5;


%parse inputs and do some error-checking
if nargin>0
    [m,n]=size(varargin);
    opts={'clims','time','scale','colormap'};
    for i=1:n;
        indi=strcmpi(varargin{i},opts);
        ind=find(indi==1);
        if isempty(ind)~=1
            switch ind
                case 1
                    clims=varargin{i+1};
                    if numel(clims)~=2
                        error('Clims should be a two-element array.')
                    end
                    
                case 2
                    t=varargin{i+1};
                    if any(isnan(t))
                        error('Time vector should not contain NaNs.')
                    end
                    if ~isnumeric(t)
                        error('Time should be entered in ML Datenum format.')
                    end
                case 3
                    scale=varargin{i+1};
                case 4
                    cmap=varargin{i+1};
                    %check size of numeric colormap input
                    if isa(cmap,'numeric')
                        [m,n]=size(cmap);
                        if n~=3
                            error('Custom colormap must have 3 columns.')
                        end
                        cmap=fliplr(cmap);
                    else
                        %if standard colormap is supplied
                        if isa(cmap,'function_handle')
                            cmap= func2str(cmap);
                        end
                        cmap=fliplr(feval(cmap));
                    end
                    
            end
        end
    end
end


[pathstr,namer] = fileparts(filename);

%get rid on nans
gind=(isfinite(lon) & isfinite(lat) & isfinite(c));
lon=lon(gind);
lat=lat(gind);
c=c(gind);


%figure out the rgb colors of each value
cvals=[-inf;linspace(clims(1),clims(2),...
    length(cmap(:,1))-2)';inf];
[n,bin]=histc(c,cvals);
colors=cmap(bin,:);

%convert to GE's hex format
rgb=cellfun(@(x)(dec2hex(floor(x.*255),2)),...
    num2cell(colors),'uni',0);

%write the GE file
header=['<?xml version="1.0" encoding="UTF-8"?>',...
    '<kml xmlns="http://www.opengis.net/kml/2.2">',...
    '<Document><name>',namer,'</name>'];
footer='</Document></kml>';

h = waitbar(0,'Creating file, Please wait...');
set(h,'name','Creating Google Earth file')

fid = fopen(filename, 'wt');
fprintf(fid, '%s \n',header);

for i=1:length(lon)
    
    %create a style to hide each point in one document
    fprintf(fid,'%s \n','<Style id="folderStyle">');
    fprintf(fid,'%s \n','<ListStyle>');
    fprintf(fid,'%s \n','<listItemType>checkHideChildren</listItemType>');
    fprintf(fid,'%s \n','</ListStyle>');
    fprintf(fid,'%s \n','</Style>');
    
    %define the point style
    fprintf(fid,'%s \n','<Style id="cpoint">');
    fprintf(fid,'%s \n','<IconStyle>');
    fprintf(fid,'%s \n',['<color>ff',[rgb{i,:}],'</color>']);
    fprintf(fid,'%s \n',['<scale>',sprintf('%.1f',scale),'</scale>']);
    fprintf(fid,'%s \n',['<Icon><href>http://maps.google.com/mapfiles/',...
        'kml/shapes/shaded_dot.png</href></Icon>']);
    fprintf(fid,'%s \n','</IconStyle>');
    fprintf(fid,'%s \n','</Style>');
    
    %add the placemark
    fprintf(fid, '%s \n','<Placemark>');
    fprintf(fid,'%s \n','<styleUrl>#cpoint</styleUrl>');
    
    %create a simple description for each point
    fprintf(fid, '%s \n','<description><![CDATA[<table width="200"></table>');
    fprintf(fid, '%s \n',['<h2>Filename: ',namer,'<br>']);
    fprintf(fid, '%s \n',['<h3>Value: ',sprintf('%.1f',c(i)),'<br>']);
    if ~isempty(t)
        fprintf(fid, '%s \n',['Time (GMT): ',datestr(t(i)),'<br>']);
    end
    fprintf(fid, '%s \n',']]></description>');
    
    
    fprintf(fid,'%s \n','<Point>');
    fprintf(fid,'%s','<coordinates>');
    fprintf(fid, ' %.6f, %.6f, %.2f', [lon(i) lat(i) c(i)]);
    fprintf(fid,'%s \n','</coordinates>');
    fprintf(fid,'%s \n','</Point>');
    
    if ~isempty(t)
        fprintf(fid,'%s \n','<TimeSpan>');
        fprintf(fid,'%s \n',['<begin>',datestr(t(1),29),...
            'T',datestr(t(1),13),'Z</begin>']);
        fprintf(fid,'%s \n',['<end>',datestr(t(end),29),...
            'T',datestr(t(end),13),'Z</end>']);
        fprintf(fid,'%s \n','</TimeSpan>');
    end
    
    
    fprintf(fid, '%s \n','</Placemark>');
    
    waitbar(i/length(lon),h,sprintf('%d%% complete...',...
        round((i/length(lon))*100)));
    
end

fprintf(fid, '%s \n','<styleUrl>#folderStyle</styleUrl>');
fprintf(fid, '%s \n',footer);

close(h);
fclose(fid);
end
%%%%----------------------------------------------------------------------
function bindata = readBin(varargin)
% READBIN - Read raw acoustic data from .BIN file
%
%   READBIN reads binary files with the extention .BIN
%   containing raw waveform acoustic data from a single-
%   beam echosounder recorded with HYPACK. A .RAW file with
%   the same name must be present in the same directory
%   as the .BIN file. The code is currently only configured
%   for a single-frequency sonar system.
%
%   INPUTS: Two optional input arguements are available
%       filename - If a string is among the input arguments,
%                  the program will use this as the file to
%                   process.  If no filename is provided,
%                   the user will be prompted to select a file.
%       plotflag - 0 (default) or 1. If plotflag is set to 1,
%                  a simple plot will be made.
%
%   OUTPUT: The data are returned in a structure with
%   the following fileds:
%       'filename' -  name of data file
%       'vtime'    -  Hypack time tag (millisecs since midnight)
%       'range'    -  range from transducer (m)
%       'vals'     -  backscatter values (unknown units)
%       'mtime'    -  Matlab datenum
%       'lat'      -  latitude (dec. degrees, WGS84)
%       'lon'      -  longitude (dec. degrees, WGS84)
%       'x'        -  projected x-coordinate
%       'y'        -  projected y-coordinate
%       'tide'     -  tide correction (m)
%       'zraw'     -  digitized depth, uncorrected
%       'zc'       -  digitized depth, corrected
%
%   EXAMPLES:
%
%   data=readBIN; %user selects file, no plot made
%
%   data=readBIN('d:\data\foo.bin'); %user supplies file, no plot
%
%   data=readBIN(1); %user selects file, plot is made
%
% SEE ALSO readRAW transectViewer

% Andrew Stevens
% astevens@usgs.gov
% 8/12/2009


%check inputs
narginchk(0,2);

if nargin>0
    for i=1:length(varargin)
        
        if ischar(varargin{i});
            
            fpath=varargin{i};
            [dname,fname] = fileparts(fpath);
            if (exist([dname,filesep,fname,'.bin'],'file')==0 || ...
                    exist([dname,filesep,fname,'.raw'],'file')==0)
                error('.BIN and/or .RAW files not found. Try again');
            end
            
        end
        
    end
end


% h=waitbar(0);

%
% set(h,'name','Reading .BIN file');

fid=fopen([dname,filesep,fname,'.bin'],'r');
fseek(fid,0,'eof');
numbytes=ftell(fid);
frewind(fid)

pos=ftell(fid);
numread=1;
while pos<numbytes
    
    dummy=fread(fid,1,'uint32'); %#ok
    vers=fread(fid,1,'uint8'); %#ok
    dep_no=fread(fid,1,'uint8'); %#ok
    timer=fread(fid,1,'uint32')./1000;
    draft=fread(fid,1,'ushort'); %#ok
    unit=fread(fid,1,'char'); %#ok
    sw(numread)=fread(fid,1,'ushort'); %#ok
    eos(numread)=fread(fid,1,'ushort'); %#ok
    num_samp(numread)=fread(fid,1,'ushort');
    if isempty(num_samp(numread))
        break
    end
    fseek(fid,16,'cof');
    samp=fread(fid,num_samp(numread),'ushort');
    
    if numread==1
        bytes_per_ping=35+(num_samp(numread)*2);
        num_pings=numbytes/bytes_per_ping;
        
        bindata=struct('filename',fname,'vtime',zeros(1,num_pings),...
            'range',zeros(num_samp(numread),num_pings),...
            'vals',zeros(num_samp(numread),num_pings));
        
    end
    
    if isempty(samp), 
        eos(numread)=eos(numread-1);
        sw(numread)=sw(numread-1);
        num_samp(numread)=num_samp(numread-1);
        break 
    end
    
    bindata.range(:,numread)=linspace(eos(numread)-sw(numread),...
        eos(numread),num_samp(numread));
    
    
    
    bindata.vtime(numread)=timer;
    bindata.vals(:,numread)=samp;
    
    
    numread=numread+1;
    pos=ftell(fid);
    
    %     waitbar(pos/numbytes,h,sprintf('%d%% complete...',...
    %         round((pos/numbytes)*100)));
    
end
fclose(fid);


if (numel(unique(sw))>1 || numel(unique(eos))>1);
    %     set(h,'name','Processing Bin File.');
    rangei=linspace(min(eos)-max(sw),...
        max(eos),max(num_samp));
    
    for i = 1:numread-1;
        bindata.vals(:,i)=interp1(bindata.range(:,i),...
            bindata.vals(:,i),rangei);
        
        %         waitbar(i/(numread-1),h,sprintf('%d%% complete...',...
        %                 round((i/(numread-1))*100)));
    end
    bindata.range=rangei;
else
    bindata.range=bindata.range(:,1);
end
% close(h);




end



%%%%-----------------------------------------------------------------------
function showOffline(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);


ind=find(isfinite(gdata.bdata.zc));

if length(gdata.hdr.linex)==2
    
    gdata.offFig=figure;
    l1=line([min(gdata.bdata.distance) max(gdata.bdata.distance)],[0 0]);
    set(l1,'color','r')
    hold on
    gdata.offPlot=plot(gdata.bdata.distance(ind),gdata.bdata.offline(ind));
    
    xlabel('Distance (m)','fontsize',14)
    ylabel('Offline (m)','fontsize',14)
    % set(gca,'xlim',[min(gdata.hdr.liney) max(gdata.hdr.liney)])
    box on
    
    guidata(hfig,gdata);
    
end
end

%%%raw data editing tools-------------------------------------------------

function cleanProfileRaw(hfig,eventdata,handles)%#ok

gdata=guidata(hfig);
set(gdata.push5,'backgroundcolor','g');

if ~isfield(gdata,'eraw');
    gdata.eraw=gdata.bdata.zraw;
    gdata.ezc=gdata.bdata.zc;
end
if gdata.numedits==0
    set(gdata.menu12,'visible','on');
    set(gdata.menu13,'visible','on');
end
gdata.numedits=gdata.numedits+1;


if gdata.pan==1;
    set(gdata.toggle1,'value',0);
    pan off
end

gdata.l1=findobj(gca,'type','line');
[pl,xs,ys] = selectdata('sel','lasso','ignore',gdata.rim);

gdata.bdata.zc(pl)=NaN;
gdata.bdata.zraw(pl)=NaN;

set(gdata.l1,'ydata',gdata.bdata.zraw)

ind=find(isfinite(gdata.bdata.zc));

gdata.xlimo=[min(gdata.bdata.distance(ind)),...
    max(gdata.bdata.distance(ind))];
gdata.ylimo=[min(gdata.bdata.zc(ind)),...
    max(gdata.bdata.zc(ind))];

gdata.edits{gdata.numedits}=pl;

waitfor(pl)
guidata(hfig,gdata);
setFocus(hfig)


if ishandle(gdata.offFig)
    figure(gdata.offFig);
    
    set(gdata.offPlot,'xdata',...
        gdata.bdata.distance(isfinite(gdata.bdata.zc)),...
        'ydata',gdata.bdata.offline(isfinite(gdata.bdata.zc)))
end

set(gdata.push5,'backgroundcolor',[0.9255    0.9137    0.8471]);

setFocus(hfig)

end

%%%%----------------------------------------------------------------------
function clearEdits(hfig,evnt) %#ok
gdata=guidata(hfig);

edits=cell2mat(gdata.edits(:));
gdata.bdata.zraw(edits)=...
    gdata.eraw(edits);
gdata.bdata.zc(edits)=...
    gdata.ezc(edits);

gdata.edits=[];
gdata.numedits=0;

if gdata.rawflag==0
    set(gdata.l1,'ydata',gdata.bdata.zc);
else
    set(gdata.l1,'ydata',gdata.bdata.zraw)
end


set(gdata.menu12,'visible','off')
set(gdata.menu13,'visible','off')

guidata(hfig,gdata);


end
%%%%----------------------------------------------------------------------
function digitizebtm(hfig,evnt) %#ok
gdata=guidata(hfig);

blank_dist=str2double(get(gdata.rawe1,'string'));
sF1=str2double(get(gdata.rawe2,'string'));
bthick=str2double(get(gdata.rawe3,'string'));

dr=diff(gdata.bindata.range(1:2));
blkind=abs(round(blank_dist/dr));
minwidth=abs(round(bthick/dr));

binlen=size(gdata.bindata.vals,2);
nanind=unique([find(isnan(gdata.bdata.zraw));...
    cell2mat(gdata.edits(:))]);
nanind=nanind(nanind<=binlen);

rdata=gdata.bindata.vals(blkind:end,nanind);

rlog=rdata>sF1;

[m,n]=size(rdata); %#ok
bdepth=nan(n,1);
for i=1:n
    [c,ind]=getchunks(rlog(:,i));
    
    ind2=rlog(ind,i)==1;
    c2=c(ind2);
    ind2=ind(ind2);
    
    
    if numel(ind2)==1
        bdepth(i)=gdata.bindata.range(ind2+blkind);
    else
        bind=ind2(find(c2>minwidth==1,1,'last'));
        if ~isempty(bind)
            bdepth(i)=gdata.bindata.range(bind+blkind);
        end
    end
end



gdata.bdata.zraw(nanind)=bdepth;

if isfield(gdata,'sv')
    if gdata.invert==1
        zraw=-gdata.bdata.zraw;
    end
    zsos=apply_sos_prof(zraw(nanind),gdata.sv.sos_orig,...
        [gdata.sv.depth gdata.sv.sos],gdata.sv.use_mean_sos,...
        gdata.sv.mean_vel)-gdata.manoff;
    if gdata.invert==1
        zsos=-zsos;
    end
    
    gdata.bdata.zc(nanind)=zsos-gdata.bdata.tide(nanind);
    
else
    gdata.bdata.zc(nanind)=((bdepth-...
        gdata.bdata.tide(nanind)).*gdata.invert)-gdata.manoff;
end



set(gdata.l1,'ydata',gdata.bdata.zraw);

gdata.xlimo=[min(gdata.bdata.distance),...
    max(gdata.bdata.distance)];
gdata.ylimo=[min(gdata.bdata.zc),...
    max(gdata.bdata.zc)];

if ~isfield(gdata,'eraw');
    gdata.eraw=gdata.bdata.zraw;
    gdata.ezc=gdata.bdata.zc;
end
if gdata.numedits==0
    set(gdata.menu12,'visible','on');
    set(gdata.menu13,'visible','on');
end
gdata.numedits=gdata.numedits+1;
gdata.edits{gdata.numedits}=nanind;


guidata(hfig,gdata);
setFocus(hfig)
end

%%%%-----------------------------------------------------------------------
function handdigitize(hfig,evnt) %#ok
gdata=guidata(hfig);

set(gdata.rawp1,'backgroundcolor','g');
hold on
button=1;
points=0;
while button==1
    points=points+1;
    [x,y,button]=ginput(1);
    xp(points)=round(x);
    yp(points)=y;
    
    if points==1;
        dh=plot(xp,yp,'ro-','markerfacecolor','r');
    else
        set(dh,'xdata',xp,'ydata',yp)
    end
end
delete(dh);

if numel(xp)>1
    if max(xp)>numel(gdata.bindata.vtime);
        xp(xp==max(xp))=numel(gdata.bdata.vtime);
    end
    if min(xp)<=0;
        xp(xp<=0)=1;
    end
    
    [xc,xind]=unique(xp);
    yc=yp(xind);
    
    dind=min(xp):max(xp);
    gdata.bdata.zraw(dind)=interp1(xc,yc,dind);
    
    set(gdata.l1,'ydata',gdata.bdata.zraw)
    
    
    if ~isfield(gdata,'eraw');
        gdata.eraw=gdata.bdata.zraw;
        gdata.ezc=gdata.bdata.zc;
    end
    if gdata.numedits==0
        set(gdata.menu12,'visible','on');
        set(gdata.menu13,'visible','on');
    end
    gdata.numedits=gdata.numedits+1;
    gdata.edits{gdata.numedits}=dind';
    
    if isfield(gdata,'sv')
        
        if gdata.invert==1
            zraw=-gdata.bdata.zraw(dind);
        end
        zsos=apply_sos_prof(zraw,gdata.sv.sos_orig,...
            [gdata.sv.depth gdata.sv.sos],gdata.sv.use_mean_sos,...
            gdata.sv.mean_vel)-gdata.manoff;
        if gdata.invert==1
            zsos=-zsos;
        end
        
        gdata.bdata.zc(dind)=zsos-gdata.bdata.tide(dind);
        
        
    else
        gdata.bdata.zc(dind)=((gdata.bdata.zraw(dind)-...
            gdata.bdata.tide(dind)).*gdata.invert)-gdata.manoff;
    end
    
    
    
    
    zreal=find(isfinite(gdata.bdata.zc));
    minr=min(zreal);
    maxr=max(zreal);
    gdata.xlimo=[min(gdata.bdata.distance(minr:maxr)),...
        max(gdata.bdata.distance(minr:maxr))];
    gdata.ylimo=[min(gdata.bdata.zc(minr:maxr)),...
        max(gdata.bdata.zc(minr:maxr))];
    
    
    guidata(hfig,gdata);
end
set(gdata.rawp1,'backgroundcolor',[0.9255    0.9137    0.8471]);
end

function [Pxx, freq] = p_welch(x,nfft,Fs,win,noverlap)
% P_WELCH Power Spectral Density estimate via Welch's method.
% [Pxx, freq] = P_WELCH(x,nfft,Fs,win,noverlap)
%
% Version from older Matlab Toolbox, revised for use in Meg
% Palmsten's routine.
%
% The single-sided PSD is returned for real signals.
% No error checking, and limited options.
%
% [Pxx,F] = P_WELCH(X,NFFT,Fs,WIN,NOVERLAP) X is divided into overlap-
% ping sections, then windowed by the WIN parameter, then zero-padded
% to length NFFT.  The magnitude squared of the length NFFT DFTs of the
% sections are averaged to form Pxx.  Pxx is length NFFT/2+1 for NFFT
% even, (NFFT+1)/2 for NFFT odd, or NFFT if the signal X is complex. If
% you specify a scalar for WIN, a Hanning window of that length is
% used.

n = length(x);           % Number of data points
nwind = length(win);     % length of window
if(nwind==1),
    window = hanning(win);
    nwind = win;
else
    window = win;
end
if n < nwind
    x(nwind)=0;
    n=nwind;
end
% Make x a column vector; do this AFTER the zero-padding in case x is a scalar.
x = x(:);

% Number of windows; (k = fix(n/nwind) for noverlap=0)
k = fix((n-noverlap)/(nwind-noverlap));

index = 1:nwind;
KMU = k*norm(window)^2; % Normalizing scale factor ==> asymptotically unbiased
% KMU = k*sum(window)^2;% alt. Nrmlzng scale factor ==> peaks are about right

% Calculate PSD
Spec = zeros(nfft,1);
for i=1:k
    xw = window.*x(index);
    index = index + (nwind - noverlap);
    Xx = abs(fft(xw,nfft)).^2;
    Spec = Spec + Xx;
end

% Select first half
if rem(nfft,2),         % nfft odd
    select = (2:(nfft+1)/2-1)';  % don't include DC or Nyquist components
    nyq    = (nfft+1)/2;         % they're included below
else
    select = (2:nfft/2)';
    nyq    = nfft/2+1;
end
% Calculate the single-sided spectrum which includes the full power
Spec = [Spec(1); 2*Spec(select); Spec(nyq)];
freq = [0; Fs*(1:nfft/2)'/nfft ];

% Normalize, and scale by 1/Fs to get
% Power per unit of frequency
Pxx = Spec / (KMU * Fs);

end

%%%%-----------------------------------------------------------------------
function  bathy_comp(hfig,evnt) %#ok
% The basic layout of this GUI was made with the help of guidegetter,
% available on the File Exchange at Mathworks.com

gdata=guidata(hfig);

gd.mind=1;
gd.fpath1=gdata.outpath;
gd.fpath2=gdata.outpath;
gd.file1=[];
gd.file2=[];

hf = figure('units','normalized','position',[0.253 0.367 0.273 0.261],...
    'menubar','none','name','bathy_comp',...
    'numbertitle','off','color',[0.925 0.914 0.847]);

gd.run = uicontrol(hf,'style','pushbutton',...
    'units','normalized','position',[0.741 0.18 0.177 0.129],...
    'string','Run','backgroundcolor',[0.925 0.914 0.847],...
    'enable','off','callback',@run_zcomp);

gd.trun=uicontrol(hf,'style','text',...
    'units','normalized','position',[-0.01 0.885 0.45 0.1],...
    'string','Program Running. Please wait',...
    'backgroundcolor',[0.925 0.914 0.847],...
    'foregroundcolor','r',...
    'visible','off');


uipanel2 = uipanel('parent',hf,'units','normalized',...
    'position',[0.0277 0.427 0.896 0.473],'title','Files');
uipanel1 = uipanel('parent',hf,'units','normalized',...
    'position',[0.0346 0.12 0.677 0.268],'title','Options');

gd.edit1 = uicontrol(uipanel1,'style','edit','units','normalized',...
    'position',[0.519 0.323 0.388 0.323],...
    'string',sprintf('%0.1f',gd.mind),...
    'backgroundcolor',[1 1 1],...
    'callback',@set_mind);
uicontrol(uipanel1,'style','text',...
    'units','normalized','position',[0.495 0.729 0.422 0.167],...
    'string','Minimum Distance (m)','backgroundcolor',...
    [0.925 0.914 0.847]);

gd.rflag = uicontrol(uipanel1,'style','checkbox','units','normalized',...
    'position',[0.0394 0.406 0.269 0.24],...
    'string','Histogram View','backgroundcolor',[0.925 0.914 0.847],...
    'value',1);
gd.pflag = uicontrol(uipanel1,'style','checkbox','units','normalized',...
    'position',[0.0394 0.125 0.269 0.24],...
    'string','Profile View','backgroundcolor',[0.925 0.914 0.847],...
    'value',0,'callback',@run_lnw_dlg);
gd.mflag = uicontrol(uipanel1,'style','checkbox','units','normalized',...
    'position',[0.0371 0.687 0.269 0.24],'string','Map View',...
    'backgroundcolor',[0.925 0.914 0.847],'value',1);


uicontrol(uipanel2,'style','pushbutton','units',...
    'normalized','position',[0.242 0.307 0.193 0.2],...
    'string','Open File(s)','backgroundcolor',[0.925 0.914 0.847],...
    'callback',@local_open2);
uicontrol(uipanel2,'style','pushbutton','units',...
    'normalized','position',[0.242 0.636 0.193 0.2],...
    'string','Open File(s)','backgroundcolor',[0.925 0.914 0.847],...
    'callback',@local_open1);

gd.string1 = uicontrol(uipanel2,'style','edit','units','normalized',...
    'position',[0.0229 0.629 0.21 0.207],...
    'string','Vessel 1','backgroundcolor',[1 1 1]);
gd.string2 = uicontrol(uipanel2,'style','edit','units','normalized',...
    'position',[0.0267 0.271 0.206 0.257],'string','Vessel 2',...
    'backgroundcolor',[1 1 1]);


gd.list2 = uicontrol(uipanel2,'style','listbox','units','normalized',...
    'position',[0.484 0.343 0.464 0.179],...
    'string','No Files Selected','backgroundcolor',[1 1 1]);
gd.list1 = uicontrol(uipanel2,'style','listbox','units','normalized',...
    'position',[0.484 0.614 0.464 0.179],...
    'string','No Files Selected','backgroundcolor',[1 1 1]);

guidata(hf,gd);
end

%%%------------------------------------------------------------------------
function local_open1(hf,evnt) %#ok

gd=guidata(hf);

[filename, pathname] = uigetfile( ...
    {'*.xyz', 'XYZ Files (*.xyz)';...
    '*.*', 'All Files (*.*)'},...
    'Select .xyz file(s)',...
    'multiselect','on',...
    gd.fpath1);

if iscell(filename)
else
    if filename==0
        return
    end
end

gd.fpath1=pathname;
gd.file1=filename;


set(gd.list1,'string',gd.file1);

if iscell(gd.file1)
    gd.files1=cellfun(@(x)([gd.fpath1,x]),gd.file1','un',0);
else
    gd.files1={[gd.fpath1,gd.file1]};
end

if ~isempty(gd.file1) && ~isempty(gd.file2)
    [~,fidx]=intersect(gd.files1,gd.files2);
    if ~isempty(fidx)
        errordlg('Replicate files chosen. Please choose another file.')
    else
        set(gd.run,'enable','on')
    end
end


guidata(hf,gd);
end

%%%------------------------------------------------------------------------
function local_open2(hf,evnt) %#ok

gd=guidata(hf);

[filename, pathname] = uigetfile( ...
    {'*.xyz', 'XYZ Files (*.xyz)';...
    '*.*', 'All Files (*.*)'},...
    'Select .xyz file(s)',...
    'multiselect','on',...
    gd.fpath2);

if iscell(filename)
else
    if filename==0
        return
    end
end

gd.fpath2=pathname;
gd.file2=filename;

set(gd.list2,'string',gd.file2);


if iscell(gd.file2)
    gd.files2=cellfun(@(x)([gd.fpath2,x]),gd.file2','un',0);
else
    gd.files2={[gd.fpath2,gd.file2]};
end



if ~isempty(gd.file1) && ~isempty(gd.file2)
    [~,fidx]=intersect(gd.files1,gd.files2);
    if ~isempty(fidx)
        errordlg('Replicate files chosen. Please choose another file.')
    else
        set(gd.run,'enable','on')
    end
end

guidata(hf,gd);
end

%%%------------------------------------------------------------------------
function set_mind(hf,evnt) %#ok

gd=guidata(hf);

mind=str2double(get(gd.edit1,'string'));

if mind<0;
    errordlg('Minimum distance should be postive');
else
    gd.mind=mind;
    set(gd.edit1,'string',sprintf('%0.1f',gd.mind));
end
guidata(hf,gd);
end

%%%%-----------------------------------------------------------------------
function run_lnw_dlg(hf,evnt) %#ok

gd=guidata(hf);

val=get(gd.pflag,'value');

if val==1
    if isfield(gd,'ldata')
        gd.ldata=lnw_dlg(gd.ldata);
    else
        gd.ldata=lnw_dlg;
    end
    set(gd.pflag,'value',1);
    guidata(hf,gd)
end
end

%%%------------------------------------------------------------------------
function run_zcomp(hf,evnt) %#ok

gd=guidata(hf);

if all([get(gd.pflag,'value');...
        get(gd.rflag,'value');...
        get(gd.mflag,'value')]==0)
    errordlg('Please select at least one graphical output.');
    return
end

set(gd.trun,'visible','on')
drawnow;

data1=cell2mat(cellfun(@(x)(load(x)),gd.files1,'un',0));
data1(any(isnan(data1)==1,2),:)=[];

data2=cell2mat(cellfun(@(x)(load(x)),gd.files2,'un',0));
data2(any(isnan(data2)==1,2),:)=[];

%compute the interpoint distances
d=ipdm([data1(:,1) data1(:,2)],[data2(:,1) data2(:,2)],'subset','nearest',...
    'result','structure');

%find the nearest unique point
subs=[d.rowindex d.columnindex];
[dsort,didx]=sort(d.distance);
subs=subs(didx,:);

[subs2,idx]=unique(subs(:,2),'first');
subs=[subs(idx,1) subs2];

%limit the pairs of points within 1 m
subs(dsort(idx)>gd.mind,:)=[];


if isempty(subs)
    errordlg('No overlapping points found.')
    set(gd.trun,'visible','off')
    return
end


x1=data1(subs(:,1),1);
y1=data1(subs(:,1),2);
z1=data1(subs(:,1),3);

z2=data2(subs(:,2),3);

dz=z2-z1;

if get(gd.mflag,'value')==1
    figure('name','Map View')
    plot(data1(:,1),data1(:,2),'.','color',[0.6 0.6 0.6])
    hold on
    plot(data2(:,1),data2(:,2),'.','color','k')
    axis equal
    
    for i = 1 : length(subs(:,1))
        plot([data2(subs(i,2),1);...
            data1(subs(i,1),1)],...
            [data2(subs(i,2),2);...
            data1(subs(i,1),2)],'ro-');
    end
    
    legstr={get(gd.string1,'string'),...
        get(gd.string2,'string'),...
        'Point Pairs'};
    legend(legstr,'location','best',...
        'box','off')
    xlabel('X','fontsize',12)
    ylabel('Y','fontsize',12)
end

if get(gd.rflag,'value')==1
    
    
    
    
    mean_dz=mean(dz);
    sk=modskill(z2,z1);
    sk.mean_offset=mean_dz;
    
    
    q5=interp1(linspace(0,100,numel(dz)),sort(dz),1);
    q95=interp1(linspace(0,100,numel(dz)),sort(dz),99);
    
    xi=linspace(q5,q95,30);
    [n,bin]=histc(dz,xi);
    
    figure('name','Histogram View')
    ax(1)=subplot(1,2,1);
    hold on
    lh(2)=plot(z1,z2,'k.');
    set(gca,'xlim',[min([z1;z2]) max([z1;z2])],...
        'ylim',[min([z1;z2]) max([z1;z2])],...
        'da',[1 1 1],'box','on');
    lh(1)=line([min([z1;z2]) max([z1;z2])],...
        [min([z1;z2]) max([z1;z2])],...
        'color','r','linewi',2);
    
    xlabel(['\bf\it\fontsize{12}',get(gd.string1,'string')])
    ylabel(['\bf\it\fontsize{12}',get(gd.string2,'string')])
    
    xl=xlim;
    yl=ylim;
    xp=repmat(xl(1)+0.025*diff(xl),1,5);
    yp=yl(1)+(fliplr((0.75:0.05:0.95)).*diff(yl));
    fields={'n';'rmse';'rmsa';'rmsu';'mean_offset'};
    data=cell(1,length(fields));
    fmt={'%d','%.2f','%.3f','%.2f','%.2f'};
    for i=1:length(fields)
        data{i}=sk.(fields{i});
    end
    tv=cellfun(@(x,y,z)(sprintf([x,' = ',y],z)),...
        fields',fmt,data,'un',0)';
    text(xp,yp,tv,'fontsize',10,...
        'fontwei','b','fontan','it',...
        'interpreter','none')
    
    
    subplot(122)
    bh=bar(xi,(n./numel(bin)).*100,'histc');
    set(bh,'facecolor',[0.6 0.6 0.6])
    
    set(gca,'xlim',[q5 q95],...
        'plotboxaspectratio',...
        get(ax(1),'plotboxaspectratio'))
    xlabel(['\bf\it\fontsize{12}',...
        get(gd.string2,'string'),' - ',...
        get(gd.string1,'string')])
    ylabel('\bf\it\fontsize{12}Percent')
end


if get(gd.pflag,'value')==1
    
    figure('Name','Profile View');
    ax(1)=subplot(3,1,1:2);
    lh(1)=plot((1:size(data1,1)),data1(:,3),'.','color',[0.6 0.6 0.6]);
    hold on
    lh(2)=plot((1:length(x1)),data1(subs(:,1),3),'ro');
    lh(3)=plot((1:size(data2,1)),data2(:,3),'.','color','k');
    ylabel('\bf\itElevation (m)')
    
    legstr={get(gd.string1,'string'),...
        get(gd.string2,'string'),...
        'Point Pairs'};
    legend([lh(1);lh(3);lh(2)],legstr,'location','best',...
        'box','off')
    
    ax(2)=subplot(313);
    lh(4)=plot((1:length(x1)),dz,'r.');
    ylabel('\bf\itElevation Difference(m)')
    
    switch gd.ldata.ref
        case 1
            set(lh(1),'xdata',data1(:,1));
            set(lh(2),'xdata',x1);
            set(lh(3),'xdata',data2(:,1));
            set(lh(4),'xdata',x1);
            
            minx=min([data1(:,1);data2(:,1)]);
            maxx=max([data1(:,1);data2(:,1)]);
            set(ax,'xlim',[minx maxx]);
            line([minx maxx],[0 0],'color','k')
            linkaxes(ax,'x')
            
            xlabel('\bf\itEasting (m)')
            
        case 2
            
            set(lh(1),'xdata',data1(:,2));
            set(lh(2),'xdata',y1);
            set(lh(3),'xdata',data2(:,2));
            set(lh(4),'xdata',y1);
            
            minx=min([data1(:,2);data2(:,2)]);
            maxx=max([data1(:,2);data2(:,2)]);
            set(ax,'xlim',[minx maxx]);
            line([minx maxx],[0 0],'color','k')
            linkaxes(ax,'x')
            
            xlabel('\bf\itNorthing (m)')
            
        case 3
            
            lnw=gd.ldata.lnw;
            lind=gd.ldata.idx;
            
            
            
            lx=lnw(lind).x;
            ly=lnw(lind).y;
            m=(ly(end)-ly(1))/...
                (lx(end)-lx(1));
            theta=atan(m);
            
            xt1 = (data1(:,1).*cos(theta) + data1(:,2).*sin(theta));
            xt2 = (data2(:,1).*cos(theta) + data2(:,2).*sin(theta));
            o1 = (x1.*cos(theta) + y1.*sin(theta));
            
            
            x0=lx.*cos(theta) + ...
                ly.*sin(theta);
            
            dist1=xt1-x0(1);
            dist2=xt2-x0(1);
            odist1=o1-x0(1);
            
            
            set(lh(1),'xdata',dist1);
            set(lh(2),'xdata',odist1);
            set(lh(3),'xdata',dist2);
            set(lh(4),'xdata',odist1);
            
            minx=min([dist1;dist2]);
            maxx=max([dist1;dist2]);
            set(ax,'xlim',[minx maxx]);
            line([minx maxx],[0 0],'color','k')
            linkaxes(ax,'x')
            
            xlabel('\bf\itCross-shore Distance (m)')
            
            
            
    end
end

set(gd.trun,'visible','off')

end

%%%-----------------------------------------------------------------------
function line_comp = lnw_dlg(varargin)

if nargin>0;
    gd.line_comp=varargin{1};
    lnwstr=gd.line_comp.lnwfile;
else
    gd.line_comp.lnw=[];
    lnwstr='Open File';
end

hf2 = figure('units','normalized','position',[0.253 0.467 0.17 0.158],...
    'menubar','none','name','lnw_dlg',...
    'numbertitle','off','color',[0.925 0.914 0.847]);


gd.close=uicontrol(hf2,'style','pushbutton','units','normalized',...
    'position',[0.341 0.0594 0.295 0.119],'string','OK',...
    'backgroundcolor',[0.925 0.914 0.847],...
    'callback',@lnw_close);
gd.lopen=uicontrol(hf2,'style','pushbutton','units',...
    'normalized','position',[0.45 0.49 0.318 0.109],...
    'string',lnwstr,'backgroundcolor',[0.925 0.914 0.847],...
    'callback',@open_lnw,'enable','off');

uicontrol(hf2,'style','text','units','normalized',...
    'position',[0.0401 0.441 0.312 0.139],'string','Line File',...
    'backgroundcolor',[0.925 0.914 0.847]);
gd.linen = uicontrol(hf2,'style','popupmenu','units','normalized',...
    'position',[0.453 0.282 0.315 0.129],'string','None',...
    'backgroundcolor',[1 1 1],'enable','off',...
    'callback',@select_line);

uicontrol(hf2,'style','text','units','normalized',...
    'position',[0.0573 0.302 0.289 0.0693],...
    'string','Line Number','backgroundcolor',[0.925 0.914 0.847]);
gd.popup1 = uicontrol(hf2,'style','popupmenu','units','normalized',...
    'position',[0.453 0.688 0.318 0.114],'string',...
    {'X';'Y';'Line File'},'backgroundcolor',[1 1 1],...
    'callback',@listmotion);
uicontrol(hf2,'style','text','units','normalized',...
    'position',[0.0458 0.649 0.312 0.139],...
    'string','X-Axis Reference','backgroundcolor',[0.925 0.914 0.847]);
guidata(hf2,gd);

uiwait
gd=guidata(hf2);
line_comp=gd.line_comp;

close(hf2)
end
%%%%----------------------------------------------------------------------
function lnw_close(hf2,evnt) %#ok
gd=guidata(hf2);
value=get(gd.popup1,'value');
switch value
    case 1
        gd.line_comp.ref=1;
    case 2
        gd.line_comp.ref=2;
    case 3
        gd.line_comp.ref=3;
        gd.line_comp.lnw=gd.line_comp.lnw;
        gd.line_comp.idx=get(gd.linen,'value');
end
guidata(hf2,gd);

uiresume
end
%%%%----------------------------------------------------------------------
function listmotion(hf2,evnt) %#ok
gd=guidata(hf2);

value=get(gd.popup1,'value');
switch value
    case 3
        set(gd.lopen,'enable','on');
        if isempty(gd.line_comp.lnw)
            set(gd.close,'enable','off');
        else
            set(gd.linen,'string',{gd.line_comp.lnw(:).name}')
            set(gd.linen,'enable','on')
        end
    otherwise
        set(gd.linen,'enable','off');
        set(gd.lopen,'enable','off');
end
end
%%%%----------------------------------------------------------------------
function open_lnw(hf2,evnt) %#ok

gd=guidata(hf2);
[filename, pathname] = ...
    uigetfile('*.lnw', 'Pick an LNW-file');
fpath=[pathname,filename];
if filename==0
    return
end
gd.line_comp.lnw=readLNW('filename',fpath);
set(gd.linen,'enable','on');
set(gd.lopen,'string',filename);
gd.line_comp.lnwfile=filename;
set(gd.linen,'string',{gd.line_comp.lnw(:).name}')
guidata(hf2,gd);
end
%%%%----------------------------------------------------------------------
function select_line(hf2,evnt) %#ok
gd=guidata(hf2);
set(gd.close,'enable','on');
guidata(hf2,gd);

end
%%%%----------------------------------------------------------------------
function d = ipdm(data1,varargin)
% ipdm: Inter-Point Distance Matrix
% usage: d = ipdm(data1)
% usage: d = ipdm(data1,data2)
% usage: d = ipdm(data1,prop,value)
% usage: d = ipdm(data1,data2,prop,value)
%
% Arguments: (input)
%  data1 - array of data points, each point is one row. p dimensional
%          data will be represented by matrix with p columns.
%          If only data1 is provided, then the distance matrix
%          is computed between all pairs of rows of data1.
%
%          If your data is one dimensional, it MUST form a column
%          vector. A row vector of length n will be interpreted as
%          an n-dimensional data set.
%
%  data2 - second array, supplied only if distances are to be computed
%          between two sets of points.
%
%
% Class support: data1 and data2 are assumed to be either
% single or double precision. I have not tested this code to
% verify its success on integer data of any class.
%
%
% Additional parameters are expected to be property/value pairs.
% Property/value pairs are pairs of arguments, the first of which
% (properties) must always be a character string. These strings
% may be shortened as long as the shortening is unambiguous.
% Capitalization is ignored. Valid properties for ipdm are:
%
%  'Metric', 'Subset', 'Limit', 'Result'
%
%  'Metric' - numeric flag - defines the distance metric used
%          metric = 2 --> (DEFAULT) Euclidean distance = 2-norm
%                         The standard distance metric.
%
%          metric = 1 --> 1-norm = sum of absolute differences
%                         Also sometimes known as the "city block
%                         metric", since this is the sum of the
%                         differences in each dimension.
%
%          metric = inf --> infinity-norm = maximum difference
%                         over all dimensions. The name refers
%                         to the limit of the p-norm, as p
%                         approaches infinity.
%
%          metric = 0 --> minimum difference over all dimensions.
%                         This is not really a useful norm in
%                         practice.
%
%          Note: while other distance metrics exist, IMHO, these
%          seemed to be the common ones.
%
%
%  'Result' - A string variable that denotes the style of returned
%          result. Valid result types are 'Array', 'Structure'.
%          Capitalization is ignored, and the string may be
%          shortened if you wish.
%
%          result = 'Array' --> (DEFAULT) A matrix of all
%                         interpoint distances will be generated.
%                         This array may be large. If this option
%                         is specified along with a minimum or
%                         maximum value, then those elements above
%                         or below the limiting values will be
%                         set as -inf or +inf, as appropriate.
%
%                         When any of 'LargestFew', 'SmallestFew',
%                         or 'NearestNeighbor' are set, then the
%                         resulting array will be a sparse matrix
%                         if 'array' is specified as the result.
%
%          result = 'Structure' --> A list of all computed distances,
%                         defined as a structure. This structure
%                         will have fields named 'rowindex',
%                         'columnindex', and 'distance'.
%
%                         This option will be useful when a subset
%                         criterion for the distances has been
%                         specified, since then the distance matrix
%                         may be very sparsely populated. Distances
%                         for pairs outside of the criterion will
%                         not be returned.
%
%
%  'Subset' - Character string, any of:
%
%          'All', 'Maximum', 'Minimum', 'LargestFew', 'SmallestFew',
%          'NearestNeighbor', 'FarthestNeighbor', or empty
%
%          Like properties, capitalization is ignored here, and
%          any unambiguous shortening of the word is acceptable.
%
%          DEFAULT = 'All'
%
%          Some interpoint distance matrices can be huge. Often
%          these matrices are too large to be fully retained in
%          memory, yet only the pair of points with the largest
%          or smallest distance may be needed. When only some
%          subset of the complete set of distances is of interest,
%          these options allow you to specify which distances will
%          be returned.
%
%          If 'result' is defined to be an array, then a sparse
%          matrix will be returned for the 'LargestFew', 'SmallestFew',
%          'NearestNeighbor', and 'FarthestNeighbor' subset classes.
%          'Minimum' and 'Maximum' will yield full matrices by
%          default. If a structure is specified, then only those
%          elements which have been identified will be returned.
%
%          Where a subset is specified, its limiting value is
%          specified by the 'Limit' property. Call that value k.
%
%
%          'All' -->     (DEFAULT) Return all interpoint distances
%
%          'Minimum' --> Only look for those distances above
%                        the cutoff k. All other distances will
%                        be returned as -inf.
%
%          'Maximum' --> Only look for those distances below
%                        the cutoff k. All other distances will
%                        be returned as +inf.
%
%          'SmallestFew' --> Only return the subset of the k
%                        smallest distances. Where only one data
%                        set is provided, only the upper triangle
%                        of the inter-point distance matrix will
%                        be generated since that matrix is symmetric.
%
%          'LargestFew' --> Only return the subset of the k
%                        largest distances. Where only one data
%                        set is provided, only the upper triangle
%                        of the inter-point distance matrix will
%                        be generated since that matrix is symmetric.
%
%          'NearestNeighbor' --> Only return the single nearest
%                        neighbor in data2 to each point in data1.
%                        No limiting value is required for this
%                        option. If multiple points have the same
%                        nearest distance, then return the first
%                        such point found. With only one input set,
%                        a point will not be its own nearest
%                        neighbor.
%
%                        Note that exact replicates in a single set
%                        will cause problems, since a sparse matrix
%                        is returned by default. Since they will have
%                        a zero distance, they will not show up in
%                        the sparse matrix. A structure return will
%                        show those points as having a zero distance
%                        though.
%
%          'FarthestNeighbor' --> Only return the single farthest
%                        neighbor to each point. No limiting value
%                        is required for this option. If multiple
%                        points have the same farthest distance,
%                        then return the first such point found.
%
%
%  'Limit' - scalar numeric value or []. Used only when some
%           Subset is specified.
%
%          DEFAULT = []
%
%
%  'ChunkSize' - allows a user with lower RAM limits
%          to force the code to only grab smaller chunks of RAM
%          at a time (where possible). This parameter is specified
%          in bytes of RAM. The default is 32 megabytes, or 2^22
%          elements in any piece of the distance matrix. Only some
%          options will break the problem into chunks, thus as long
%          as a full matrix is expected to be returned, there seems
%          no reason to break the problem up into pieces.
%
%          DEFAULT = 2^25
%
%
% Arguments: (output)
%  d     - array of interpoint distances, or a struct wth the
%          fields {'rowindex', 'columnindex', 'distance'}.
%
%          d(i,j) represents the distance between point i
%          (from data1) and point j (from data2).
%
%          If only one (n1 x p) array is supplied, then d will
%          be an array of size == [n1,n1].
%
%          If two arrays (of sizes n1 x p and n2 x p) then d
%          will be an array of size == [n1,n2].
%
%
% Efficiency considerations:
%  Where possible, this code will use bsxfun to compute its
%  distances.
%
%
% Example:
%  Compute the interpoint distances between all pairs of points
%  in a list of 5 points, in 2 dimensions and using Euclidean
%  distance as the distance metric.
%
%  A = randn(5,2);
%  d = ipdm(A,'metric',2)
%  d =
%            0       2.3295       3.2263       2.0263       2.8244
%       2.3295            0       1.1485      0.31798       1.0086
%       3.2263       1.1485            0       1.4318       1.8479
%       2.0263      0.31798       1.4318            0       1.0716
%       2.8244       1.0086       1.8479       1.0716            0
%
% (see the demo file for many other examples)
%
% See also: pdist
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 2/26/08

% Default property values
params.Metric = 2;
params.Result = 'array';
params.Subset = 'all';
params.Limit = [];
params.ChunkSize = 2^25;

% untangle the arguments
if nargin<1
    % if called with no arguments, then the user probably
    % needs help. Give it to them.
    help ipdm
    return
end

% were two sets of data provided?
pvpairs = {};
if nargin==1
    % only 1 set of data provided
    dataflag = 1;
    data2 = [];
else
    if ischar(varargin{1})
        dataflag = 1;
        data2 = [];
        pvpairs = varargin;
    else
        dataflag = 2;
        data2 = varargin{1};
        if nargin>2
            pvpairs = varargin(2:end);
        end
    end
end

% get data sizes for later
[n1,dim] = size(data1);
if dataflag == 2
    n2 = size(data2,1);
end

% Test the class of the input variables
if ~(isa(data1,'double') || isa(data1,'single')) || ...
        ((dataflag == 2) && ~(isa(data2,'double') || isa(data2,'single')))
    error('data points must be either single or double precision variables.')
end

% do we need to process any property/value pairs?
if nargin>2
    params = parse_pv_pairs(params,pvpairs);
    
    % check for problems in the properties
    
    % was a legal Subset provided?
    if ~isempty(params.Subset) && ~ischar(params.Subset)
        error('If provided, ''Subset'' must be character')
    elseif isempty(params.Subset)
        params.Subset = 'all';
    end
    valid = {'all','maximum','minimum','largestfew','smallestfew', ...
        'nearestneighbor','farthestneighbor'};
    ind = find(strncmpi(params.Subset,valid,length(params.Subset)));
    if (length(ind)==1)
        params.Subset = valid{ind};
    else
        error(['Invalid Subset: ',params.Subset])
    end
    
    % was a limit provided?
    if ~ismember(params.Subset,{'all','nearestneighbor','farthestneighbor'}) && ...
            isempty(params.Limit)
        error('No limit provided, but a Subset that requires a limit value was specified')
    end
    % check the limit values for validity
    if length(params.Limit)>1
        error('Limit must be scalar or empty')
    end
    
    switch params.Subset
        case {'largestfew', 'smallestfew'}
            % must be at least 1, and an integer
            if (params.Limit<1) || (round(params.Limit)~=params.Limit)
                error('Limit must be a positive integer for LargestFew or NearestFew')
            end
    end
    
    % was a legal Result provided?
    if isempty(params.Result)
        params.result = 'Array';
    elseif ~ischar(params.Result)
        error('If provided, ''Result'' must be character or empty')
    end
    valid = {'array','structure'};
    ind = find(strncmpi(params.Result,valid,length(params.Result)));
    if (length(ind)==1)
        params.Result = valid{ind};
    else
        error(['Invalid Result: ',params.Subset])
    end
    
    % check for the metric
    if isempty(params.Metric)
        params.Metric = 2;
    elseif (length(params.Metric)~=1) || ~ismember(params.Metric,[0 1 2 inf])
        error('If supplied, ''Metric'' must be a scalar, and one of [0 1 2 inf]')
    end
end % if nargin>2

% If Metric was given as 2, but the dimension is only 1, then it will
% be slightly faster (and equivalent) to use the 1-norm Metric.
if (dim == 1) && (params.Metric == 2)
    params.Metric = 1;
end

% Can we use bsxfun to compute the interpoint distances?
% Older Matlab releases will not have bsxfun, but if it is
% around, it will ne both faster and less of a memory hog.
params.usebsxfun = (5==exist('bsxfun','builtin'));

% check for dimension mismatch if 2 sets
if (dataflag==2) && (size(data2,2)~=dim)
    error('If 2 point sets provided, then both must have the same number of columns')
end

% Total number of distances to compute, in case I must do it in batches
if dataflag==1
    n2 = n1;
end
ntotal = n1*n2;

% FINALLY!!! Compute inter-point distances
switch params.Subset
    case 'all'
        % The complete set of interpoint distances. There is no need
        % to break this into chunks, since we must return all distances.
        % If that is too much to compute in memory, then it will fail
        % anyway when we try to store the result. bsxfun will at least
        % do the computation efficiently.
        
        % One set or two?
        if dataflag == 1
            d = distcomp(data1,data1,params);
        else
            d = distcomp(data1,data2,params);
        end
        
        % Must we return it as a struct?
        if params.Result(1) == 's'
            [rind,cind] = ndgrid(1:size(d,1),1:size(d,2));
            ds.rowindex = rind(:);
            ds.columnindex = cind(:);
            ds.distance = d(:);
            d = ds;
        end
        
    case {'minimum' 'maximum'}
        % There is no reason to break this into pieces if the result
        % sill be filled in the end with +/- inf. Only break it up
        % if the final result is a struct.
        if ((ntotal*8)<=params.ChunkSize) || (params.Result(1) == 'a')
            % its small enough to do it all at once
            
            % One set or two?
            if dataflag == 1
                d = distcomp(data1,data1,params);
            else
                d = distcomp(data1,data2,params);
            end
            
            % Must we return it as a struct?
            if params.Result(1) == 'a'
                % its an array, fill the unwanted distances with +/- inf
                if params.Subset(2) == 'i'
                    % minimum
                    d(d<=params.Limit) = -inf;
                else
                    % maximum
                    d(d>=params.Limit) = +inf;
                end
            else
                % a struct will be returned
                if params.Subset(2) == 'i'
                    % minimum
                    [dist.rowindex,dist.columnindex] = find(d>=params.Limit);
                else
                    % maximum
                    [dist.rowindex,dist.columnindex] = find(d<=params.Limit);
                end
                dist.distance = d(dist.rowindex + n1*(dist.columnindex-1));
                d = dist;
            end
            
        else
            % we need to break this into chunks. This branch
            % will always return a struct.
            
            % this is the number of rows of data1 that we will
            % process at a time.
            bs = floor(params.ChunkSize/(8*n2));
            bs = min(n1,max(1,bs));
            
            % Accumulate the result into a cell array. Do it this
            % way because we don't know in advance how many elements
            % that we will find satisfying the minimum or maximum
            % limit specified.
            accum = cell(0,1);
            
            % now loop over the chunks
            batch = 1:bs;
            while ~isempty(batch)
                
                % One set or two?
                if dataflag == 1
                    dist = distcomp(data1(batch,:),data1,params);
                else
                    dist = distcomp(data1(batch,:),data2,params);
                end
                
                % big or small as requested
                if ('i'==params.Subset(2))
                    % minimum value specified
                    [I,J,V] = find(dist>=params.Limit);
                else
                    % maximum limit
                    [I,J] = find(dist<=params.Limit);
                    I = I(:);
                    J = J(:);
                    V = dist(I + (J-1)*length(batch));
                    I = I + (batch(1)-1);
                end
                
                % and stuff them into the cell structure
                if ~isempty(V)
                    accum{end+1,1} = [I,J,V(:)]; %#ok
                end
                
                % increment the batch
                batch = batch + bs;
                if batch(end)>n1
                    batch(batch>n1) = [];
                end
                
            end
            
            % convert the cells into one flat array
            accum = cell2mat(accum);
            
            if isempty(accum)
                d.rowindex = [];
                d.columnindex = [];
                d.distance = [];
            else
                % we found something
                
                % sort on the second column, to put them in a reasonable order
                accum = sortrows(accum,[2 1]);
                
                d.rowindex = accum(:,1);
                d.columnindex = accum(:,2);
                d.distance = accum(:,3);
            end
            
        end
        
    case {'smallestfew' 'largestfew'}
        % find the k smallest/largest distances. k is
        % given by params.Limit
        
        % if only 1 set, params.Limit must be less than n*(n-1)/2
        if dataflag == 1
            params.Limit = min(params.Limit,n1*(n1-1)/2);
        end
        
        % is this a large problem?
        if ((ntotal*8) <= params.ChunkSize)
            % small potatoes
            
            % One set or two?
            if dataflag == 1
                dist = distcomp(data1,data1,params);
                % if only one data set, set the diagonal and
                % below that to +/- inf so we don't find it.
                temp = find(tril(ones(n1,n1),0));
                if params.Subset(1) == 's'
                    dist(temp) = inf;
                else
                    dist(temp) = -inf;
                end
            else
                dist = distcomp(data1,data2,params);
            end
            
            % sort the distances to find those we need
            if ('s'==params.Subset(1))
                % smallestfew
                [val,tags] = sort(dist(:),'ascend');
            else
                % largestfew
                [val,tags] = sort(dist(:),'descend');
            end
            val = val(1:params.Limit);
            tags = tags(1:params.Limit);
            
            % recover the row and column index from the linear
            % index returned by sort in tags.
            [d.rowindex,d.columnindex] = ind2sub([n1,size(dist,2)],tags);
            
            % create the matrix as a sparse one or a struct?
            if params.Result(1)=='a'
                % its an array, so make the array sparse.
                d = sparse(d.rowindex,d.columnindex,val,n1,size(dist,2));
            else
                % a structure
                d.distance = val;
            end
            
        else
            % chunks
            
            % this is the number of rows of data1 that we will
            % process at a time.
            bs = floor(params.ChunkSize/(8*n2));
            bs = min(n1,max(1,bs));
            
            % We need to find the extreme cases. There are two possible
            % algorithms, depending on how many total elements we will
            % search for.
            % 1. Only a very few total elements.
            % 2. A relatively large number of total elements, forming
            %    a significant fraction of the total set.
            %
            % Case #1 would suggest to retain params.Limit numberr of
            % elements from each batch, then at the end, sort them all
            % to find the best few. Case #2 will result in too many
            % elements to retain, so we must distinguish between these
            % alternatives.
            if (8*params.Limit*n1/bs) <= params.ChunkSize
                % params.Limit is small enough to fall into case #1.
                
                % Accumulate the result into a cell array. Do it this
                % way because we don't know in advance how many elements
                % that we will find satisfying the minimum or maximum
                % limit specified.
                accum = cell(0,1);
                
                % now loop over the chunks
                batch = (1:bs)';
                while ~isempty(batch)
                    % One set or two?
                    if dataflag == 1
                        dist = distcomp(data1(batch,:),data1,params);
                        k = find(tril(ones(length(batch),n2),batch(1)-1));
                        if ('s'==params.Subset(1))
                            dist(k) = inf;
                        else
                            dist(k) = -inf;
                        end
                    else
                        dist = distcomp(data1(batch,:),data2,params);
                    end
                    
                    % big or small as requested, keeping only the best
                    % params.Limit number of elements
                    if ('s'==params.Subset(1))
                        % minimum value specified
                        [tags,tags] = sort(dist(:),1,'ascend');
                        tags = tags(1:bs);
                        [I,J] = ndgrid(batch,1:n2);
                        ijv = [I(tags),J(tags),dist(tags)];
                    else
                        % maximum limit
                        [tags,tags] = sort(dist(:),1,'descend');
                        tags = tags(1:bs);
                        [I,J] = ndgrid(batch,1:n2);
                        ijv = [I(tags),J(tags),dist(tags)];
                    end
                    % and stuff them into the cell structure
                    accum{end+1,1} = ijv; %#ok
                    
                    % increment the batch
                    batch = batch + bs;
                    if batch(end)>n1
                        batch(batch>n1) = [];
                    end
                end
                
                % convert the cells into one flat array
                accum = cell2mat(accum);
                
                % keep only the params.Limit best of those singled out
                accum = sortrows(accum,3);
                if ('s'==params.Subset(1))
                    % minimum value specified
                    accum = accum(1:params.Limit,:);
                else
                    % minimum value specified
                    accum = accum(end + 1 - (1:params.Limit),:);
                end
                d.rowindex = accum(:,1);
                d.columnindex = accum(:,2);
                d.distance = accum(:,3);
                
                % create the matrix as a sparse one or a struct?
                if params.Result(1)=='a'
                    % its an array, so make the array sparse.
                    d = sparse(d.rowindex,d.columnindex,d.distance,n1,size(dist,2));
                end
                
            else
                % params.Limit forces us into the domain of case #2.
                % Here we cannot retain params.Limit elements from each chunk.
                % so we will grab each chunk and append it to the best elements
                % found so far, then filter out the best after each chunk is
                % done. This may be slower than we want, but its the only way.
                ijv = zeros(0,3);
                
                % loop over the chunks
                batch = (1:bs)';
                while ~isempty(batch)
                    % One set or two?
                    if dataflag == 1
                        dist = distcomp(data1(batch,:),data1,params);
                        k = find(tril(ones(length(batch),n2),batch(1)-1));
                        if ('s'==params.Subset(1))
                            dist(k) = inf;
                        else
                            dist(k) = -inf;
                        end
                    else
                        dist = distcomp(data1(batch,:),data2,params);
                    end
                    
                    [I,J] = ndgrid(batch,1:n2);
                    ijv = [ijv;[I(:),J(:),dist(:)]]; %#ok
                    
                    % big or small as requested, keeping only the best
                    % params.Limit number of elements
                    if size(ijv,1) > params.Limit
                        if ('s'==params.Subset(1))
                            % minimum value specified
                            [tags,tags] = sort(ijv(:,3),1,'ascend');
                        else
                            [tags,tags] = sort(ijv(:,3),1,'ascend');
                        end
                        ijv = ijv(tags(1:params.Limit),:);
                    end
                    
                    % increment the batch
                    batch = batch + bs;
                    if batch(end)>n1
                        batch(batch>n1) = [];
                    end
                end
                
                % They are fully trimmed down. stuff a structure
                d.rowindex = ijv(:,1);
                d.columnindex = ijv(:,2);
                d.distance = ijv(:,3);
                
                % create the matrix as a sparse one or a struct?
                if params.Result(1)=='a'
                    % its an array, so make the array sparse.
                    d = sparse(d.rowindex,d.columnindex,d.distance,n1,size(dist,2));
                end
                
            end
            
        end
        
    case {'nearestneighbor' 'farthestneighbor'}
        % find the closest/farthest neighbor for every point
        
        % is this a large problem? Or a 1-d problem?
        if dim == 1
            % its a 1-d nearest/farthest neighbor problem. we can
            % special case these easily enough, and all the distance
            % metric options are the same in 1-d.
            
            % first split it into the farthest versus nearest cases.
            if params.Subset(1) == 'f'
                % farthest away
                
                % One set or two?
                if dataflag == 1
                    [d2min,minind] = min(data1);
                    [d2max,maxind] = max(data1);
                else
                    [d2min,minind] = min(data2);
                    [d2max,maxind] = max(data2);
                end
                
                d.rowindex = (1:n1)';
                d.columnindex = repmat(maxind,n1,1);
                d.distance = repmat(d2max,n1,1);
                
                % which endpoint was further away?
                k = abs((data1 - d2min)) >= abs((data1 - d2max));
                if any(k)
                    d.columnindex(k) = minind;
                    d.distance(k) = d2min;
                end
                
            else
                % nearest. this is mainly a sort and some fussing around.
                d.rowindex = (1:n1)';
                d.columnindex = ones(n1,1);
                d.distance = zeros(n1,1);
                
                % One set or two?
                if dataflag == 1
                    % if only one data point, then we are done
                    if n1 == 2
                        % if exactly two data points, its trivial
                        d.columnindex = [2 1];
                        d.distance = repmat(abs(diff(data1)),2,1);
                    elseif n1>2
                        % at least three points. do a sort.
                        [sorted_data,tags] = sort(data1);
                        
                        % handle the first and last points separately
                        d.columnindex(tags(1)) = tags(2);
                        d.distance(tags(1)) = sorted_data(2) - sorted_data(1);
                        d.columnindex(tags(end)) = tags(end-1);
                        d.distance(tags(end)) = sorted_data(end) - sorted_data(end-1);
                        
                        ind = (2:(n1-1))';
                        
                        d1 = sorted_data(ind) - sorted_data(ind-1);
                        d2 = sorted_data(ind+1) - sorted_data(ind);
                        
                        k = d1 < d2;
                        d.distance(tags(ind(k))) = d1(k);
                        d.columnindex(tags(ind(k))) = tags(ind(k)-1);
                        k = ~k;
                        d.distance(tags(ind(k))) = d2(k);
                        d.columnindex(tags(ind(k))) = tags(ind(k)+1);
                    end % if n1 == 2
                else
                    % Two sets of data. still really a sort and some fuss.
                    if n2 == 1
                        % there is only one point in data2
                        d.distance = abs(data1 - data2);
                        % d.columnindex is already set correctly
                    else
                        % At least two points in data2
                        % We need to sort all the data points together, but also
                        % know which points from each set went where. ind12 and
                        % bool12 will help keep track.
                        ind12 = [1:n1,1:n2]';
                        bool12 = [zeros(n1,1);ones(n2,1)];
                        [sorted_data,tags] = sort([data1;data2]);
                        
                        ind12 = ind12(tags);
                        bool12 = bool12(tags);
                        
                        % where did each point end up after the sort?
                        loc1 = find(~bool12);
                        loc2 = find(bool12);
                        
                        % for each point in data1, what is the (sorted) data2
                        % element which appears most nearly to the left of it?
                        cs = cumsum(bool12);
                        leftelement = cs(loc1);
                        
                        % any points which fell below the minimum element in data2
                        % will have a zero for the index of the element on their
                        % left. fix this.
                        leftelement = max(1,leftelement);
                        
                        % likewise, any point greater than the max in data2 will
                        % have an n2 in left element. this too will be a problem
                        % later, so fix it.
                        leftelement = min(n2-1,leftelement);
                        
                        % distance to the left hand element
                        dleft = abs(sorted_data(loc1) - sorted_data(loc2(leftelement)));
                        dright = abs(sorted_data(loc1) - sorted_data(loc2(leftelement+1)));
                        
                        % find the points which are closer to the left element in data2
                        k = (dleft < dright);
                        d.distance(ind12(loc1(k))) = dleft(k);
                        d.columnindex(ind12(loc1(k))) = ind12(loc2(leftelement(k)));
                        k = ~k;
                        d.distance(ind12(loc1(k))) = dright(k);
                        d.columnindex(ind12(loc1(k))) = ind12(loc2(leftelement(k)+1));
                        
                    end % if n2 == 1
                end % if dataflag == 1
            end % if params.Subset(1) == 'f'
            
            % create the matrix as a sparse one or a struct?
            if params.Result(1)=='a'
                % its an array, so make the array sparse.
                d = sparse(d.rowindex,d.columnindex,d.distance,n1,n2);
            end
            
        elseif (ntotal>1000) && (((params.Metric == 0) && (params.Subset(1) == 'n')) || ...
                ((params.Metric == inf) && (params.Subset(1) == 'f')))
            % nearest/farthest neighbour in n>1 dimensions, but for an
            % infinity norm metric. Reduce this to a sequence of
            % 1-d problems, each of which will be faster in general.
            % do this only if the problem is moderately large, since
            % we must overcome the extra overhead of the recursive
            % calls to ipdm.
            
            % do the first dimension
            if dataflag == 1
                d = ipdm(data1(:,1),data1(:,1),'subset',params.Subset,'metric',params.Metric,'result','struct');
            else
                d = ipdm(data1(:,1),data2(:,1),'subset',params.Subset,'metric',params.Metric,'result','struct');
            end
            
            % its slightly different for nearest versus farthest here
            % now, loop over dimensions
            for i = 2:dim
                if dataflag == 1
                    di = ipdm(data1(:,i),data1(:,i),'subset',params.Subset,'metric',params.Metric,'result','struct');
                else
                    di = ipdm(data1(:,i),data2(:,i),'subset',params.Subset,'metric',params.Metric,'result','struct');
                end
                
                % did any of the distances change?
                if params.Metric == 0
                    % the 0 norm, with nearest neighbour, so take the
                    % smallest distance in any dimension.
                    k = d.distance > di.distance;
                else
                    % inf norm. so take the largest distance across dimensions
                    k = d.distance < di.distance;
                end
                
                if any(k)
                    d.distance(k) = di.distance(k);
                    d.columnindex(k) = di.columnindex(k);
                end
            end
            
            % create the matrix as a sparse one or a struct?
            if params.Result(1)=='a'
                % its an array, so make the array sparse.
                d = sparse(d.rowindex,d.columnindex,d.distance,n1,n2);
            end
            
        elseif ((ntotal*8) <= params.ChunkSize)
            % None of the other special cases apply, so do it using brute
            % force for the small potatoes problem.
            
            % One set or two?
            if dataflag == 1
                dist = distcomp(data1,data1,params);
            else
                dist = distcomp(data1,data2,params);
            end
            
            % if only one data set and if a nearest neighbor
            % problem, set the diagonal to +inf so we don't find it.
            if (dataflag==1) && (n1>1) && ('n'==params.Subset(1))
                diagind = (1:n1) + (0:n1:(n1^2-1));
                dist(diagind) = +inf;
            end
            
            if ('n'==params.Subset(1))
                % nearest
                [val,j] = min(dist,[],2);
            else
                % farthest
                [val,j] = max(dist,[],2);
            end
            
            % create the matrix as a sparse one or a struct?
            if params.Result(1)=='a'
                % its an array, so make the array sparse.
                d = sparse((1:n1)',j,val,n1,size(dist,2));
            else
                % a structure
                d.rowindex = (1:n1)';
                d.columnindex = j;
                d.distance = val;
            end
            
        else
            
            % break it into chunks
            bs = floor(params.ChunkSize/(8*n2));
            bs = min(n1,max(1,bs));
            
            % pre-allocate the result
            d.rowindex = (1:n1)';
            d.columnindex = zeros(n1,1);
            d.distance = zeros(n1,1);
            
            % now loop over the chunks
            batch = 1:bs;
            while ~isempty(batch)
                
                % One set or two?
                if dataflag == 1
                    dist = distcomp(data1(batch,:),data1,params);
                else
                    dist = distcomp(data1(batch,:),data2,params);
                end
                
                % if only one data set and if a nearest neighbor
                % problem, set the diagonal to +inf so we don't find it.
                if (dataflag==1) && (n1>1) && ('n'==params.Subset(1))
                    diagind = 1:length(batch);
                    diagind = diagind + (diagind-2+batch(1))*length(batch);
                    dist(diagind) = +inf;
                end
                
                % big or small as requested
                if ('n'==params.Subset(1))
                    % nearest
                    [val,j] = min(dist,[],2);
                else
                    % farthest
                    [val,j] = max(dist,[],2);
                end
                
                % and stuff them into the result structure
                d.columnindex(batch) = j;
                d.distance(batch) = val;
                
                % increment the batch
                batch = batch + bs;
                if batch(end)>n1
                    batch(batch>n1) = [];
                end
                
            end
            
            % did we need to return a struct or an array?
            if params.Result(1) == 'a'
                % an array. make it a sparse one
                d = sparse(d.rowindex,d.columnindex,d.distance,n1,n2);
            end
            
        end % if dim == 1
        
end  % switch params.Subset

% End of mainline
end
% ======================================================
% begin subfunctions
% ======================================================
function d = distcomp(set1,set2,params)
% Subfunction to compute all distances between two sets of points
dim = size(set1,2);
% can we take advantage of bsxfun?
% Note: in theory, there is no need to loop over the dimensions. We
% could Just let bsxfun do ALL the work, then wrap a sum around the
% outside. In practice, this tends to create large intermediate
% arrays, especially in higher numbers of dimensions. Its also when
% we might gain here by use of a vectorized code. This will only be
% a serious gain when the number of points is relatively small and
% the dimension is large.
if params.usebsxfun
    % its a recent enough version of matlab that we can
    % use bsxfun at all.
    n1 = size(set1,1);
    n2 = size(set2,1);
    if (dim>1) && ((n1*n2*dim)<=params.ChunkSize)
        % its a small enough problem that we might gain by full
        % use of bsxfun
        switch params.Metric
            case 2
                d = sum(bsxfun(@minus,reshape(set1,[n1,1,dim]),reshape(set2,[1,n2,dim])).^2,3);
            case 1
                d = sum(abs(bsxfun(@minus,reshape(set1,[n1,1,dim]),reshape(set2,[1,n2,dim]))),3);
            case inf
                d = max(abs(bsxfun(@minus,reshape(set1,[n1,1,dim]),reshape(set2,[1,n2,dim]))),[],3);
            case 0
                d = min(abs(bsxfun(@minus,reshape(set1,[n1,1,dim]),reshape(set2,[1,n2,dim]))),[],3);
        end
    else
        % too big, so that the ChunkSize will have been exceeded, or just 1-d
        if params.Metric == 2
            d = bsxfun(@minus,set1(:,1),set2(:,1)').^2;
        else
            d = abs(bsxfun(@minus,set1(:,1),set2(:,1)'));
        end
        for i=2:dim
            switch params.Metric
                case 2
                    d = d + bsxfun(@minus,set1(:,i),set2(:,i)').^2;
                case 1
                    d = d + abs(bsxfun(@minus,set1(:,i),set2(:,i)'));
                case inf
                    d = max(d,abs(bsxfun(@minus,set1(:,i),set2(:,i)')));
                case 0
                    d = min(d,abs(bsxfun(@minus,set1(:,i),set2(:,i)')));
            end
        end
    end
else
    % Cannot use bsxfun. Sigh. Do things the hard (and slower) way.
    n1 = size(set1,1);
    n2 = size(set2,1);
    if params.Metric == 2
        % Note: While some people might use a different Euclidean
        % norm computation based on expanding the square of the
        % difference of two numbers, that computation is inherantly
        % inaccurate when implemented in floating point arithmetic.
        % While it might be faster, I won't use it here. Sorry.
        d = (repmat(set1(:,1),1,n2) - repmat(set2(:,1)',n1,1)).^2;
    else
        d = abs(repmat(set1(:,1),1,n2) - repmat(set2(:,1)',n1,1));
    end
    for i=2:dim
        switch params.Metric
            case 2
                d = d + (repmat(set1(:,i),1,n2) - repmat(set2(:,i)',n1,1)).^2;
            case 1
                d = d + abs(repmat(set1(:,i),1,n2) - repmat(set2(:,i)',n1,1));
            case inf
                d = max(d,abs(repmat(set1(:,i),1,n2) - repmat(set2(:,i)',n1,1)));
            case 0
                d = min(d,abs(repmat(set1(:,i),1,n2) - repmat(set2(:,i)',n1,1)));
        end
    end
end
% if 2 norm, then we must sqrt at the end
if params.Metric==2
    d = sqrt(d);
end
end
% ==============================================================
%    end main ipdm
% ==============================================================

%%%%%----------------------------------------------------------------------
function skill=modskill(pred,obs,dim)
% MSKILL - Calculate model performance or "skill" measures
%
%   SKILL =  MSKILL(PRED,OBS) calculates model performance measures
%   suggested by Willmott (1982) betweeen predicted or modeled
%   parameter, PRED, and observed value, OBS. NaNs are removed.
%   PRED and OBS should have the same size and shape. If a matrix
%   is provided for PRED and OBS, then the model skill will be
%   computed for each column-wise (unless an optional DIM argument
%   is provided).
%
%   SKILL = MSKILL(...,DIM) - computes the model skill along the
%   dimension DIM.  MSKILL is only valid for 2D matrices, so DIM
%   should either be 1 (default) or 2.
%
%   OUTPUT - Function returns a (possibly n dimensional structure):
%
%       'n'         -  number of elements
%       'meanp'     -  mean of predicted values
%       'stdp'      -  standard deviation of predicted values
%       'meano'     -  mean of observations
%       'stdo'      -  standard deviation of predicted values
%       'slope'     -  slope of least squares linear line
%       'intercept' -  intercept of least squares linear line
%       'rmsa'      -  systematic error
%       'rmsu'      -  unsystematic error
%       'd'         -  index of agreement
%
% REFERENCE
%   Willmot, C.J., (1982) Some comments on the evaluation
%      of model performance, Bulletin of American Meteorological
%      Society, vol. 63, pp 1309-1313.
%
% EXAMPLE
% %Create some data
%  x =0:0.1:2*pi;
%  pred = sin(x);
%  obs = pred+0.2*rand(1,numel(pred));
%
% %calculate the skill of the estimate
%  skill = modskill(pred,obs);


%check inputs
if (~isequal(numel(pred),numel(obs)) || ...
        ~isequal(size(pred),size(obs)));
    error('Inputs should be vectors or matrices of the same size');
end
if ~exist('dim','var')
    dim=1;
else
    if dim==0 || dim>2
        error('Dimension can only be 1 or 2');
    end
end
[m,n,z]=size(obs);
if z>1
    error('MSKILL only works on 2-D matrices');
end

if m==1 || n==1 %vector case
    skill = calc_error(pred,obs);
else %matrix case
    predc=num2cell(pred,dim);
    obsc=num2cell(obs,dim);
    skill=cellfun(@(x,y)(calc_error(x,y)),predc,obsc);
end
end
%engine-------------------------------------------------------------
function skill=calc_error(pred,obs)
%remove nans
idx = isnan(pred) | isnan(obs);
pred(idx)=[];
obs(idx)=[];

%summary measures
skill.n=numel(pred);
skill.meanp=mean(pred);
skill.stdp=std(pred);
skill.meano=mean(obs);
skill.stdo=std(obs);

%linear correlation
p=polyfit(obs,pred,1);
skill.slope=p(1);
skill.intercept=p(2);

%difference measures
phat=polyval(p,obs);
skill.rmse=sqrt((sum((pred-obs).^2)/skill.n)); %root-mean square error
skill.rmsa=sqrt((sum((phat-obs).^2))./skill.n); %systematic error
skill.rmsu=sqrt((sum((pred-phat).^2))./skill.n); %unsystematic error
skillf = @(x,y)(1 -(sum((x-y).^2) / ...
    sum((abs(x-mean(y))+ abs(y-mean(y))).^2)));

skill.d=skillf(pred,obs);
end

%%%------------------------------------------------------------------------
function run_tidedlg(hfig,evnt) %#ok
gdata=guidata(hfig);
gdata.tideopt=tidedlg(gdata.tideopt);
guidata(hfig,gdata);

applytidecorr(hfig);

gdata=guidata(hfig);
set(gdata.l1,'ydata',gdata.bdata.zc);

end

%%%------------------------------------------------------------------------
function applytidecorr(hfig,evnt) %#ok
gdata=guidata(hfig);

if strcmpi(gdata.tideopt.method,'none')
    %     gdata.bdata.tide(gdata.tideopt.badtide)=NaN;
    %     gdata.bdata.zc=((gdata.bdata.zraw-...
    %         gdata.bdata.tide).*...
    %         gdata.invert)-gdata.manoff;
    %     if isfield(gdata,'l1')
    %         if ishandle(gdata.l1)
    %             set(gdata.l1,'ydata',gdata.bdata.zc);
    %         end
    %     end
    %     guidata(hfig,gdata);
    return
    
end

%fill gaps less than max value
if isempty(gdata.tideopt.badtide);
    gaps=isnan(gdata.bdata.tide(1:end-1));
    gdata.tideopt.badtide=gaps;
else
    
    gaps=gdata.tideopt.badtide;
    gdata.bdata.tide(gaps)=NaN;
    gdata.bdata.zc=((gdata.bdata.zraw-...
        gdata.bdata.tide).*...
        gdata.invert)-gdata.manoff;
end

if isempty(gaps)~=1;
    [c,ind]=getchunks(gaps,'-full');
    bind=find(isnan(gdata.bdata.tide(ind))==1);
    bstart=ind(bind)-1;
    
    bend=bsxfun(@plus,bstart,c(bind)+1);
    
    
    bend(bstart==0)=[];
    bstart(bstart==0)=[];
    bstart(bend==numel(gdata.bdata.distance))=[];
    bend(bend==numel(gdata.bdata.distance))=[];
    
    dist=cellfun(@(x,y)(gdata.bdata.distance(y)-...
        gdata.bdata.distance(x)),...
        num2cell(bstart),num2cell(bend));
    
    
    fillGaps=find(dist<gdata.tideopt.maxgap)  ;
    fillStart=bstart(fillGaps);
    fillEnd=bend(fillGaps);
    
    switch gdata.tideopt.method
        case 'linear'
            for i=1:numel(fillGaps)
                
                gdata.bdata.tide(fillStart(i):fillEnd(i))=...
                    interp1([gdata.bdata.distance(fillStart(i));...
                    gdata.bdata.distance(fillEnd(i))],...
                    [gdata.bdata.tide(fillStart(i));...
                    gdata.bdata.tide(fillEnd(i))],...
                    gdata.bdata.distance(fillStart(i):fillEnd(i)));
                
                gdata.bdata.zc(fillStart(i):fillEnd(i))=...
                    ((gdata.bdata.zraw(fillStart(i):fillEnd(i))-...
                    gdata.bdata.tide(fillStart(i):fillEnd(i))).*...
                    gdata.invert)-gdata.manoff;
            end
        case 'mean'
            
            
            for i=1:numel(fillGaps)
                ind=isfinite(gdata.bdata.tide);
                p=polyfit(gdata.bdata.distance(ind)',...
                    gdata.bdata.tide(ind),1);
                
                gdata.bdata.tide(fillStart(i):fillEnd(i))=...
                    polyval(p,gdata.bdata.distance(fillStart(i):fillEnd(i)));
                gdata.bdata.zc(fillStart(i):fillEnd(i))=...
                    ((gdata.bdata.zraw(fillStart(i):fillEnd(i))-...
                    gdata.bdata.tide(fillStart(i):fillEnd(i))).*...
                    gdata.invert)-gdata.manoff;
            end
    end
    
end

guidata(hfig,gdata);
end

%%%------------------------------------------------------------------------
function [tide_opt] = tidedlg(varargin)

if nargin>0
    gd.tide=varargin{1};
else
    gd.tide.maxgap=20;
end

hf = figure('units','normalized','position',[0.253 0.434 0.167 0.191],...
    'menubar','none','name','Tide Options','numbertitle','off',...
    'color',[0.925 0.914 0.847]);

gd.checkbox1 = uicontrol(hf,'style','checkbox','units','normalized',...
    'position',[0.0962 0.791 0.647 0.115],...
    'string','Fill Missing Tide Corrections',...
    'backgroundcolor',[0.925 0.914 0.847],...
    'callback',@fillmotion);



uipanel1 = uipanel('parent',hf,'units','normalized',...
    'position',[0.0845 0.225 0.787 0.516],'title','Options');

gd.checkbox2 = uicontrol(uipanel1,'style','checkbox','units','normalized',...
    'position',[0.0714 0.413 0.88 0.229],...
    'string','Use Mean Tide','backgroundcolor',[0.925 0.914 0.847],...
    'enable','off','callback',@meanmotion);

gd.checkbox3 = uicontrol(uipanel1,'style','checkbox','units','normalized',...
    'position',[0.0714 0.679 0.88 0.229],...
    'string','Linear Interpolation',...
    'backgroundcolor',[0.925 0.914 0.847],'enable','off',...
    'callback',@linmotion);

gd.text1 = uicontrol(uipanel1,'style','text','units','normalized',...
    'position',[0.447 0.101 0.496 0.239],...
    'string','Max. Gap (m)','backgroundcolor',[0.925 0.914 0.847],...
    'horizontalalign','left','enable','off');
gd.edit1 = uicontrol(uipanel1,'style','edit','units','normalized',...
    'position',[0.0714 0.0917 0.312 0.248],...
    'string',sprintf('%0.1f',gd.tide.maxgap),'backgroundcolor',[1 1 1],...
    'enable','off','callback',@editmotion);


gd.push1 = uicontrol(hf,'style','pushbutton','units','normalized',...
    'position',[0.373 0.0574 0.207 0.119],'string','Done',...
    'backgroundcolor',[0.925 0.914 0.847],'callback',@closetide);

idx=strcmpi(gd.tide.method,{'linear';'mean'});
if any(idx)
    if idx(1)
        set(gd.checkbox3,'value',1);
    else
        set(gd.checkbox2,'value',1)
    end
    set(gd.checkbox1,'value',1);
    set(gd.checkbox2,'enable','on')
    set(gd.checkbox3,'enable','on')
    set(gd.text1,'enable','on')
    set(gd.edit1,'enable','on')
end


guidata(hf,gd);

uiwait
gd=guidata(hf);
tide_opt=gd.tide;
close(hf)
end
%%%------------------------------------------------------------------------
function fillmotion(hf,evnt) %#ok

gd=guidata(hf);
val=get(gd.checkbox1,'value');
switch val
    case 1
        set(gd.checkbox2,'enable','on')
        set(gd.checkbox3,'enable','on')
        set(gd.text1,'enable','on')
        set(gd.edit1,'enable','on')
        set(gd.push1,'enable','off')
    case 0
        set(gd.checkbox2,'enable','off')
        set(gd.checkbox3,'enable','off')
        set(gd.text1,'enable','off')
        set(gd.edit1,'enable','off')
        set(gd.push1,'enable','on')
end
end
%%%------------------------------------------------------------------------

function meanmotion(hf,evnt) %#ok
gd=guidata(hf);
val=get(gd.checkbox2,'value');
switch val
    case 1
        set(gd.checkbox3,'enable','off')
        set(gd.push1,'enable','on')
        
    case 0
        set(gd.checkbox3,'enable','on')
        set(gd.push1,'enable','off')
        
end
end

%%%------------------------------------------------------------------------
function linmotion(hf,evnt) %#ok
gd=guidata(hf);
val=get(gd.checkbox3,'value');
switch val
    case 1
        set(gd.checkbox2,'enable','off')
        set(gd.push1,'enable','on')
        
    case 0
        set(gd.checkbox2,'enable','on')
        set(gd.push1,'enable','off')
        
end
end
%%%------------------------------------------------------------------------
function editmotion(hf,evnt) %#ok
gd=guidata(hf);
val=str2double(get(gd.edit1,'string'));

if isnan(val) || val<0
    errordlg('Value should be a positive value.')
    set(gd.edit1,'string',sprintf('%0.1f',gd.maxgap))
else
    gd.maxgap=val;
    set(gd.edit1,'string',sprintf('%0.1f',gd.maxgap))
    guidata(hf,gd)
end
end
%%%------------------------------------------------------------------------
function closetide(hf,evnt) %#ok

gd=guidata(hf);

val0=get(gd.checkbox1,'value');
if val0==1
    
    val=get(gd.checkbox3,'value');
    switch val
        case 1
            gd.tide.method='linear';
        case 0
            gd.tide.method='mean';
    end
    gd.tide.maxgap=str2double(get(gd.edit1,'string'));
else
    gd.tide.method='none';
    gd.tide.maxgap=gd.tide.maxgap;
end

guidata(hf,gd)
uiresume

end

%%%%----------------------------------------------------------------------
function res = uical(sDate,lang)
% UICAL - Calendar date picker
%
%   Version : 1.1
%   Created : 10/08/2007
%   Modified: 14/08/2007
%   Author  : Thomas Montagnon (The MathWorks France)
%
%   RES = UICAL() displays the calendar in English with the today's date
%   selected. RES is the serial date number corresponding to the date the user
%   selected.
%
%   RES = UICAL(SDATE) displays the calendar in English with SDATE date
%   selected. SDATE can be either a serial date number or a 3-elements vector
%   containing year, month and day in this order.
%
%   RES = UICAL(SDATE,LANG) displays the calendar in the language specified by
%   LANG with SDATE date selected.
%   Default available languages are English and French but you can easily add
%   your own by editing the SET_LANGUAGES nested function at the bottom of this
%   file.
%
%
%   Customization Tips:
%   You can customize the colors of the calendar by editing the "Colors" section
%   at the begining of the main function. See comments there for more details.
%
%
%   Examples:
%
%     myDate = uical()
%
%     myDate = uical(now+7)
%
%     myDate = uical([1980 7 24],'fr')
%
%   See also CALENDAR.
%


% ------------------------------------------------------------------------------
% --- INITIALIZATIONS                                                        ---
% ------------------------------------------------------------------------------

% Colors
figColor     = get(0,'DefaultUicontrolBackgroundColor'); % Figure background color
colorNoDay   = [0.95  0.95  0.95]; % Background color of the cells that are not days of the selected month
colorDayB    = [1.00  1.00  1.00]; % Background color of the cells that are day of the selected month
colorDayF    = [0.00  0.00  0.00]; % Foreground color of the cells that are day of the selected month
colorDayNB   = [0.30  0.30  0.30]; % Background color of the column headers
colorDayNF   = [1.00  1.00  1.00]; % Foreground color of the column headers
colorSelDayB = [0.70  0.00  0.00]; % Background color of the selected day
colorSelDayF = [1.00  1.00  1.00]; % Foreground color of the selected day

% Default input arguments
sDateDef = now;
langDef  = 'en';

% Use default values if input does not exist
switch nargin
    case 0
        sDate = sDateDef;
        lang  = langDef;
    case 1
        lang = langDef;
end

% Check input argument validity
if ~isnumeric(sDate)
    error('MYCAL:DateFormat:WrongClass','First input must be numeric');
end
switch numel(sDate)
    case 1
        sDate = datevec(sDate);
    case 3
        if sDate(1) < 0
            error('MYCAL:DateFormat:WrongYearVal','First element of the first input must be a valid year number');
        end
        if (sDate(2) > 12) && (sDate(2) < 1)
            error('MYCAL:DateFormat:WrongMonthVal','Second element of the first input must be a valid month number');
        end
        if (sDate(3) > 31) && (sDate(3) < 1)
            error('MYCAL:DateFormat:WrongDayVal','Third element of the first input must be a valid day number');
        end
    otherwise
        error('MYCAL:DateFormat:WrongVal','First input must be a numeric scalar or a 3-elements vector');
end


% Language dependent strings
set_language(lang);

% Initializes output
res = sDate(1:3);

% Dimensions
dayH   = 24;
dayW   = 30;
ctrlH  = 20;
figH   = 7 * (dayH - 1) + ctrlH;
figW   = 7 * (dayW-1);
ctrlYW = 60;
ctrlMW = 80;
ctrlCW = figW - ctrlYW - ctrlMW;

daysNy = figH - ctrlH - dayH + 1;


% ------------------------------------------------------------------------------
% --- UICONTROL CREATION                                                     ---
% ------------------------------------------------------------------------------


% Create figure
handles.FgCal = figure( ...
    'Visible', 'off', ...
    'Tag', 'FgCal', ...
    'Name', '', ...
    'Units', 'pixels', ...
    'Position', [50 50 figW figH], ...
    'Toolbar', 'none', ...
    'MenuBar', 'none', ...
    'NumberTitle', 'off', ...
    'Color', figColor, ...
    'CloseRequestFcn',@FgCal_CloseRequestFcn,...
    'WindowStyle','modal');

% Move the GUI to the center of the screen
movegui(handles.FgCal,'center')

% Columns Headers containing initials of the week days
for dayNidx=1:7
    
    daysNx = (dayNidx - 1) * (dayW - 1);
    
    handles.EdDayN(dayNidx) = uicontrol( ...
        'Parent', handles.FgCal, ...
        'Tag', 'EdDay', ...
        'Style', 'edit', ...
        'Units', 'pixels', ...
        'Position', [daysNx daysNy dayW dayH], ...
        'ForegroundColor', colorDayNF, ...
        'BackgroundColor', colorDayNB, ...
        'String', daysN{dayNidx}, ...
        'HorizontalAlignment', 'center', ...
        'Enable','inactive');
    
end

% Days UI controls
for dayIdx=1:42
    
    % X and Y Positions
    [i,j] = ind2sub([6,7],dayIdx);
    
    dayX = (j - 1) * (dayW - 1);
    dayY = (dayH - 1) * 6 - i * (dayH - 1);
    
    handles.EdDay(dayIdx) = uicontrol( ...
        'Parent', handles.FgCal, ...
        'Tag', 'EdDay', ...
        'Style', 'edit', ...
        'Units', 'pixels', ...
        'Position', [dayX dayY dayW dayH], ...
        'BackgroundColor', colorDayB, ...
        'ForegroundColor', colorDayF, ...
        'String', '', ...
        'HorizontalAlignment', 'center', ...
        'Enable','inactive');
    
end

% Listbox containing the list of months
handles.PuMonth = uicontrol( ...
    'Parent', handles.FgCal, ...
    'Tag', 'PuMonth', ...
    'Style', 'popupmenu', ...
    'Units', 'pixels', ...
    'Position', [ctrlYW-2 figH-ctrlH+1 ctrlMW+2 ctrlH], ...
    'BackgroundColor', [1 1 1], ...
    'String', monthsN, ...
    'Value', res(2), ...
    'Callback',@set_cal);

% Edit control which enables you to enter a year number
handles.EdYear = uicontrol( ...
    'Parent', handles.FgCal, ...
    'Tag', 'EdYear', ...
    'Style', 'edit', ...
    'Units', 'pixels', ...
    'Position', [0 figH-ctrlH ctrlYW-1 ctrlH+1], ...
    'BackgroundColor', [1 1 1], ...
    'String', res(1), ...
    'Callback',@set_cal);

% Selection button
handles.PbChoose = uicontrol( ...
    'Parent', handles.FgCal, ...
    'Tag', 'PbChoose', ...
    'Style', 'pushbutton', ...
    'Units', 'pixels', ...
    'Position', [ctrlYW+ctrlMW figH-ctrlH+1 ctrlCW ctrlH], ...
    'String', 'OK', ...
    'Callback','uiresume');

% Display calendar for the default date
set_cal();

% Make the calendar visible
set(handles.FgCal,'Visible','on')

% Wait for user action
uiwait(handles.FgCal);

% Convert date to serial date number
res = datenum(res);

% Close the calendar figure
delete(handles.FgCal);


% ------------------------------------------------------------------------------
% --- CALLBACKS                                                              ---
% ------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
    function FgCal_CloseRequestFcn(varargin)
        % Callback executed when the user click on the close button of the figure.
        % This means he wants to cancel date selection so function returns the
        % intial date (the one used when we opened the calendar)
        
        % Set the output to the intial date value
        res = sDate(1:3);
        
        % End execution of the window
        uiresume;
        
    end


%-------------------------------------------------------------------------------
    function EdDay_ButtonDownFcn(varargin)
        % Callback executed when the user click on day.
        % Updates the RES variable containing the currently selected date and then
        % update the calendar.
        
        res(1) = str2double(get(handles.EdYear,'String'));
        res(2) = get(handles.PuMonth,'Value');
        res(3) = str2double(get(varargin{1},'String')); % Number of the selected day
        
        set_cal();
        
    end


%-------------------------------------------------------------------------------
    function set_cal(varargin)
        % Displays the calendar according to the selected date stored in RES
        
        % Get selected Year and Month
        year   = str2double(get(handles.EdYear,'String'));
        res(2) = get(handles.PuMonth,'value');
        
        % Check Year value (keep previous value if the new one is wrong)
        if ~isnan(year)
            res(1) = abs(round(year)); % ensure year is a positive integer
        end
        set(handles.EdYear,'String',res(1))
        
        % Get the matrix of the calendar for selected month and year then convert it
        % into a cell array
        c = calendar(res(1),res(2));
        v = mat2cell(c,ones(1,6),ones(1,7));
        
        % Cell array of indices used to index the vector of handles
        i = mat2cell((1:42)',ones(1,42),1);
        
        % Set String property for all cells of the calendar
        cellfun(@(i,x) set(handles.EdDay(i),'string',x),i,v(:))
        
        % Change properties of the "non-day" cells of the calendar
        set(handles.EdDay(c==0), ...
            'ButtonDownFcn'  , '', ...
            'BackgroundColor', colorNoDay, ...
            'string'         , '')
        
        % Change the properties of the calendar's cells containing existing days
        set(handles.EdDay(c~=0), ...
            'ButtonDownFcn'  , @EdDay_ButtonDownFcn, ...
            'BackgroundColor', colorDayB, ...
            'ForegroundColor', colorDayF, ...
            'FontWeight'     ,'normal')
        
        % Highlight the selected day
        set(handles.EdDay(c==res(3)), ...
            'BackgroundColor', colorSelDayB, ...
            'ForegroundColor', colorSelDayF, ...
            'FontWeight'     ,'bold')
        
        % Update the name of the figure to reflect the selected day
        set(handles.FgCal,'Name',sprintf('%u/%u/%u',fliplr(res)))
        
        % Give focus to the "OK" button
        uicontrol(handles.PbChoose);
        
    end


%-------------------------------------------------------------------------------
    function set_language(lang)
        % Sets language dependent strings used in the calendar
        % You can add languages by adding cases below.
        
        switch lang
            
            case 'en'
                
                daysN   = {'S','M','T','W','T','F','S'}; % First day is always Sunday
                monthsN = {'January','February','March','April','May','June',...
                    'July','August','September','October','November','December'};
                
            case 'fr'
                
                daysN   = {'D','L','M','M','J','V','S'};
                monthsN = {'Janvier','Février','Mars','Avril','Mai','Juin',...
                    'Juillet','Août','Septembre','Octobre','Novembre','Décembre'};
                
            otherwise
                
                % If language is not recognized then use the English strings
                
                daysN   = {'S','M','T','W','T','F','S'};
                monthsN = {'January','February','March','April','May','June',...
                    'July','August','September','October','November','December'};
                
        end
        
    end

end



