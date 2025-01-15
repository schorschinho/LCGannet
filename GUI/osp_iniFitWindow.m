function osp_iniFitWindow(gui)
%% osp_iniFitWindow
%   This function creates the initial fitting window in the gui.
%
%   USAGE:
%       osp_iniFitWindow(gui);
%
%   INPUT:      gui      = gui class containing all handles and the MRSCont
%
%   OUTPUT:     Changes in gui parameters and MRSCont are written into the
%               gui class
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-01-16)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-01-16: First version of the code.


%%% 1. GET HANDLES %%%
% This function initializes the fit tab
MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class
gui.layout.tabs.TabEnables{3} = 'on';
gui.layout.tabs.Selection  = 3;
gui.layout.EmptyFitPlot = 0;

%%% 2. CREATING SUB TABS FOR THIS TAB %%%%
% In this case one tab fo each fit (off,sum,diff1,diff2,ref,water)
gui.layout.fitTab.TabWidth   = 115;
for t = 1 : gui.fit.Number %Create tabs depending on the number of fits
    gui.layout.(['fitTab' gui.fit.Names{t}]) = uix.VBox('Parent', gui.layout.fitTab,...
        'BackgroundColor',gui.colormap.Background,'Spacing',5);
    gui.layout.fitTabhandles{t} = ['fitTab' gui.fit.Names{t}];
end
gui.layout.fitTab.TabTitles  = gui.fit.Names;

%%% 3. FILLING INFO PANEL FOR THIS TAB %%%%
% All the information from the Raw data is read out here
for t = 1 : gui.fit.Number %Loop over fits
    Selection = gui.fit.Names{t};
    % Parameter shown in the info panel on top
    gui.upperBox.fit.box{t} = uix.HBox('Parent', gui.layout.(gui.layout.fitTabhandles{t}),'BackgroundColor',gui.colormap.Background,'Spacing',5);
    if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
        gui.upperBox.fit.upperLeftButtons = uix.Panel('Parent', gui.upperBox.fit.box{t}, ...
            'Padding', 5, 'Title', ['Navigate voxel'],...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
            'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
        gui.controls.Buttonbox = uix.HBox('Parent',gui.upperBox.fit.upperLeftButtons, 'BackgroundColor',gui.colormap.Background);
        gui.controls.navigate_RawTab = uix.Grid('Parent',gui.controls.Buttonbox,'BackgroundColor',gui.colormap.Background);
        gui.controls.text_x = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','X:',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.controls.text_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Y:',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.controls.text_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Z:',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.controls.b_left_x = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
        gui.controls.b_left_y = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
        gui.controls.b_left_z = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
        set(gui.controls.b_left_x,'Callback',{@osp_onLeftX,gui});
        set(gui.controls.b_left_y,'Callback',{@osp_onLeftY,gui});
        set(gui.controls.b_left_z,'Callback',{@osp_onLeftZ,gui});
        if gui.info.nXvoxels <= 1
            gui.controls.b_left_x.Enable = 'off';
        end
        if gui.info.nYvoxels <= 1
            gui.controls.b_left_y.Enable = 'off';
        end
        if gui.info.nZvoxels <= 1
            gui.controls.b_left_z.Enable = 'off';
        end
        gui.controls.text_act_x = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.controls.text_act_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.controls.text_act_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.controls.b_right_x = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
        gui.controls.b_right_y = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
        gui.controls.b_right_z = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
        set(gui.controls.b_right_x,'Callback',{@osp_onRightX,gui});
        set(gui.controls.b_right_y,'Callback',{@osp_onRightY,gui});
        set(gui.controls.b_right_z,'Callback',{@osp_onRightZ,gui});
        if gui.info.nXvoxels <= 1
            gui.controls.b_right_x.Enable = 'off';
        end
        if gui.info.nYvoxels <= 1
            gui.controls.b_right_y.Enable = 'off';
        end
        if gui.info.nZvoxels <= 1
            gui.controls.b_right_z.Enable = 'off';
        end
        set( gui.controls.navigate_RawTab, 'Widths', [-20 -30 -20 -30], 'Heights', [-33 -33 -33] );
    end


    gui.upperBox.fit.upperLeftButtons = uix.Panel('Parent', gui.upperBox.fit.box{t}, ...
        'Padding', 5, 'Title', ['Navigate model'],...
        'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
        'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
    gui.controls.Buttonbox = uix.HBox('Parent',gui.upperBox.fit.upperLeftButtons, 'BackgroundColor',gui.colormap.Background);
    gui.controls.navigate_RawTab = uix.Grid('Parent',gui.controls.Buttonbox,'BackgroundColor',gui.colormap.Background);

    gui.controls.text_x = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Exp:',...
        'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
    gui.controls.text_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Spec:',...
        'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
    gui.controls.text_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Basis:',...
        'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
    gui.controls.b_left_x = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
    gui.controls.b_left_y = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
    gui.controls.b_left_z = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
    set(gui.controls.b_left_x,'Callback',{@osp_onLeftX,gui});
    set(gui.controls.b_left_y,'Callback',{@osp_onLeftY,gui});
    set(gui.controls.b_left_z,'Callback',{@osp_onLeftZ,gui});


    gui.controls.text_act_x = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
        'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
    gui.controls.text_act_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
        'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
    gui.controls.text_act_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
        'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

    gui.controls.b_right_x = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
    gui.controls.b_right_y = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
    gui.controls.b_right_z = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
    set(gui.controls.b_right_x,'Callback',{@osp_onRightX,gui});
    set(gui.controls.b_right_y,'Callback',{@osp_onRightY,gui});
    set(gui.controls.b_right_z,'Callback',{@osp_onRightZ,gui});

    gui.controls.b_left_x.Enable = 'off';
    gui.controls.b_left_y.Enable = 'off';
    gui.controls.b_left_z.Enable = 'off';
    gui.controls.b_right_x.Enable = 'off';
    gui.controls.b_right_y.Enable = 'off';
    gui.controls.b_right_z.Enable = 'off';

    buttonString = [num2str(MRSCont.nDatasets(2) > 1) num2str(size(MRSCont.fit.results.(Selection).fitParams,3)>1) num2str(size(MRSCont.fit.results.(Selection).fitParams,1)>1)];
    switch buttonString
            case '001'
                gui.controls.b_left_z.Enable = 'on';
                gui.controls.b_right_z.Enable = 'on';
            case '010'
                gui.controls.b_left_y.Enable = 'on';
                gui.controls.b_right_y.Enable = 'on';
            case '100'
                gui.controls.b_left_x.Enable = 'on';
                gui.controls.b_right_x.Enable = 'on';
            case '011'
                gui.controls.b_left_y.Enable = 'on';
                gui.controls.b_right_y.Enable = 'on';
                gui.controls.b_left_z.Enable = 'on';
                gui.controls.b_right_z.Enable = 'on';
           case '101'
                gui.controls.b_left_x.Enable = 'on';
                gui.controls.b_right_x.Enable = 'on';
                gui.controls.b_left_z.Enable = 'on';
                gui.controls.b_right_z.Enable = 'on';
           case '110'
                gui.controls.b_left_x.Enable = 'on';
                gui.controls.b_right_x.Enable = 'on';
                gui.controls.b_left_y.Enable = 'on';
                gui.controls.b_right_y.Enable = 'on';
            case '111'
                gui.controls.b_left_x.Enable = 'on';
                gui.controls.b_left_y.Enable = 'on';
                gui.controls.b_left_z.Enable = 'on';
                gui.controls.b_right_x.Enable = 'on';
                gui.controls.b_right_y.Enable = 'on';
                gui.controls.b_right_z.Enable = 'on';
    end
    set( gui.controls.navigate_RawTab, 'Widths', [-20 -30 -20 -30], 'Heights', [-33 -33 -33]);

    gui.upperBox.fit.Info{t} = uix.Panel('Parent',  gui.upperBox.fit.box{t}, ...
        'Padding', 5, 'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}],...
        'FontName', gui.font,'HighlightColor', gui.colormap.Foreground,'BackgroundColor',...
        gui.colormap.Background,'ForegroundColor',gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);

    gui.upperBox.fit.Cosmet{t} = uix.Panel('Parent',  gui.upperBox.fit.box{t}, ...
        'Padding', 5, 'Title', 'Cosmetics',...
        'FontName', gui.font,'HighlightColor', gui.colormap.Foreground,'BackgroundColor',...
        gui.colormap.Background,'ForegroundColor',gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);

    gui.controls.Cosmetic_Grid = uix.Grid('Parent',gui.upperBox.fit.Cosmet{t},'BackgroundColor',gui.colormap.Background);
    gui.controls.text_zoom = uicontrol(gui.controls.Cosmetic_Grid,'Style','text','String','Zoom',...
        'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

    gui.controls.text_lb = uicontrol(gui.controls.Cosmetic_Grid,'Style','text','String','LB (Hz)',...
        'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
    gui.controls.s_zoom = uicontrol('Parent',gui.controls.Cosmetic_Grid,'Style','Slider');
    gui.controls.s_lb = uicontrol('Parent',gui.controls.Cosmetic_Grid,'Style','Slider');
    set(gui.controls.s_zoom,'Min', 1, 'Max', 3, 'Value', 2.75,'Tooltip', 'Zoom factor', 'SliderStep', [1/16,1]);
    set(gui.controls.s_lb,'Min', 0, 'Max', 8, 'Value', 0, 'Tooltip', ['Exponential Linebroadening (Hz) ' num2str(gui.controls.s_lb.Value)], 'SliderStep', [1/16,1]);
    set(gui.controls.s_zoom,'Callback',{@osp_onZoom,gui});
    set(gui.controls.s_lb,'Callback',{@osp_onLB,gui});
    set( gui.controls.Cosmetic_Grid, 'Widths', [-0.3 -0.5]);

    gui.upperBox.fit.upperButtons = uix.Panel('Parent', gui.upperBox.fit.box{t}, ...
        'Padding', 5, 'Title', ['Save'],...
        'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
        'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
    gui.controls.b_save_fitTab{t} = uicontrol('Parent',gui.upperBox.fit.upperButtons,'Style','PushButton');
    [img, ~, ~] = imread('Printer.png', 'BackgroundColor', gui.colormap.Background);
    [img2] = imresize(img, 0.1);
    set(gui.controls.b_save_fitTab{t},'CData', img2, 'TooltipString', 'Create EPS figure from current file');
    set(gui.controls.b_save_fitTab{t},'Callback',{@osp_onPrint,gui});
    if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
        set(gui.upperBox.fit.box{t}, 'Width', [-0.12 -0.78 -0.1]);
    else
            set(gui.upperBox.fit.box{t}, 'Width', [-0.15 -0.65 -0.1 -0.1]);
    end
    % Creates layout for plotting and data control
    gui.Plot.fit{t} = uix.HBox('Parent', gui.layout.(gui.layout.fitTabhandles{t}), ...
        'Padding', 5,'BackgroundColor',gui.colormap.Background);
    set(gui.layout.(gui.layout.fitTabhandles{t}), 'Heights', [-0.1 -0.9]);
    if  ~strcmp (MRSCont.opts.fit.style, 'Concatenated') ||  strcmp(gui.fit.Names{t}, 'ref') || strcmp(gui.fit.Names{t}, 'w') %Is not concateneted or is reference/water fit
        gui.fit.Style = gui.fit.Names{t};
    else %Is concatenated and not water/reference
        gui.fit.Style = 'conc';
    end


    if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
         switch MRSCont.opts.fit.method
                case 'LCModel'
                    if strcmp(gui.fit.Names{t}, 'ref') || strcmp(gui.fit.Names{t}, 'w')
                        RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.h2oarea .* MRSCont.fit.scale{1,gui.controls.Selected};
                    else
                        RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{1,gui.controls.Selected};
                        CRLB    = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.CRLB;
                    end
                case 'Osprey'
                    if strcmp(gui.fit.Names{t}, 'ref') || strcmp(gui.fit.Names{t}, 'w')
                        RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,end}.ampl .* MRSCont.fit.scale{1,gui.controls.Selected};
                    else
                        RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,end}.ampl .* MRSCont.fit.scale{1,gui.controls.Selected};
                        CRLB    = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,end}.CRLB;
                    end
        end
        ph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,end}.ph0;
        ph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,end}.ph1;
        if ~strcmp(gui.fit.Names{t}, 'ref') && ~strcmp(gui.fit.Names{t}, 'w')
            refShift = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,end}.refShift;
            refFWHM = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,end}.refFWHM;
        end
    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
         switch MRSCont.opts.fit.method
                case 'LCModel'
                    if strcmp(gui.fit.Names{t}, 'ref') || strcmp(gui.fit.Names{t}, 'w')
                        RawAmpl = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{gui.controls.Selected}.h2oarea .* MRSCont.fit.scale{gui.controls.Selected};
                    else
                        RawAmpl = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
                        CRLB    = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{gui.controls.Selected}.CRLB;
                    end
                case 'Osprey'
                    if strcmp(gui.fit.Names{t}, 'ref') || strcmp(gui.fit.Names{t}, 'w')
                        RawAmpl = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
                    else
                        RawAmpl = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
                        CRLB    = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{gui.controls.Selected}.CRLB;
                    end
        end
        ph0 = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph0;
        ph1 = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph1;
        if ~strcmp(gui.fit.Names{t}, 'ref') && ~strcmp(gui.fit.Names{t}, 'w')
            refShift = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refShift;
            refFWHM = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refFWHM;
        end
    else
        switch MRSCont.opts.fit.method
                case 'LCModel'
                    if strcmp(gui.fit.Names{t}, 'ref') || strcmp(gui.fit.Names{t}, 'w')
                        RawAmpl = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{gui.controls.Selected}.h2oarea .* MRSCont.fit.scale{gui.controls.Selected};
                    else
                        RawAmpl = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
                        CRLB    = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{gui.controls.Selected}.CRLB;
                    end
                case 'Osprey'
                    if strcmp(gui.fit.Names{t}, 'ref') || strcmp(gui.fit.Names{t}, 'w')
                        RawAmpl = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
                    else
                        RawAmpl = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
                        CRLB    = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{gui.controls.Selected}.CRLB;
                    end
        end
        ph0 = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph0;
        ph1 = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph1;
        if ~strcmp(gui.fit.Names{t}, 'ref') && ~strcmp(gui.fit.Names{t}, 'w')
            refShift = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refShift;
            refFWHM = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refFWHM;
        end
    end

    % For this visualization, we will have to make a few
    % distinctions upfront since the modeling algorithms (LCModel
    % vs. Osprey) do not always return the same kinds of data, or they
    % return them in different formats.
    switch MRSCont.opts.fit.method
        case 'LCModel'
            % Number of metabolites and lipid/MM basis functions
            basisNames = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.name;
            nLip    = sum(~cellfun(@isempty, strfind(basisNames, 'Lip')));
            nMM     = sum(~cellfun(@isempty, strfind(basisNames, 'MM')));
            nMMLip  = nLip + nMM;
            nMets   = length(basisNames) - nMMLip;
            nComb   = sum(~cellfun(@isempty, strfind(basisNames, '_')));
            % No info panel string for the water fit range
            waterFitRangeString = '';
            % Where are the metabolite names stored?
            basisSetNames = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.name;
            % Smaller fonts for the results
            resultsFontSize = 8;
        case 'Osprey'
            % Number of metabolites and lipid/MM basis functions
            nMets   = MRSCont.fit.basisSet.nMets;
            nMMLip  = MRSCont.fit.basisSet.nMM;
            % Additional info panel string for the water fit range
            waterFitRangeString = ['Fitting range: ' num2str(MRSCont.opts.fit.rangeWater(1)) ' to ' num2str(MRSCont.opts.fit.rangeWater(2)) ' ppm'];
            % Where are the metabolite names stored?
            if strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')
                basisSetNames = MRSCont.fit.resBasisSet.(gui.fit.Style).(['np_sw_' num2str(round(MRSCont.processed.(gui.fit.Style){1}.sz(1))) '_' num2str(round(MRSCont.processed.(gui.fit.Style){1}.spectralwidth))]){1,1}.name;
            else if strcmp(gui.fit.Style, 'conc')
                    basisSetNames = MRSCont.fit.resBasisSet.(gui.fit.Style).(['np_sw_' num2str(round(MRSCont.processed.metab{1}.sz(1))) '_' num2str(round(MRSCont.processed.metab{1}.spectralwidth))]){1,1}.name;
                else
                    basisSetNames = MRSCont.fit.resBasisSet.(gui.fit.Style).(['np_sw_' num2str(round(MRSCont.processed.metab{1}.sz(1))) '_' num2str(round(MRSCont.processed.metab{1}.spectralwidth))]){1,end}.name;
                end
            end
             % Larger fonts for the results
            resultsFontSize = 8;
    end


    % Get parameter from file to fill the info panel
    if  ~strcmp (Selection, 'ref') && ~strcmp (Selection, 'w') %Metabolite data?
        StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' gui.fit.Names{t},...
            '\nFitting range: ' num2str(MRSCont.opts.fit.range(1)) ' to ' num2str(MRSCont.opts.fit.range(2)) ' ppm; Baseline knot spacing: ' num2str(MRSCont.opts.fit.bLineKnotSpace) ' ppm; ph0: ' num2str(ph0,'%1.2f'),...
            'deg; ph1: ' num2str(ph1,'%1.2f') 'deg; refShift: ' num2str(refShift,'%1.2f') ' Hz; refFWHM: ' num2str(refFWHM,'%1.2f')...
            ' ppm\nNumber of metabolites: ' num2str(nMets) '; Number of MM/lipids: ' num2str(nMMLip) ...
            ' scale: '  num2str(MRSCont.fit.scale{gui.controls.Selected})];
    else if strcmp (Selection, 'ref') %Reference data?
            StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' gui.fit.Names{t},...
                '\n' waterFitRangeString];
        else %Is water data
            StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' gui.fit.Names{t},...
                '\n' waterFitRangeString];
        end
    end

    %%% 4. FILLING FITTED AMPLITUDE PANEL %%%
    % Creates the panel on the right side with the fitted ammplitudes
    gui.InfoText.fit{t}  = uicontrol('Parent',gui.upperBox.fit.Info{t},'style','text',...
        'FontSize', 12, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
    gui.Results.fit{t} = uix.Panel('Parent', gui.Plot.fit{t},...
        'Title', ['Raw Amplitudes'],'FontName', gui.font,'HighlightColor', gui.colormap.Foreground,...
        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);

    if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
        if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
            if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                NameText    = [''];
                RawAmplText = [''];
                CRLBText    = [''];
                for m = 1 : length(RawAmpl) %Names and Amplitudes
                    NameText = [NameText, [basisSetNames{m} ' \n']];
                    RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    if isinf(CRLB(m))
                        CRLBText = [CRLBText, [num2str(round(CRLB(m),1), '%1.3g') '\n']]; 
                    else if CRLB(m) > 999
                            CRLBText = [CRLBText, [num2str(Inf, '%1.3g') '\n']]; 
                        else
                            CRLBText = [CRLBText, [num2str(round(CRLB(m),1), '%1.3g') '%%\n']];
                        end
                    end

                end
            else %Water/reference fit but this should never happen in this loop
                NameText = ['Water: ' ];
                RawAmplText = [num2str(RawAmpl,'%1.2e')];
            end
            set(gui.Results.fit{t}, 'Title', ['Raw Amplitudes']);
            gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
            gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                'FontSize', resultsFontSize, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                'FontSize', resultsFontSize, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            gui.Results.FitTextCRLB  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', resultsFontSize, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(CRLBText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
            if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                switch MRSCont.opts.fit.method
                    case 'Osprey'
                        if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results.ref.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                        else
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results.w.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                        end
                    case 'LCModel'
                end
                NameText = [''];
                RawAmplText = [''];
                CRLBText    = [''];
                for m = 1 : length(RawAmpl) %Names and Amplitudes
                    NameText = [NameText, [basisSetNames{m} ' \n']];
                    RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    if isinf(CRLB(m))
                        CRLBText = [CRLBText, [num2str(round(CRLB(m),1), '%1.3g') '\n']]; 
                    else if CRLB(m) > 999
                            CRLBText = [CRLBText, [num2str(Inf, '%1.3g') '\n']]; 
                        else
                            CRLBText = [CRLBText, [num2str(round(CRLB(m),1), '%1.3g') '%%\n']];
                        end
                    end

                end

                set(gui.Results.fit{t}, 'Title', ['Raw Water Ratio']);
                gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', resultsFontSize, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', resultsFontSize, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.Results.FitTextCRLB  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', resultsFontSize, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(CRLBText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            else %Water/reference fit
                NameText = ['Water: ' ];
                RawAmplText = [num2str(RawAmpl,'%1.2e')];
                set(gui.Results.fit{t}, 'Title', ['Raw Amplitudes']);
                gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            end
        end
    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
        if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
            if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                NameText = [''];
                RawAmplText = [''];
                for m = 1 : length(RawAmpl) %Names and Amplitudes
                    NameText = [NameText, [MRSCont.fit.resBasisSet{1,gui.controls.act_x}.(gui.fit.Style){1,MRSCont.info.A.unique_ndatapoint_indsort(gui.controls.Selected)}.name{m} ': \n']];
                    RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                end
            else %Water/reference fit but this should never happen in this loop
                NameText = ['Water: ' ];
                RawAmplText = [num2str(RawAmpl,'%1.2e')];
            end
            set(gui.Results.fit{t}, 'Title', ['Raw Amplitudes']);
            gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
            gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
            if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                    RawAmpl = RawAmpl ./ (MRSCont.fit.results{1,gui.controls.act_x}.ref.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                else
                    RawAmpl = RawAmpl ./ (MRSCont.fit.results{1,gui.controls.act_x}.w.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                end

 %%% 4. FILLING FITTED AMPLITUDE PANEL %%%
 % Creates the panel on the right side with the fitted ammplitudes
            gui.InfoText.fit{t}  = uicontrol('Parent',gui.upperBox.fit.Info,'style','text',...
                                        'FontSize', 12, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            gui.Results.fit{t} = uix.Panel('Parent', gui.Plot.fit,...
                                       'Title', ['Raw Amplitudes'],'FontName', gui.font,'HighlightColor', gui.colormap.Foreground,...
                                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);

            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                    if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                        NameText = [''];
                        RawAmplText = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fit.Style){1,MRSCont.info.A.unique_ndatapoint_indsort(gui.controls.Selected)}.name{m} ': \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                        end
                    else %Water/reference fit but this should never happen in this loop
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                    end
                    set(gui.Results.fit{t}, 'Title', ['Raw Amplitudes']);
                        gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                        gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                    if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                        if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results.ref.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                        else
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results.w.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                        end
                        NameText = [''];
                        RawAmplText = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fit.Style).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                        end
                        set(gui.Results.fit{t}, 'Title', ['Raw Water Ratio']);
                        gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                        gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    else %Water/reference fit
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                       set(gui.Results.fit{t}, 'Title', ['Raw Amplitudes']);
                       gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                       gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                       'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                       gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                       'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    end
                end
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                    if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                        NameText = [''];
                        RawAmplText = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            NameText = [NameText, [MRSCont.fit.resBasisSet{1,gui.controls.act_x}.(gui.fit.Style).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                        end
                    else %Water/reference fit but this should never happen in this loop
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                    end
                    set(gui.Results.fit{t}, 'Title', ['Raw Amplitudes']);
                        gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                        gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                    if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                        if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results{1,gui.controls.act_x}.ref.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                        else
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results{1,gui.controls.act_x}.w.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                        end
                        NameText = [''];
                        RawAmplText = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            NameText = [NameText, [MRSCont.fit.resBasisSet{1,gui.controls.act_x}.(gui.fit.Style).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                        end
                        set(gui.Results.fit{t}, 'Title', ['Raw Water Ratio']);
                        gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                        gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    else %Water/reference fit
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                       set(gui.Results.fit{t}, 'Title', ['Raw Amplitudes']);
                       gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                       gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                       'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                       gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                       'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    end
                end
            else
                if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                    if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                        NameText = [''];
                        RawAmplText = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            NameText = [NameText, [MRSCont.fit.resBasisSet{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                        end
                    else %Water/reference fit but this should never happen in this loop
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                    end
                    set(gui.Results.fit{t}, 'Title', ['Raw Amplitudes']);
                        gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                        gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                    if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                        if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.ref.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                        else
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.w.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                        end
                        NameText = [''];
                        RawAmplText = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            try
                                NameText = [NameText, [MRSCont.fit.resBasisSet{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                            catch
                                NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fit.Style){1,1}.name{m} ': \n']];
                            end
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                        end
                        set(gui.Results.fit{t}, 'Title', ['Raw Water Ratio']);
                        gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                        gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    else %Water/reference fit
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                       set(gui.Results.fit{t}, 'Title', ['Raw Amplitudes']);
                       gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                       gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                       'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                       gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                       'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    end
                end
            end
            set(gui.Results.fit{t}, 'Title', ['Raw Amplitudes']);
            gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
            gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
            if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                    RawAmpl = RawAmpl ./ (MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.ref.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                else
                    RawAmpl = RawAmpl ./ (MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.w.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                end
                NameText = [''];
                RawAmplText = [''];
                for m = 1 : length(RawAmpl) %Names and Amplitudes
                    try
                        NameText = [NameText, [MRSCont.fit.resBasisSet{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style){1,MRSCont.info.A.unique_ndatapoint_indsort(gui.controls.Selected)}.name{m} ': \n']];
                    catch
                        NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fit.Style){1,1}.name{m} ': \n']];
                    end
                    RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                end
                set(gui.Results.fit{t}, 'Title', ['Raw Water Ratio']);
                gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            else %Water/reference fit
                NameText = ['Water: ' ];
                RawAmplText = [num2str(RawAmpl,'%1.2e')];
                set(gui.Results.fit{t}, 'Title', ['Raw Amplitudes']);
                gui.Results.FitText = uix.HBox('Parent', gui.Results.fit{t}, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            end
        end
        end
    end
    %%%  5. VISUALIZATION PART OF THIS TAB %%%
    %osp_plotFit is used to visualize the fits (off,diff1,diff2,sum,ref,water)
    temp = figure( 'Visible', 'off' );
    if ~((isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI))
        temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,[gui.controls.act_x gui.controls.act_y gui.controls.act_z],gui.fit.Names{t});
    elseif MRSCont.flags.isPRIAM
        temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,gui.controls.act_x,Selection); %Create figure
    else
        temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,[gui.controls.act_x gui.controls.act_y],Selection); %Create figure
    end
    ViewAxes = gca();
    set(ViewAxes, 'Parent', gui.Plot.fit{t} );
    close( temp );

    set(gui.Plot.fit{t},'Widths', [-0.16 -0.84]);
    set(gui.Plot.fit{t}.Children(2), 'Units', 'normalized');
    set(gui.Plot.fit{t}.Children(2), 'OuterPosition', [0.17,0.02,0.75,0.98])
end
h = findall(groot,'Type','figure');
for ff = 1 : length(h)
    if ~(strcmp(h(ff).Tag, 'Osprey') ||  strcmp(h(ff).Tag, 'TMWWaitbar'))
        close(h(ff))
    end
end
setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
end
