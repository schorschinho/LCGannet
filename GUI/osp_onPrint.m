function osp_onPrint( ~, ~ ,gui)
%% osp_onPrint
%   Callback function on print figure button click.
%
%
%   USAGE:
%       osp_onPrint( ~, ~ ,gui);
%
%   INPUT:      gui      = gui class containing all handles and the MRSCont
%
%   OUTPUT:     Changes in gui parameters and MRSCont are written into the
%               gui class
%
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
%%% 1. GET DATA %%%
    MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class
    selectedTab = get(gui.layout.tabs, 'Selection');
    screenSize      = [1, 1, 1300, 1100];
    canvasSize      = screenSize;
    canvasSize(4)   = screenSize(4) * 0.7;
    canvasSize(3)   = canvasSize(4) * (11/8.5);
    canvasSize(2)   = (screenSize(4) - canvasSize(4))/2;
    canvasSize(1)   = (screenSize(3) - canvasSize(3))/2;
    out = figure('NumberTitle', 'off', 'Visible', 'on', 'Menu', 'none','Position', canvasSize,...
                    'ToolBar', 'none', 'HandleVisibility', 'off', 'Renderer', 'painters', 'Color', gui.colormap.Background);

    Title = MRSCont.ver.Osp;

    Frame = uix.Panel('Parent',out, 'Padding', 1, 'Title', Title,...
                                 'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                 'HighlightColor', gui.colormap.Background, 'ShadowColor', gui.colormap.Background);
    input_figure = uix.VBox('Parent', Frame,  'BackgroundColor',gui.colormap.Background, 'Spacing', 5);
    box = uix.HBox('Parent', input_figure,'BackgroundColor',gui.colormap.Background, 'Spacing',6);
    Info = uix.Panel('Parent',box, 'Padding', 5, 'Title', MRSCont.files{gui.controls.Selected},...
                                 'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
    LogoFig = figure('Visible','off');
    [I, map] = imread('osprey.gif','gif');
    axes(LogoFig, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
    ViewAxes = gca();
    set( ViewAxes,'Box','off','XColor', 'none','YColor','none', 'Parent', box );

    set(box, 'Width', [-0.9 -0.1]);
    switch selectedTab
        case 1 %Load tab
            outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyLoad');
            [~,filename,~]  = fileparts(MRSCont.files{gui.controls.Selected});

                % Grid for Plot and Data control sliders
                if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) % HBox for HERMES/HERCULES
                    Plot = uix.HBox('Parent', input_figure, 'BackgroundColor',gui.colormap.Background, 'Units', 'normalized');
                else
                    Plot= uix.VBox('Parent', input_figure, 'BackgroundColor',gui.colormap.Background, 'Units', 'normalized');
                end
                InfoText  = uicontrol('Parent',Info,'style','text',...
                                              'FontSize', 12, 'FontName', gui.font,'ForegroundColor', gui.colormap.Foreground,...
                                              'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
            % Get parameter from file to fill the info panel
            switch gui.load.Names.Spec{gui.load.Selected}
                case 'metabolites'
                StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) '\naverages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
                             '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                             num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
                case 'MM'
                        StatText = ['MM Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.spectralwidth) ' Hz'...   %re_mm
                             '\nraw subspecs: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...   %re_mm
                             num2str(MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];   %re_mm
                case 'reference'
                StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.spectralwidth) ' Hz'...   %re_mm
                             '\nraw subspecs: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...   %re_mm
                             num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];   %re_mm
               case 'MM reference'
                        StatText = ['MM reference Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                             num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_mm_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_mm_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
                case 'water'
                    StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                             num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
                end

            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                set(InfoText, 'String',sprintf(StatText))
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                StatText = [StatText '; Voxel ' num2str(gui.controls.act_x)];
                set(InfoText, 'String',sprintf(StatText))
            end
     %%% 4. VISUALIZATION PART OF THIS TAB %%%
     %osp_plotLoad is used to visualize the raw data. Number of subplots
     %depends on the number of subspectra of the seuqence
            Exp = gui.controls.act_Exp;
            switch gui.load.Names.Spec{gui.load.Selected}
            case 'metabolites'
                if Exp > max(max(MRSCont.opts.MultipleSpectra.metab))
                    Exp = 1;
                    gui.controls.act_Exp = 1;
                end
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mets',Exp);
                    VoxelIndex = 1;
                elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mets',gui.controls.act_x);
                    VoxelIndex = gui.controls.act_x;
                end
                outputFile      = [filename '_Voxel_' num2str(VoxelIndex) '_Exp_' num2str(Exp) '_OspreyLoad_mets.pdf'];
                if MRSCont.flags.isUnEdited % One window for UnEdited
                    drawnow
                    set( temp.Children, 'Parent', Plot );
                end
                if MRSCont.flags.isMEGA %Two windows for MEGA
                    set( temp.Children(2), 'Parent', Plot );
                    set( temp.Children(1), 'Parent', Plot );
                    set(Plot,'Heights', [-0.49 -0.49]);
                    set(Plot.Children(2), 'Units', 'normalized')
                    set(Plot.Children(2), 'OuterPosition', [0,0.5,1,0.5])
                    set(Plot.Children(1), 'Units', 'normalized')
                    set(Plot.Children(1), 'OuterPosition', [0,0,1,0.5])
                end
                if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) %Four windows for HERMES/HERCULES
                   multiACload = uix.VBox('Parent', Plot, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                        multiAload = uix.VBox('Parent', multiACload,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                        multiCload = uix.VBox('Parent', multiACload,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    multiBDload = uix.VBox('Parent', Plot,'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                        multiBload = uix.VBox('Parent', multiBDload, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                        multiDload = uix.VBox('Parent', multiBDload, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    set( temp.Children(1), 'Parent', multiDload );
                    set( temp.Children(1), 'Parent', multiCload );
                    set( temp.Children(1), 'Parent', multiBload );
                    set( temp.Children(1), 'Parent', multiAload );
                    set(Plot,'Width', [-0.49 -0.49]);
                    set(multiDload.Children(1), 'Units', 'normalized')
                    set(multiDload.Children(1), 'OuterPosition', [0,0,1,1])
                    set(multiCload.Children(1), 'Units', 'normalized')
                    set(multiCload.Children(1), 'OuterPosition', [0,0,1,1])
                    set(multiBload.Children(1), 'Units', 'normalized')
                    set(multiBload.Children(1), 'OuterPosition', [0,0,1,1])
                    set(multiAload.Children(1), 'Units', 'normalized')
                    set(multiAload.Children(1), 'OuterPosition', [0,0,1,1])

                end
                case 'MM'
                    if Exp > max(MRSCont.opts.MultipleSpectra.mm)
                        Exp = 1;
                        gui.controls.act_Exp = 1;
                    end
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm',Exp);
                         VoxelIndex = 1;
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm',Exp,gui.controls.act_x);
                         VoxelIndex = gui.controls.act_x;
                          else
                            temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm',Exp,[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                end
                outputFile      = [filename '_Voxel_' num2str(VoxelIndex) '_Exp_' num2str(Exp) '_OspreyLoad_MM.pdf'];
                if MRSCont.flags.isUnEdited % One window for UnEdited
                    drawnow
                    set( temp.Children, 'Parent', Plot );
                end
                if MRSCont.flags.isMEGA %Two windows for MEGA
                    set( temp.Children(2), 'Parent', Plot );
                    set( temp.Children(1), 'Parent', Plot );
                    set(Plot,'Heights', [-0.49 -0.49]);
                    set(Plot.Children(2), 'Units', 'normalized')
                    set(Plot.Children(2), 'OuterPosition', [0,0.5,1,0.5])
                    set(Plot.Children(1), 'Units', 'normalized')
                    set(Plot.Children(1), 'OuterPosition', [0,0,1,0.5])
                end
                if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) %Four windows for HERMES/HERCULES
                    gui.layout.multiACload = uix.VBox('Parent', Plot, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                        gui.layout.multiAload = uix.VBox('Parent', gui.layout.multiACload,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                        gui.layout.multiCload = uix.VBox('Parent', gui.layout.multiACload,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    gui.layout.multiBDload = uix.VBox('Parent', Plot,'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                        gui.layout.multiBload = uix.VBox('Parent', gui.layout.multiBDload, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                        gui.layout.multiDload = uix.VBox('Parent', gui.layout.multiBDload, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    set( temp.Children(1), 'Parent', gui.layout.multiDload );
                    set( temp.Children(1), 'Parent', gui.layout.multiCload );
                    set( temp.Children(1), 'Parent', gui.layout.multiBload );
                    set( temp.Children(1), 'Parent', gui.layout.multiAload );
                    set(Plot,'Width', [-0.49 -0.49]);
                    set(gui.layout.multiDload.Children(1), 'Units', 'normalized')
                    set(gui.layout.multiDload.Children(1), 'OuterPosition', [0,0,1,1])
                    set(gui.layout.multiCload.Children(1), 'Units', 'normalized')
                    set(gui.layout.multiCload.Children(1), 'OuterPosition', [0,0,1,1])
                    set(gui.layout.multiBload.Children(1), 'Units', 'normalized')
                    set(gui.layout.multiBload.Children(1), 'OuterPosition', [0,0,1,1])
                    set(gui.layout.multiAload.Children(1), 'Units', 'normalized')
                    set(gui.layout.multiAload.Children(1), 'OuterPosition', [0,0,1,1])

                end
            case 'reference'
                if Exp > max(MRSCont.opts.MultipleSpectra.ref)
                    Exp = 1;
                    gui.controls.act_Exp = 1;
                end
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref',Exp);
                        VoxelIndex = 1;
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref',Exp,gui.controls.act_x);
                        VoxelIndex = gui.controls.act_x;
                      else
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref',Exp,[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                end
                drawnow
                set( temp.Children, 'Parent', Plot );
                outputFile      = [filename '_Voxel_' num2str(VoxelIndex) '_Exp_' num2str(Exp) '_OspreyLoad_ref.pdf'];

             case 'MM reference'
                if Exp > max(MRSCont.opts.MultipleSpectra.mm_ref)
                    Exp = 1;
                    gui.controls.act_Exp = 1;
                end
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm_ref',Exp);
                        VoxelIndex = 1;
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm_ref',Exp,gui.controls.act_x);
                        VoxelIndex = gui.controls.act_x;
                      else
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm_ref',Exp,[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                end
                drawnow
               set( temp.Children, 'Parent', Plot );
                outputFile      = [filename '_Voxel_' num2str(VoxelIndex)  '_Exp_' num2str(Exp)  '_OspreyLoad_mm_ref.pdf'];

             case 'water'
                if Exp > max(MRSCont.opts.MultipleSpectra.w)
                    Exp = 1;
                    gui.controls.act_Exp = 1;
                end
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',Exp);
                    VoxelIndex = 1;
                elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',Exp,gui.controls.act_x);
                    VoxelIndex = gui.controls.act_x;
              else
                temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',Exp,[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                end
                drawnow
                set( temp.Children, 'Parent', Plot );
                outputFile      = [filename '_Voxel_' num2str(VoxelIndex)  '_Exp_' num2str(Exp) '_OspreyLoad_w.pdf'];
            end
            set(input_figure, 'Heights', [-0.1 -0.9]);
            % Get rid of the Load figure
            close( temp );

        case 2 %Process Tab
            outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyProcess');
            [~,filename,~]  = fileparts(MRSCont.files{gui.controls.Selected});
            Selection = gui.process.TabTitles{gui.process.Selected};
            Exp = gui.controls.act_x;
            SubSpec = gui.controls.act_y;
            Plot = uix.HBox('Parent', input_figure, ...
                'Padding', 5,'BackgroundColor', gui.colormap.Background);
            set(input_figure, 'Heights', [-0.11 -0.89]);

            % Get parameter from file to fill the info panel
            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                if strcmp(gui.process.TabTitles{gui.process.Selected},'metab')
                     StatText = ['SNR(' gui.process.(gui.process.TabTitles{gui.process.Selected}).SNR{SubSpec} '): '  num2str(MRSCont.QM.SNR.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected,SubSpec)) '; FWHM (' gui.process.(gui.process.TabTitles{gui.process.Selected}).SNR{SubSpec} '): '...
                                num2str(MRSCont.QM.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected,SubSpec)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected,SubSpec)/MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM.freqShift.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected,SubSpec)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM.drift.pre.AvgDeltaCr.A(gui.controls.Selected)*MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM.drift.post.AvgDeltaCr.A(gui.controls.Selected)*MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)/1e6) ' Hz'];
                else
                    StatText = ['SNR(' gui.process.(gui.process.TabTitles{gui.process.Selected}).SNR{1} '): ' num2str(MRSCont.QM.SNR.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) '; FWHM (' gui.process.(gui.process.TabTitles{gui.process.Selected}).SNR{1} '): '...
                            num2str(MRSCont.QM.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)/MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)*1e6))...
                            ' Hz / ppm'];
                end
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                if (strcmp(gui.process.Names{gui.process.Selected},'A') || strcmp(gui.process.Names{gui.process.Selected},'B') || strcmp(gui.process.Names{gui.process.Selected},'C') || strcmp(gui.process.Names{gui.process.Selected},'D') || strcmp(gui.process.Names{gui.process.Selected},'diff1') || strcmp(gui.process.Names{gui.process.Selected},'diff2') || strcmp(gui.process.Names{gui.process.Selected},'sum') )
                    StatText = ['Metabolite Data -> SNR(' gui.process.SNR{gui.process.Selected} '): '  num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                                num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM{1,gui.controls.act_x}.freqShift.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM{1,gui.controls.act_x}.drift.pre.AvgDeltaCr.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM{1,gui.controls.act_x}.drift.post.AvgDeltaCr.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq/1e6) ' Hz'];
                else if (strcmp(gui.process.Names{gui.process.Selected},'ref') || strcmp(gui.process.Names{gui.process.Selected},'mm'))
                StatText = ['Reference Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                            num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                            ' Hz / ppm'];
                    else
                        StatText = ['Water Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                                    num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '/' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                                    ' Hz / ppm'];
                    end
                end
            end

            if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                StatText = ['Voxel ' num2str(gui.controls.act_x) ': ' StatText];
            end
            InfoText  = uicontrol('Parent',Info,'style','text','FontSize', 12, 'FontName', gui.font,...
                                         'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                         'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

 %%% 4. VISUALIZATION PART OF THIS TAB %%%
 %osp_plotProcess is used to visualize the processed spectra
            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                temp = osp_plotProcess(MRSCont, gui.controls.Selected,Selection,SubSpec,Exp); % Create figure
                VoxelIndex = 1;
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                temp = osp_plotProcess(MRSCont, gui.controls.Selected,Selection,SubSpec,Exp,gui.controls.act_x); %Create figure
                VoxelIndex = gui.controls.act_x;
            end
            %Subplots are distributed here
                proSpecs = uix.VBox('Parent', Plot, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    proPre = uix.VBox('Parent', proSpecs,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    proPost = uix.VBox('Parent', proSpecs,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                proOut = uix.VBox('Parent', Plot,'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    proDrift = uix.VBox('Parent', proOut, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    proAlgn = uix.VBox('Parent', proOut, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);

            set( temp.Children(1), 'Parent', proDrift );
            set( temp.Children(1), 'Parent', proAlgn );
            set( temp.Children(1), 'Parent', proPost );
            set( temp.Children(1), 'Parent', proPre );
            close( temp );
%%% 5. DATA CONTROLS FOR THIS TAB %%%
            set(Plot,'Widths', [-0.49 -0.49]);
            set(proPre.Children(1), 'Units', 'normalized')
            set(proPre.Children(1), 'OuterPosition', [0,0,1,1])
            set(proPost.Children(1), 'Units', 'normalized')
            set(proPost.Children(1), 'OuterPosition', [0,0,1,1])
            set(proDrift.Children(1), 'Units', 'normalized')
            set(proDrift.Children(1), 'OuterPosition', [0,0,1,1])
            set(proAlgn.Children(1), 'Units', 'normalized')
            set(proAlgn.Children(1), 'OuterPosition', [0,0,1,1])

            outputFile      = [filename '_Voxel_' num2str(VoxelIndex) '_OspreyProcess_' Selection '_' MRSCont.processed.(Selection){gui.controls.Selected}.names{SubSpec} '.pdf'];
        case 3 %Fit
             outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyFit');
            % For this visualization, we will have to make a few
            % distinctions upfront since the modeling algorithms (LCModel
            % vs. Osprey) do not always return the same kinds of data, or they
            % return them in different formats.

            basis = gui.controls.act_z;
            subspectrum = gui.controls.act_y;
            switch MRSCont.opts.fit.method
                case 'LCModel'
                    % Number of metabolites and lipid/MM basis functions
                    basisNames = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,subspectrum}.name;
                    nLip    = sum(~cellfun(@isempty, strfind(basisNames, 'Lip')));
                    nMM     = sum(~cellfun(@isempty, strfind(basisNames, 'MM')));
                    nMMLip  = nLip + nMM;
                    nMets   = length(basisNames) - nMMLip;
                    nComb   = sum(~cellfun(@isempty, strfind(basisNames, '_')));
                    % No info panel string for the water fit range
                    waterFitRangeString = '';
                    % Where are the metabolite names stored?
                    basisSetNames = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,subspectrum}.name;
                    if subspectrum == 1                 %Works only or MEGA-PRESS and unedited
                        subSpecName = 'A';
                    else
                        subSpecName = 'diff1';
                    end
                    % Smaller fonts for the results
                    resultsFontSize = 6;
                case 'Osprey'

                    % Additional info panel string for the water fit range
                    waterFitRangeString = ['Fitting range: ' num2str(MRSCont.opts.fit.rangeWater(1)) ' to ' num2str(MRSCont.opts.fit.rangeWater(2)) ' ppm'];
                    % Where are the metabolite names stored?
                    if strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')
                        basisSet = MRSCont.fit.resBasisSet.(gui.fit.Style).(['np_sw_' num2str(round(MRSCont.processed.metab{gui.controls.Selected}.sz(1))) '_' num2str(round(MRSCont.processed.metab{gui.controls.Selected}.spectralwidth))]){1};
                        basisSetNames = basisSet.name;
                        subSpecName = gui.fit.Style;
                    else if strcmp(gui.fit.Style, 'conc')
                            basisSet = MRSCont.fit.resBasisSet.(gui.fit.Style).(['np_sw_' num2str(round(MRSCont.processed.metab{gui.controls.Selected}.sz(1))) '_' num2str(round(MRSCont.processed.metab{gui.controls.Selected}.spectralwidth))]){basis,1};
                            basisSetNames = basisSet.name;
                            subSpecName = basisSet.names{1};
                            else
                            basisSet = MRSCont.fit.resBasisSet.(gui.fit.Style).(['np_sw_' num2str(round(MRSCont.processed.metab{gui.controls.Selected}.sz(1))) '_' num2str(round(MRSCont.processed.metab{gui.controls.Selected}.spectralwidth))]){basis,1,subspectrum};
                            basisSetNames = basisSet.name;
                            subSpecName = basisSet.names{1};
                        end
                    end
                % Number of metabolites and lipid/MM basis functions
                nMets   = MRSCont.fit.basisSet.nMets;
                nMMLip  = MRSCont.fit.basisSet.nMM;
                 % Larger fonts for the results
                resultsFontSize = 11;
            end
            [~,filename,~]  = fileparts(MRSCont.files{gui.controls.Selected});
            Selection = gui.fit.Names{gui.fit.Selected};
            Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
            set(input_figure, 'Heights', [-0.12 -0.88]);
            if  ~strcmp (MRSCont.opts.fit.style, 'Concatenated') ||  strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') || strcmp(gui.fit.Names{gui.fit.Selected}, 'w') %Is not concateneted or is reference/water fit
                gui.fit.Style = gui.fit.Names{gui.fit.Selected};
                switch MRSCont.opts.fit.method
                    case 'LCModel'
                        if strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') || strcmp(gui.fit.Names{gui.fit.Selected}, 'w')
                            RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,subspectrum}.h2oarea .* MRSCont.fit.scale{gui.controls.Selected};
                        else
                            RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,subspectrum}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
                            CRLB    = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,subspectrum}.CRLB;
                        end
                    case 'Osprey'
                        RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
                end
            else %Is concatenated and not water/reference
                gui.fit.Style = 'conc';
            end
            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                ph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.ph0;
                ph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.ph1;
                if ~strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') && ~strcmp(gui.fit.Names{gui.fit.Selected}, 'w')
                    refShift = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.refShift;
                    refFWHM = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.refFWHM;
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        iniph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph0;
                        iniph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph1;
                    case 'LCModel'
                        iniph0 = nan;
                        iniph1 = nan;
                    end
                end
                RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                if ~strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') && ~strcmp(gui.fit.Names{gui.fit.Selected}, 'w')
                    refShift = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refShift;
                    refFWHM = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refFWHM;
                    ph0 = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph0;
                    ph1 = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph1;
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        iniph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph0;
                        iniph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph1;
                    case 'LCModel'
                        iniph0 = nan;
                        iniph1 = nan;
                    end
                end
                RawAmpl = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
            else
               if ~strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') && ~strcmp(gui.fit.Names{gui.fit.Selected}, 'w')
                    refShift = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refShift;
                    refFWHM = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refFWHM;
                    ph0 = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph0;
                    ph1 = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph1;
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        iniph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph0;
                        iniph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph1;
                    case 'LCModel'
                        iniph0 = nan;
                        iniph1 = nan;
                    end
                end
                RawAmpl = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
            end
            % Get parameter from file to fill the info panel
            if  ~strcmp (Selection, 'ref') && ~strcmp (Selection, 'w') %Metabolite data
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        iniph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.prelimParams.ph0;
                        iniph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.prelimParams.ph1;
                    case 'LCModel'
                        iniph0 = nan;
                        iniph1 = nan;
                    end
                    StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' gui.fit.Names{gui.fit.Selected},...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.range(1)) ' to ' num2str(MRSCont.opts.fit.range(2)) ' ppm; Baseline knot spacing: ' num2str(MRSCont.opts.fit.bLineKnotSpace) ' ppm; ph0: ' num2str(ph0,'%1.2f'),...
                        'deg; ph1: ' num2str(ph1,'%1.2f') 'deg; refShift: ' num2str(refShift,'%1.2f') ' Hz; refFWHM: ' num2str(refFWHM,'%1.2f')...
                        ' ppm\nNumber of metabolites: ' num2str(nMets) '; Number of MM/lipids: ' num2str(nMMLip) ...
                        ' scale: '  num2str(MRSCont.fit.scale{gui.controls.Selected})];
                elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        iniph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.prelimParams.ph0;
                        iniph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.prelimParams.ph1;
                    case 'LCModel'
                        iniph0 = nan;
                        iniph1 = nan;
                    end
                    StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' gui.fit.Names{gui.fit.Selected},...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.range(1)) ' to ' num2str(MRSCont.opts.fit.range(2)) ' ppm; Baseline knot spacing: ' num2str(MRSCont.opts.fit.bLineKnotSpace) ' ppm; ph0: ' num2str(ph0,'%1.2f'),...
                        'deg; ph1: ' num2str(ph1,'%1.2f') 'deg; refShift: ' num2str(refShift,'%1.2f') ' Hz; refFWHM: ' num2str(refFWHM,'%1.2f')...
                        ' ppm\nNumber of metabolites: ' num2str(nMets) '; Number of MM/lipids: ' num2str(nMMLip) ...
                        ' scale: '  num2str(MRSCont.fit.scale{gui.controls.Selected})];
                else
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        iniph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.prelimParams.ph0;
                        iniph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.prelimParams.ph1;
                    case 'LCModel'
                        iniph0 = nan;
                        iniph1 = nan;
                    end
                    StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' gui.fit.Names{gui.fit.Selected},...
                            '\nFitting range: ' num2str(MRSCont.opts.fit.range(1)) ' to ' num2str(MRSCont.opts.fit.range(2)) ' ppm; Baseline knot spacing: ' num2str(MRSCont.opts.fit.bLineKnotSpace) ' ppm; ph0: ' num2str(ph0,'%1.2f'),...
                            'deg; ph1: ' num2str(ph1,'%1.2f') 'deg; refShift: ' num2str(refShift,'%1.2f') ' Hz; refFWHM: ' num2str(refFWHM,'%1.2f')...
                            ' ppm\nNumber of metabolites: ' num2str(nMets) '; Number of MM/lipids: ' num2str(nMMLip) ...
                            ' scale: '  num2str(MRSCont.fit.scale{gui.controls.Selected})];
                end

            else if strcmp (Selection, 'ref') %Reference data?
            StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' Selection,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.rangeWater(1)) ' to ' num2str(MRSCont.opts.fit.rangeWater(2)) ' ppm'];
                else %Is water data
                    StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' Selection,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.rangeWater(1)) ' to ' num2str(MRSCont.opts.fit.rangeWater(2)) ' ppm'];
                end
            end
 %%% 4. FILLING FITTED AMPLITUDE PANEL %%%
 % Creates the panel on the right side with the fitted ammplitudes
            InfoText  = uicontrol('Parent',Info,'style','text',...
                                        'FontSize', 12, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            Results = uix.Panel('Parent', Plot,...
                                       'Title', ['Raw Amplitudes'],'FontName', gui.font,'HighlightColor', gui.colormap.Foreground,...
                                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                    if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                        NameText = [''];
                        RawAmplText = [''];
                        CRLBText    = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            NameText = [NameText, [basisSetNames{m} ' \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                            if strcmp(MRSCont.opts.fit.method, 'LCModel')
                                CRLBText = [CRLBText, [num2str(CRLB(m), '%i') '%%\n']];
                            end
                        end
                    else %Water/reference fit but this should never happen in this loop
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                    end
                    set(Results, 'Title', ['Raw Amplitudes']);
                        FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                        FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        if strcmp(MRSCont.opts.fit.method, 'LCModel')
                            FitTextCRLB  = uicontrol('Parent',FitText,'style','text',...
                            'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(CRLBText),...
                            'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        end
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
                            if strcmp(MRSCont.opts.fit.method, 'LCModel')
                                CRLBText = [CRLBText, [num2str(CRLB(m), '%i') '%%\n']];
                            end
                        end
                        set(Results, 'Title', ['Raw Water Ratio']);
                        FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                        FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        if strcmp(MRSCont.opts.fit.method, 'LCModel')
                            FitTextCRLB  = uicontrol('Parent',FitText,'style','text',...
                            'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(CRLBText),...
                            'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        end
                    else %Water/reference fit
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                       set(Results, 'Title', ['Raw Amplitudes']);
                       FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                       FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                       'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                       FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
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
                            NameText = [NameText, [MRSCont.fit.resBasisSet{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                        end
                    else %Water/reference fit but this should never happen in this loop
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                    end
                    set(Results, 'Title', ['Raw Amplitudes']);
                        FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                        FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
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
                        set(Results, 'Title', ['Raw Water Ratio']);
                        FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                        FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    else %Water/reference fit
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                       set(Results, 'Title', ['Raw Amplitudes']);
                       FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                       FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                       'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                       FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
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
                    set(Results, 'Title', ['Raw Amplitudes']);
                        FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                        FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
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
                        set(Results, 'Title', ['Raw Water Ratio']);
                        FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                        FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                        FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    else %Water/reference fit
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                       set(Results, 'Title', ['Raw Amplitudes']);
                       FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                       FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                       'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                       FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                       'FontSize', 11, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    end
                end
            end

%%%  5. VISUALIZATION PART OF THIS TAB %%%
%osp_plotFit is used to visualize the fits (off,diff1,diff2,sum,ref,water)
            temp = figure( 'Visible', 'off' );
            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,[gui.controls.act_x gui.controls.act_y gui.controls.act_z],Selection);
                VoxelIndex = 1;
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,gui.controls.act_x,Selection);
                VoxelIndex = gui.controls.act_x;
            else
                temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,[gui.controls.act_x,gui.controls.act_y],Selection);
                VoxelIndex = gui.controls.act_x;
            end
            ViewAxes = gca();
            set(ViewAxes, 'Parent', Plot );
            close( temp );

            set(Plot,'Widths', [-0.16 -0.84]);
            set(Plot.Children(2), 'Units', 'normalized');
            set(Plot.Children(2), 'OuterPosition', [0.17,0.02,0.75,0.98])
            outputFile      = [filename '_Voxel_' num2str(VoxelIndex) '_OspreyFit_' gui.fit.Style '_' subSpecName '_basis_' num2str(gui.controls.act_z) '.pdf'];
        case 4 %Coreg/Seg
            outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyCoregSeg');
            [~,filename,~]  = fileparts(MRSCont.files{gui.controls.Selected});

            % Creates layout for plotting and data control
            Plot = uix.HBox('Parent', input_figure,'BackgroundColor',gui.colormap.Background);
            set(input_figure, 'Heights', [-0.1 -0.9]);
            % Get parameter from file to fill the info panel

            StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '\nraw subspecs: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];

           InfoText  = uicontrol('Parent',Info,'style','text',...
                'FontSize', 12, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
            'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

        %%% 2.VISUALIZATION PART OF THIS TAB %%%
        % In this case osp_plotCoreg or osp_plotSegment is used to visualize the
        % coregistration or the segmentation
            Results = uix.VBox('Parent', Plot,'BackgroundColor',gui.colormap.Background);
            temp = figure( 'Visible', 'off' );
            warning('off','MATLAB:handle_graphics:exceptions:SceneNode');
            if MRSCont.flags.didSeg %Did segment. In this case coreg has already been performed. Visualize both
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                    temp = osp_plotCoreg(MRSCont, gui.controls.Selected);
                else
                    temp = osp_plotCoreg(MRSCont, gui.controls.Selected,gui.controls.act_x);
                end
                ViewAxes = gca();
                set(ViewAxes, 'Parent', Results );
                colormap(Results.Children,'gray')
                close( temp );
                temp = figure( 'Visible', 'off' );
                if ~isfield(MRSCont.flags,'isPRIAM')  && ~MRSCont.flags.isPRIAM
                    temp = osp_plotSegment(MRSCont, gui.controls.Selected);
                else
                    temp = osp_plotSegment(MRSCont, gui.controls.Selected,gui.controls.act_x);
                end
                ViewAxes = gca();
                set(ViewAxes, 'Parent', Results );
                colormap(Results.Children(1),'gray');
                close( temp );
            else % Only coreg has been run
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                    temp = osp_plotCoreg(MRSCont, gui.controls.Selected);
                else
                    temp = osp_plotCoreg(MRSCont, gui.controls.Selected,gui.controls.act_x);
                end
                ViewAxes = gca();
                set(ViewAxes, 'Parent', Results );
                colormap(Results.Children,'gray');
                close( temp );
            end
            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                VoxelIndex = 1;
            else
                VoxelIndex = gui.controls.act_x;
            end
            warning('on','MATLAB:handle_graphics:exceptions:SceneNode');
            outputFile      = [filename '_Voxel_' num2str(VoxelIndex) '_OspreyCoregSeg.pdf'];
        case 6 %Overview
            ovSelection = get(gui.layout.overviewTab, 'Selection');
            set(Info,'Title', 'Descriptive Information');
            groupString = '';
            for g = 1 : MRSCont.overview.NoGroups
                groupString = [groupString MRSCont.overview.groupNames{g} ' with ' num2str(sum(MRSCont.overview.groups(:) == g)) ' subjects; '];
            end
            StatText = ['Sequence: ' gui.load.Names.Seq '; Number of subjects: ' num2str(MRSCont.nDatasets) '; Number of Groups: ' num2str(MRSCont.overview.NoGroups) '\n'...
                         'Distribution: ' groupString '\n'];

           InfoText  = uicontrol('Parent',Info,'style','text',...
                'FontSize', 12, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
            'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
            set(input_figure, 'Heights', [-0.1 -0.9]);
            switch ovSelection
                case 1 %SpecOverview
                    Selection = gui.controls.pop_specsOvPlot.String(gui.overview.Selected.Spec);
                    outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview','Individual');
                    if gui.controls.GM == 0
                        outputFile  = [Selection{1} '.pdf'];
                        for g = 1 :  gui.overview.Number.Groups %Loop over groups
                            temp = osp_plotOverviewSpec(MRSCont, Selection{1},g, gui.layout.shiftind,'Frequency (ppm)','','',gui.controls.act_z);
                            if g == 1
                                fig_hold = temp;
                                ViewAxes=get(temp,'Children');
                                drawnow

                            else
                                ViewAxes=get(temp,'Children');
                                drawnow
                                copyobj(ViewAxes.Children, fig_hold.Children(1));
                                close(temp);
                            end
                        end
                        set(fig_hold.Children, 'Parent', Plot );
                    else
                        outputFile  = [Selection{1} ' basis_' num2str(gui.controls.act_z) ' Grand_mean.pdf'];
                       fig_hold = osp_plotOverviewSpec(MRSCont, Selection{1},'GMean', gui.layout.shiftind,'Frequency (ppm)','','',gui.controls.act_z);
                        ViewAxes=get(fig_hold,'Children');
                        drawnow
                       set(fig_hold.Children, 'Parent', Plot );
                    end
                    close(fig_hold);
                case 2 %MeanOverview
                    Selection = gui.controls.pop_meanOvPlot.String(gui.overview.Selected.Spec);
                    outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview', 'Mean');
                    if gui.controls.GM == 0
                        outputFile  = [Selection{1} '.pdf'];
                        for g = 1 :  gui.overview.Number.Groups
                            if gui.overview.Number.Groups > 1
                                temp = osp_plotMeanSpec(MRSCont, Selection{1},g,1,1/gui.overview.Number.Groups);
                                if g == 1
                                    fig_hold = temp;
                                else
                                    ax=get(temp,'Children');
                                    lines=get(ax,'Children');
                                    copyobj(lines, fig_hold.Children(1));
                                    close(temp);
                                end
                            else
                                fig_hold = osp_plotMeanSpec(MRSCont, Selection{1},g);
                            end
                        end
                    else
                        if ~(isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM)
                            outputFile  = [Selection{1} ' basis_' num2str(gui.controls.act_z)  ' Grand_mean.pdf'];
                            fig_hold = osp_plotMeanSpec(MRSCont, Selection{1},'GMean', 1,0,'Frequency (ppm)','','',gui.controls.act_z);
                        else
                            outputFile  = [Selection{1} ' basis_' num2str(gui.controls.act_z) 'Grand_mean_Dual_Voxel.pdf'];
                            fig_hold = osp_plotMeanSpec(MRSCont, Selection{1},1, 0.01,10,'Frequency (ppm)','','',gui.controls.act_z);
                        end
                    end
                    drawnow
                    set(fig_hold.Children, 'Parent', Plot );
                    set(out.Children.Children(1).Children(1).Children,'Children',flipud(out.Children.Children(1).Children(1).Children.Children));
                    close(fig_hold);

                case 4 %Raincloud plot
                    outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview', 'Raincloud');
                    Selection = gui.quant.popMenuNames{gui.quant.Selected.Quant};
                    if ~strcmp(Selection,'Quality')
                        split_Selection = strsplit(Selection,'-');
                        ind = find(strcmp(MRSCont.overview.FitSpecNamesStruct.(split_Selection{1})(1,:),split_Selection{2}));
                        if strcmp(split_Selection{3},'AlphaCorrWaterScaled') || strcmp(split_Selection{3},'AlphaCorrWaterScaledGroupNormed')
                            metab = 'GABA';
                        else
                            metab = MRSCont.quantify.names.(split_Selection{1}){gui.controls.act_z,ind}{gui.overview.Selected.Metab};
                        end
                        if  ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                            if ~gui.controls.GM
                                fig_hold = osp_plotRaincloud(MRSCont,split_Selection{2},split_Selection{3},metab,'Raincloud plot',0,gui.controls.act_x,gui.controls.act_z);
                            else
                                fig_hold = osp_plotRaincloud(MRSCont,split_Selection{2},split_Selection{3},metab,'Raincloud plot',1,gui.controls.act_x,gui.controls.act_z);
                            end
                            VoxelIndex = gui.controls.act_x;
                        else
                            if ~gui.controls.GM
                                fig_hold = osp_plotRaincloud(MRSCont,split_Selection{2},split_Selection{3},metab,'Raincloud plot',0,1,gui.controls.act_z);
                            else
                                fig_hold = osp_plotRaincloud(MRSCont,split_Selection{2},split_Selection{3},metab,'Raincloud plot',1,1,gui.controls.act_z);
                            end
                            VoxelIndex = 1;
                        end
                        split_Selection{4}=['basis ' num2str(gui.controls.act_z)];
                    else
                       quality = {'SNR','FWHM','freqShift'};
                       if ~gui.controls.GM
                            fig_hold = osp_plotRaincloud(MRSCont,'Quality','Quality',quality{gui.overview.Selected.Metab},'Raincloud plot');
                       else
                            fig_hold = osp_plotRaincloud(MRSCont,'Quality','Quality',quality{gui.overview.Selected.Metab},'Raincloud plot',1);
                       end
                       split_Selection{2}='Quality';
                        split_Selection{1}='Spectral';
                        split_Selection{3}=quality{gui.overview.Selected.Metab};
                        split_Selection{4}='';
                        metab = '';
                        VoxelIndex = 1;
                    end
                    delete( fig_hold.Children(1));
                    set( fig_hold.Children, 'Parent', Plot );
                    set(out.Children.Children(1).Children(1).Children,'Children',flipud(out.Children.Children(1).Children(1).Children.Children));
                    close(fig_hold);
                    if ~gui.controls.GM
                        outputFile  = ['Voxel_' num2str(VoxelIndex) '_' metab '_' split_Selection{1} '_' split_Selection{2} '_' split_Selection{3}  '_' split_Selection{4} '.pdf'];
                    else
                        outputFile  = ['Voxel_' num2str(VoxelIndex) '_' metab '_' split_Selection{1} '_' split_Selection{2} '_' split_Selection{3}  '_' split_Selection{4} '_Grand_mean.pdf'];
                    end
                case 5 %Correlation plot
                    outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview', 'Correlation');
                    Selection = gui.quant.popMenuNames{gui.quant.Selected.Quant};
                    split_Selection = strsplit(Selection,'-');
                    ind = find(strcmp(MRSCont.overview.FitSpecNamesStruct.(split_Selection{1})(1,:),split_Selection{2}));
                    MRSCont.flags.isGUI =0;
                    if strcmp(split_Selection{3},'AlphaCorrWaterScaled') || strcmp(split_Selection{3},'AlphaCorrWaterScaledGroupNormed')
                        metab = 'GABA';
                    else
                        metab = MRSCont.quantify.names.(split_Selection{1}){gui.controls.act_z,ind}{gui.overview.Selected.Metab};
                    end
                    if  ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        VoxelIndex = 1;
                        if gui.overview.Selected.CorrChoice == 1
                            switch gui.overview.Selected.Corr
                                case 1
                                fig_hold = osp_plotScatter(MRSCont,split_Selection{2},split_Selection{3},metab,MRSCont.QM.SNR.metab(1,:,ind)',gui.overview.Names.QM{gui.overview.Selected.Corr},1,gui.controls.act_z);
                                outputFile  = ['Voxel_' num2str(VoxelIndex) '_' metab '_' split_Selection{2} '_' split_Selection{3} '_basis' num2str(gui.controls.act_z) '_'  gui.overview.Names.QM{gui.overview.Selected.Corr} '.pdf'];
                                case 2
                                fig_hold = osp_plotScatter(MRSCont,split_Selection{2},split_Selection{3},metab,MRSCont.QM.FWHM.metab(1,:,ind)',gui.overview.Names.QM{gui.overview.Selected.Corr},1,gui.controls.act_z);
                                outputFile  = ['Voxel_' num2str(VoxelIndex) '_' metab '_' split_Selection{3} '_' split_Selection{3} '_basis' num2str(gui.controls.act_z) '_'  gui.overview.Names.QM{gui.overview.Selected.Corr} '.pdf'];
                            end
                        else if gui.overview.Selected.CorrChoice == 2
                            fig_hold = osp_plotScatter(MRSCont,split_Selection{2},split_Selection{3},metab,MRSCont.quantify.names.(split_Selection{1}){1,ind}{gui.overview.Selected.Corr},MRSCont.quantify.names.(split_Selection{1}){1,ind}{gui.overview.Selected.Corr},1,gui.controls.act_z);
                            outputFile  = ['Voxel_' num2str(VoxelIndex) '_' metab '_' split_Selection{2} '_' split_Selection{3} '_basis' num2str(gui.controls.act_z) '_'  MRSCont.quantify.names.(split_Selection{1}){1,ind}{gui.overview.Selected.Metab} '.pdf'];
                            else
                                fig_hold = osp_plotScatter(MRSCont,split_Selection{2},split_Selection{3},metab,gui.overview.CorrMeas{gui.overview.Selected.Corr},gui.overview.Names.Corr{gui.overview.Selected.Corr},1,gui.controls.act_z);
                                outputFile  = ['Voxel_' num2str(VoxelIndex) '_' metab '_' split_Selection{2} '_' split_Selection{3} '_basis' num2str(gui.controls.act_z) '_'  gui.overview.Names.Corr{gui.overview.Selected.Corr} '.pdf'];
                            end
                        end
                    else
                        VoxelIndex = gui.controls.act_x;
                        if gui.overview.Selected.CorrChoice == 1
                            switch gui.overview.Selected.Corr
                                case 1
                                fig_hold = osp_plotScatter(MRSCont,split_Selection{1},split_Selection{2},metab,MRSCont.QM.SNR.A',gui.overview.Names.QM{gui.overview.Selected.Corr},gui.controls.act_x);
                                outputFile  = ['Voxel_' num2str(VoxelIndex) '_' metab '_' split_Selection{1} '_' split_Selection{2} '_'  gui.overview.Names.QM{gui.overview.Selected.Corr} '.pdf'];
                                case 2
                                fig_hold = osp_plotScatter(MRSCont,split_Selection{1},split_Selection{2},metab,MRSCont.QM.FWHM.A',gui.overview.Names.QM{gui.overview.Selected.Corr},gui.controls.act_x);
                                outputFile  = ['Voxel_' num2str(VoxelIndex) '_' metab '_' split_Selection{1} '_' split_Selection{2} '_'  gui.overview.Names.QM{gui.overview.Selected.Corr} '.pdf'];
                            end
                        else if gui.overview.Selected.CorrChoice == 2
                            fig_hold = osp_plotScatter(MRSCont,split_Selection{1},gui.quant.Names.Quants{gui.quant.Selected.Quant},MRSCont.quantify.metabs.(split_Selection{1}){gui.overview.Selected.Metab},metab,metab,gui.controls.act_x);
                            outputFile  = ['Voxel_' num2str(VoxelIndex) '_' metab '_' split_Selection{1} '_' split_Selection{2} '_'  MRSCont.quantify.metabs.(split_Selection{1}){gui.overview.Selected.Metab} '.pdf'];
                            else
                               fig_hold = osp_plotScatter(MRSCont,split_Selection{1},split_Selection{2},metab,gui.overview.CorrMeas{gui.overview.Selected.Corr},gui.overview.Names.Corr{gui.overview.Selected.Corr},gui.controls.act_x);
                                outputFile  = ['Voxel_' num2str(VoxelIndex) '_' metab '_' split_Selection{1} '_' split_Selection{2} '_'  gui.overview.Names.Corr{gui.overview.Selected.Corr} '.pdf'];
                            end
                        end
                    end
                    delete( fig_hold.Children(1));
                    delete( fig_hold.Children(1));
                    set(fig_hold.Children, 'Parent', Plot );
                    set(out.Children.Children.Children(1).Children,'Children',flipud(out.Children.Children.Children(1).Children.Children));
                    close(fig_hold);
                    MRSCont.flags.isGUI =1;
                    end

    end
    set(out,'Renderer','painters','Menu','none','Toolbar','none');
    setappdata(gui.figure,'MRSCont',MRSCont);   % Write MRSCont into hidden container in gui class
%% Clean up and save

% Save the figure to the output folder
% Determine output folder
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
fig_pos = out.PaperPosition;
out.PaperSize = [fig_pos(3) fig_pos(4)];

saveas(out,fullfile(outputFolder,outputFile),'pdf');
close(out);
end
